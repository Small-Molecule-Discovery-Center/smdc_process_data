"""Insert mixtures from Excel/CSV source files and optionally rewrite the source file.

Supported source formats
------------------------
1. Excel workbooks containing component rows. The default component ID columns are
   ``SMDC`` and ``LOT_ID``.
2. CSV files like ``new_source.csv`` where component IDs are in ``ICL_IC_ID`` and
   ``lot_id``.

For either format, rows with a non-empty ``mixture_part_of`` value are grouped into
mixtures. Existing mixtures are looked up in ``INV_COMPOUND`` and
``INV_COMPOUND_MIXTURE``; those tables are not modified. A new row is inserted into
``INV_COMP_LOT`` with the next ``ICL_LOT_NUM``, followed by new rows in
``INV_COMPOUND_LOT_MIXTURE``. When an output path is supplied, the returned source
contains all non-mixture rows plus one newly coded mixture row per mixture; the original
``mixture_part`` rows are omitted.

The mixture output row is based on the first component row in the group, with these
fields replaced from the newly inserted database records:
- ``lot_id`` -> new mixture ICL/lot ID
- ``ICL_IC_ID`` -> new mixture IC/compound ID
- ``alias`` -> ``mixture_part_of``
- ``ICL_LOT_NUM`` -> inserted mixture lot number
- ``compoundindex`` -> blank
- ``mixture_part_of`` -> blank
- ``mixture_part_count`` -> blank

Usage
-----
CSV source:
  python insert_mixtures_dynamic_fields.py \
      --source new_source.csv \
      --output new_source_with_mixtures.csv \
      --host localhost \
      --db QB3_HiTS_nightly

Excel source:
  python insert_mixtures_dynamic_fields.py \
      --source mixture1.xlsx \
      --output mixture1_with_mixtures.csv \
      --host localhost \
      --db QB3_HiTS_nightly

Python:
  from insert_mixtures_dynamic_fields import insert_mixtures
  result_df = insert_mixtures(
      "new_source.csv",
      output_path="new_source_with_mixtures.csv",
      host="localhost",
      db="QB3_HiTS_nightly",
  )
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import openpyxl
import pandas as pd
from sqlalchemy import MetaData, Table, create_engine, insert, select, update
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session

sys.path.insert(0, str(Path(__file__).parent))
from mysql_cnx import create_mysql_cnx


MIXTURE_COMPOUND_CLASS = "hits.compound.Mix"
DEFAULT_MIXTURE_LOT_CLASS = "hits.lot.Mixture"
DEFAULT_SCREENING_CENTER = "SMDC"
DEFAULT_STATUS_ID = 1


@dataclass(frozen=True)
class SourceRow:
    raw: Dict[str, Any]

    def get(self, key: str, default: Any = None) -> Any:
        value = self.raw.get(key, default)
        if value == "":
            return default
        return value


def _clean(v: Any) -> Any:
    """Normalize empty/NA values so pandas NaN or literal CSV "nan" is not treated as data."""
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    if isinstance(v, str):
        stripped = v.strip()
        if not stripped or stripped.lower() in {"nan", "none", "null", "<na>"}:
            return None
        return stripped
    return v


def _to_int(v: Any) -> Optional[int]:
    v = _clean(v)
    if v is None:
        return None
    if isinstance(v, bool):
        return int(v)
    if isinstance(v, int):
        return v
    if isinstance(v, float):
        return int(v)
    return int(str(v).strip())


def _to_float(v: Any) -> Optional[float]:
    v = _clean(v)
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    return float(str(v).strip())


def _to_date(v: Any):
    v = _clean(v)
    if v is None:
        return None
    if hasattr(v, "date") and not isinstance(v, str):
        return v.date() if hasattr(v, "hour") else v
    if isinstance(v, str):
        return datetime.fromisoformat(v.strip()).date()
    return v


def _to_datetime(v: Any) -> datetime:
    v = _clean(v)
    if v is None:
        return datetime.now(timezone.utc).replace(tzinfo=None)
    if isinstance(v, datetime):
        return v.replace(tzinfo=None)
    if isinstance(v, str):
        for fmt in (
            "%m/%d/%y %H:%M",
            "%m/%d/%Y %H:%M",
            "%Y-%m-%d %H:%M:%S",
        ):
            try:
                return datetime.strptime(v.strip(), fmt)
            except ValueError:
                pass
    raise ValueError(f"Could not parse datetime value: {v!r}")



def _normalize_field_name(name: str) -> str:
    """Normalize a field name for case/punctuation-insensitive matching."""
    return "".join(ch for ch in str(name).upper() if ch.isalnum())


def _source_candidates_for_db_field(db_field: str) -> List[str]:
    """Return source-field candidates for an INV_COMP_LOT column.

    INV_COMP_LOT columns conventionally use an ICL_ prefix while source files
    generally contain the same logical field without that prefix. Matching also
    tolerates punctuation/case differences.
    """
    db_field = str(db_field)
    candidates = [db_field]
    if db_field.upper().startswith("ICL_"):
        candidates.append(db_field[4:])
    return candidates


def _build_dynamic_field_map(
    inv_comp_lot: Table,
    source_fields: Sequence[str],
) -> Dict[str, str]:
    """Map INV_COMP_LOT column names to source columns when unambiguous."""
    normalized_source: Dict[str, List[str]] = {}
    for field in source_fields:
        normalized_source.setdefault(_normalize_field_name(field), []).append(field)

    field_map: Dict[str, str] = {}
    for column in inv_comp_lot.columns:
        db_field = column.name
        matches: List[str] = []
        for candidate in _source_candidates_for_db_field(db_field):
            matches.extend(normalized_source.get(_normalize_field_name(candidate), []))

        # De-duplicate while preserving source order.
        unique_matches = list(dict.fromkeys(matches))
        if len(unique_matches) == 1:
            field_map[db_field] = unique_matches[0]
        elif len(unique_matches) > 1:
            # Prefer the exact normalized DB name over the prefix-stripped name.
            exact = normalized_source.get(_normalize_field_name(db_field), [])
            if len(exact) == 1:
                field_map[db_field] = exact[0]

    return field_map


def _convert_for_column(value: Any, column: Any) -> Any:
    """Convert a reconciled source value to a value suitable for a DB column."""
    value = _clean(value)
    if value is None:
        return None
    try:
        py_type = column.type.python_type
    except (AttributeError, NotImplementedError):
        return value

    if py_type is int:
        return _to_int(value)
    if py_type is float:
        return _to_float(value)
    if py_type is datetime:
        return _to_datetime(value)
    if py_type.__name__ == "date":
        return _to_date(value)
    return value


def _dynamic_lot_payload(
    inv_comp_lot: Table,
    group: Sequence[SourceRow],
    source_fields: Sequence[str],
    *,
    mixture_ic_id: int,
    assigned_lot_num: int,
    mixture_alias: str,
    mixture_lot_class: str,
    screening_center: str,
    ins_by: int,
    now: datetime,
) -> Dict[str, Any]:
    """Build an INV_COMP_LOT payload from reconciled source metadata.

    A source value is copied only when all non-null values within the mixture
    group agree. This prevents silently selecting one component's metadata when
    the source contains conflicting values.
    """
    field_map = _build_dynamic_field_map(inv_comp_lot, source_fields)
    payload: Dict[str, Any] = {}

    # DB-generated / controlled fields are always set explicitly.
    controlled = {
        "ICL_ID": None,
        "ICL_IC_ID": mixture_ic_id,
        "ICL_LOT_NUM": assigned_lot_num,
        "ICL_ALIAS": mixture_alias,
        "ICL_STATUS_ID": DEFAULT_STATUS_ID,
        "ICL_INS_BY": ins_by,
        "ICL_INS_DATE": now,
        "class": mixture_lot_class,
        "ICL_SCREENING_CENTER": screening_center,
    }

    for column in inv_comp_lot.columns:
        db_field = column.name
        if db_field in controlled:
            value = controlled[db_field]
            if value is not None:
                payload[db_field] = value
            continue
        if db_field not in field_map:
            continue

        source_field = field_map[db_field]
        reconciled = _require_same(group, source_field)
        if reconciled is None:
            continue
        payload[db_field] = _convert_for_column(reconciled, column)

    return payload

def load_rows_from_excel(excel_path: str, sheet_name: Optional[str] = None) -> List[SourceRow]:
    wb = openpyxl.load_workbook(excel_path, data_only=True)
    ws = wb[sheet_name] if sheet_name else wb[wb.sheetnames[0]]
    headers = [ws.cell(1, c).value for c in range(1, ws.max_column + 1)]
    headers = [x.upper() for x in headers]
    rows: List[SourceRow] = []
    for r in range(2, ws.max_row + 1):
        vals = [ws.cell(r, c).value for c in range(1, ws.max_column + 1)]
        row = {headers[i]: vals[i] for i in range(len(headers)) if headers[i] is not None}
        if all(_clean(v) is None for v in row.values()):
            continue
        rows.append(SourceRow(row))
    return rows


def load_rows_from_csv(csv_path: str) -> Tuple[List[SourceRow], List[str]]:
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file has no header: {csv_path}")
        fieldnames = list(reader.fieldnames)
        fieldnames = [x.upper() for x in fieldnames]
        rows = []
        for raw in reader:
            row = {k: v for k, v in raw.items() if k is not None}
            if all(_clean(v) is None for v in row.values()):
                continue
            rows.append(SourceRow(row))
    return rows, fieldnames


def load_rows_from_dataframe(df: pd.DataFrame) -> Tuple[List[SourceRow], List[str]]:
    """Convert a pandas DataFrame to the internal row representation."""
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"Expected pandas.DataFrame, got {type(df).__name__}")
    fieldnames = [str(c) for c in df.columns]
    fieldnames = [x.upper() for x in fieldnames]
    normalized = df.copy()
    normalized.columns = fieldnames
    normalized = normalized.where(pd.notna(normalized), None)
    rows = [
        SourceRow(record)
        for record in normalized.to_dict(orient="records")
        if any(_clean(v) is not None for v in record.values())
    ]
    return rows, fieldnames


def _is_csv(path: str) -> bool:
    return Path(path).suffix.lower() == ".csv"


def load_rows_from_file(path: str, sheet_name: Optional[str] = None) -> Tuple[List[SourceRow], Optional[List[str]]]:
    """Load either CSV or Excel input.

    Returns rows plus CSV fieldnames (or None for Excel).
    """
    if _is_csv(path):
        return load_rows_from_csv(path)
    return load_rows_from_excel(path, sheet_name), None


def group_rows_by_mixture(rows: Iterable[SourceRow]) -> Dict[str, List[SourceRow]]:
    grouped: Dict[str, List[SourceRow]] = {}
    for row in rows:
        mixture_alias = _clean(row.get("MIXTURE_PART_OF"))
        if mixture_alias is None:
            continue
        grouped.setdefault(str(mixture_alias), []).append(row)
    return grouped


def _require_same(group: Sequence[SourceRow], key: str) -> Any:
    values = [_clean(r.get(key.upper())) for r in group]
    non_null = [v for v in values if v is not None]
    if not non_null:
        return None
    if any(v != non_null[0] for v in non_null[1:]):
        raise ValueError(f"Inconsistent {key} values in one mixture group: {values}")
    return non_null[0]


def _next_lot_num(session: Session, inv_comp_lot: Table, mixture_ic_id: int) -> int:
    stmt = select(inv_comp_lot.c.ICL_LOT_NUM).where(inv_comp_lot.c.ICL_IC_ID == mixture_ic_id)
    existing = [row[0] for row in session.execute(stmt).all() if row[0] is not None]
    return (max(existing) + 1) if existing else 1


def _find_existing_lot_for_constituents(
    session: Session,
    inv_comp_lot: Table,
    inv_compound_lot_mixture: Table,
    mixture_ic_id: int,
    constituent_lot_ids: Sequence[int],
) -> Optional[Tuple[int, int]]:
    """Return (ICL_ID, ICL_LOT_NUM) when this exact constituent-lot tuple already exists.

    Constituent IDs are compared as a sorted tuple so row insertion order does not matter,
    while duplicate IDs are still preserved as part of the exact tuple.
    """
    target = tuple(sorted(int(x) for x in constituent_lot_ids))
    if not target:
        return None

    stmt = (
        select(
            inv_comp_lot.c.ICL_ID,
            inv_comp_lot.c.ICL_LOT_NUM,
            inv_compound_lot_mixture.c.ICLM_CONSTITUENT_ID,
        )
        .select_from(
            inv_comp_lot.join(
                inv_compound_lot_mixture,
                inv_comp_lot.c.ICL_ID == inv_compound_lot_mixture.c.ICLM_MIXTURE_ID,
            )
        )
        .where(inv_comp_lot.c.ICL_IC_ID == mixture_ic_id)
    )

    grouped: Dict[int, List[int]] = {}
    lot_nums: Dict[int, int] = {}
    for lot_id, lot_num, constituent_id in session.execute(stmt).all():
        if lot_id is None or constituent_id is None:
            continue
        lot_id = int(lot_id)
        grouped.setdefault(lot_id, []).append(int(constituent_id))
        lot_nums[lot_id] = int(lot_num) if lot_num is not None else 0

    for lot_id, constituent_ids in grouped.items():
        if tuple(sorted(constituent_ids)) == target:
            return lot_id, lot_nums[lot_id]
    return None


def _find_mixture_compound(
    session: Session,
    inv_compound: Table,
    mixture_alias: str,
) -> Optional[int]:
    stmt = select(inv_compound.c.IC_ID).where(
        inv_compound.c.IC_ALIAS_ID == mixture_alias
    )
    matches = [row[0] for row in session.execute(stmt).all() if row[0] is not None]
    if not matches:
        return None
    if len(set(matches)) > 1:
        raise ValueError(
            f"Multiple INV_COMPOUND rows found for mixture alias {mixture_alias!r}: {matches}"
        )
    return int(matches[0])


def _get_mixture_component_ids(
    session: Session,
    inv_compound_mixture: Table,
    mixture_ic_id: int,
) -> List[int]:
    stmt = select(inv_compound_mixture.c.ICM_PURE_COMPOUND_ID).where(
        inv_compound_mixture.c.ICM_MIXTURE_ID == mixture_ic_id
    )
    return [
        int(row[0])
        for row in session.execute(stmt).all()
        if row[0] is not None
    ]


def _insert_mixture_compound(
    session: Session,
    inv_compound: Table,
    group: Sequence[SourceRow],
    mixture_alias: str,
    screening_center: str,
    ins_by: int,
    sync_ic_compound_id: bool,
) -> int:
    now = datetime.now()
    compound_payload = {
        "IC_COMPOUND_ID": None,
        "IC_ALIAS_ID": mixture_alias,
        "IC_SMILES": None,
        "IC_INCHI": None,
        "IC_INCHI_KEY": None,
        "IC_SEQUENCE": None,
        "IC_HELM": None,
        "IC_NOTES": None,
        "IC_STEREO": None,
        "IC_MIXTURE": None,
        "IC_INS_BY": ins_by,
        "IC_INS_DATE": now,
        "class": MIXTURE_COMPOUND_CLASS,
        "IC_SCREENING_CENTER": screening_center,
    }
    compound_result = session.execute(
        insert(inv_compound).values(**compound_payload)
    )
    mixture_ic_id = int(compound_result.inserted_primary_key[0])

    if sync_ic_compound_id:
        session.execute(
            update(inv_compound)
            .where(inv_compound.c.IC_ID == mixture_ic_id)
            .values(IC_COMPOUND_ID=mixture_ic_id)
        )

    return mixture_ic_id

def insert_one_mixture(
    session: Session,
    tables: Dict[str, Table],
    group: Sequence[SourceRow],
    *,
    sync_ic_compound_id: bool = True,
    mixture_lot_class: str = DEFAULT_MIXTURE_LOT_CLASS,
    mixture_lot_num: Optional[int] = None,
    mixture_alias: Optional[str] = None,
    component_compound_id_col: str = "SMDC",
    component_lot_id_col: str = "LOT_ID",
    concentration_col: str = "Concentration_uM",
    volume_col: str = "Volume_uL",
    ins_by: int = 569,
) -> Dict[str, int]:
    if not group:
        raise ValueError("Empty mixture group")

    rep = group[0]
    mixture_alias = str(mixture_alias or _clean(rep.get("MIXTURE_PART_OF")) or "").strip()
    if not mixture_alias:
        raise ValueError("mixture alias is required; expected mixture_part_of to be populated")

    inv_compound = tables["INV_COMPOUND"]
    inv_comp_lot = tables["INV_COMP_LOT"]
    inv_compound_mixture = tables["INV_COMPOUND_MIXTURE"]
    inv_compound_lot_mixture = tables["INV_COMPOUND_LOT_MIXTURE"]

    screening_center = _clean(rep.get("SCREENING_CENTER")) or DEFAULT_SCREENING_CENTER

    # Path 1: reuse an existing mixture compound.
    mixture_ic_id = _find_mixture_compound(session, inv_compound, mixture_alias)
    if mixture_ic_id is not None:
        existing_component_ids = _get_mixture_component_ids(
            session, inv_compound_mixture, mixture_ic_id
        )
        if not existing_component_ids:
            raise ValueError(
                f"Mixture IC_ID={mixture_ic_id} exists but has no rows in INV_COMPOUND_MIXTURE"
            )
    else:
        # Path 2: this is a new mixture, so create the compound and its component links.
        mixture_ic_id = _insert_mixture_compound(
            session,
            inv_compound,
            group,
            mixture_alias,
            screening_center,
            ins_by,
            sync_ic_compound_id,
        )
        existing_component_ids = []
        for row in group:
            component_compound_id = _to_int(row.get(component_compound_id_col.upper()))
            if component_compound_id is None:
                raise ValueError(
                    f"Missing component compound id in column {component_compound_id_col!r}"
                )
            existing_component_ids.append(component_compound_id)
            session.execute(
                insert(inv_compound_mixture).values(
                    ICM_MIXTURE_ID=mixture_ic_id,
                    ICM_PURE_COMPOUND_ID=component_compound_id,
                )
            )

    source_component_ids: List[int] = []
    for row in group:
        component_compound_id = _to_int(row.get(component_compound_id_col.upper()))
        if component_compound_id is None:
            raise ValueError(
                f"Missing component compound id in column {component_compound_id_col!r}"
            )
        source_component_ids.append(component_compound_id)

    if sorted(source_component_ids) != sorted(existing_component_ids):
        raise ValueError(
            f"Mixture {mixture_alias!r} does not match INV_COMPOUND_MIXTURE for "
            f"IC_ID={mixture_ic_id}. Source component IDs={sorted(source_component_ids)}, "
            f"database component IDs={sorted(existing_component_ids)}"
        )

    source_component_lot_ids: List[int] = []
    for row in group:
        component_lot_id = _to_int(row.get(component_lot_id_col.upper()))
        if component_lot_id is None:
            raise ValueError(f"Missing component lot id in column {component_lot_id_col!r}")
        source_component_lot_ids.append(component_lot_id)

    # Reuse an existing lot when its constituent-lot tuple is exactly the same.
    existing_lot = _find_existing_lot_for_constituents(
        session,
        inv_comp_lot,
        inv_compound_lot_mixture,
        mixture_ic_id,
        source_component_lot_ids,
    )
    if existing_lot is not None:
        existing_icl_id, existing_lot_num = existing_lot
        return {
            "MIXTURE_IC_ID": mixture_ic_id,
            "MIXTURE_ICL_ID": existing_icl_id,
            "MIXTURE_LOT_NUM": int(existing_lot_num),
            "INSERTED": 0,
        }

    # No lot with this exact constituent-lot tuple exists, so create a new lot.
    next_lot_num = _next_lot_num(session, inv_comp_lot, mixture_ic_id)
    assigned_lot_num = next_lot_num
    if mixture_lot_num is not None and next_lot_num == 1:
        assigned_lot_num = mixture_lot_num

    now = datetime.now()
    source_fields = list(group[0].raw.keys())
    lot_payload = _dynamic_lot_payload(
        inv_comp_lot,
        group,
        source_fields,
        mixture_ic_id=mixture_ic_id,
        assigned_lot_num=assigned_lot_num,
        mixture_alias=mixture_alias,
        mixture_lot_class=mixture_lot_class,
        screening_center=screening_center,
        ins_by=ins_by,
        now=now,
    )
    lot_result = session.execute(insert(inv_comp_lot).values(**lot_payload))
    mixture_icl_id = int(lot_result.inserted_primary_key[0])

    for row, component_lot_id in zip(group, source_component_lot_ids):
        session.execute(
            insert(inv_compound_lot_mixture).values(
                ICLM_MIXTURE_ID=mixture_icl_id,
                ICLM_CONSTITUENT_ID=component_lot_id,
                ICLM_FINAL_CONCENTRATION=_to_float(row.get(concentration_col.upper())),
                ICLM_VOLUME_UL=_to_float(row.get(volume_col.upper())),
            )
        )

    return {
        "MIXTURE_IC_ID": mixture_ic_id,
        "MIXTURE_ICL_ID": mixture_icl_id,
        "MIXTURE_LOT_NUM": int(assigned_lot_num),
        "INSERTED": 1,
    }


def build_tables(engine: Engine) -> Dict[str, Table]:
    metadata = MetaData()
    inv_compound = Table("INV_COMPOUND", metadata, autoload_with=engine)
    inv_comp_lot = Table("INV_COMP_LOT", metadata, autoload_with=engine)
    inv_compound_mixture = Table("INV_COMPOUND_MIXTURE", metadata, autoload_with=engine)
    inv_compound_lot_mixture = Table("INV_COMPOUND_LOT_MIXTURE", metadata, autoload_with=engine)
    return {
        "INV_COMPOUND": inv_compound,
        "INV_COMP_LOT": inv_comp_lot,
        "INV_COMPOUND_MIXTURE": inv_compound_mixture,
        "INV_COMPOUND_LOT_MIXTURE": inv_compound_lot_mixture,
    }


def _make_output_row(
    representative: SourceRow,
    result: Dict[str, int],
    *,
    mixture_alias: str,
    output_fields: Sequence[str],
    source_format: str,
) -> Dict[str, Any]:
    """Create one source row representing the newly inserted mixture."""
    row = dict(representative.raw)
    row["LOT_ID"] = result["MIXTURE_ICL_ID"]
    row["ICL_IC_ID"] = result["MIXTURE_IC_ID"]
    row["ALIAS"] = mixture_alias
    row["ICL_LOT_NUM"] = result["MIXTURE_LOT_NUM"]
    row["COMPOUNDINDEX"] = None
    row["MIXTURE_PART_OF"] = None
    row["MIXTURE_PART_COUNT"] = None

    # Excel output is intentionally not supported: the returned source is CSV so it
    # can be consumed directly by the same pipeline that produced new_source.csv.
    if source_format != "csv":
        raise ValueError("Rewritten output is only supported as CSV")
    return {field: row.get(field, "") for field in output_fields}


def build_rewritten_csv_rows(
    rows: Sequence[SourceRow],
    grouped: Dict[str, List[SourceRow]],
    results: Dict[str, Dict[str, int]],
    fieldnames: Sequence[str],
) -> List[Dict[str, Any]]:
    """Keep ordinary rows and append one coded mixture row per mixture."""
    mixture_row_ids = {id(member) for group in grouped.values() for member in group}
    output_rows: List[Dict[str, Any]] = []

    for row in rows:
        if id(row) in mixture_row_ids:
            continue
        output_rows.append({field: row.raw.get(field, "") for field in fieldnames})

    # User requested the newly coded mixture rows to be appended to the source.
    for mixture_alias, group in grouped.items():
        output_rows.append(
            _make_output_row(
                group[0],
                results[mixture_alias],
                mixture_alias=mixture_alias,
                output_fields=fieldnames,
                source_format="csv",
            )
        )
    return output_rows


def write_csv(path: str, fieldnames: Sequence[str], rows: Sequence[Dict[str, Any]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(output)


def insert_mixtures(
    source,
    engine: Optional = None,
    host: Optional[str] = "localhost",
    db: Optional[str] = "QB3_HiTS_nightly",
    sheet_name: Optional[str] = None,
    output_path: Optional[str] = None,
    sync_ic_compound_id: bool = True,
    mixture_lot_class: str = DEFAULT_MIXTURE_LOT_CLASS,
    mixture_lot_num: Optional[int] = None,
    compound_id_col: Optional[str] = None,
    lot_id_col: Optional[str] = None,
    concentration_col: str = "Concentration_mM",
    volume_col: str = "Volume_uL",
    ins_by: int = 569,
) -> pd.DataFrame:
    """Insert mixtures from either a DataFrame or a CSV/Excel path.

    Returns a DataFrame containing the non-mixture rows plus one newly coded
    row per mixture. If ``output_path`` is provided, that returned DataFrame is
    also written to CSV or Excel based on the output suffix.
    """
    if isinstance(source, pd.DataFrame):
        rows, fieldnames = load_rows_from_dataframe(source)
        source_is_csv_like = True
    elif isinstance(source, (str, Path)):
        rows, loaded_fields = load_rows_from_file(str(source), sheet_name)
        if loaded_fields is None:
            # Preserve Excel columns for the returned DataFrame.
            fieldnames = list(rows[0].raw.keys()) if rows else []
            source_is_csv_like = False
        else:
            fieldnames = loaded_fields
            source_is_csv_like = True
    else:
        raise TypeError("source must be a pandas DataFrame or a CSV/Excel path")

    grouped = group_rows_by_mixture(rows)
    print("Number of mixtures to process: ", len(grouped))
    if not grouped:
        raise ValueError("No rows with mixture_part_of values were found")

    # new_source-style DataFrames/CSVs use these columns. Excel retains the
    # legacy defaults unless explicitly overridden.
    if compound_id_col is None:
        compound_id_col = "ICL_IC_ID" if "ICL_IC_ID" in fieldnames else "SMDC"
    if lot_id_col is None:
        lot_id_col = "lot_id" if "lot_id" in fieldnames else "LOT_ID"

    if engine is None:
        engine = create_mysql_cnx(host=host, db=db)
    with engine.connect() as conn:
        print(conn.exec_driver_sql("select 'connection ok'").scalar())

    tables = build_tables(engine)
    results: Dict[str, Dict[str, int]] = {}

    with Session(engine) as session:
        try:
            for mixture_alias, group in grouped.items():
                result = insert_one_mixture(
                    session, tables, group,
                    sync_ic_compound_id=sync_ic_compound_id,
                    mixture_lot_class=mixture_lot_class,
                    mixture_lot_num=mixture_lot_num,
                    mixture_alias=mixture_alias,
                    component_compound_id_col=compound_id_col,
                    component_lot_id_col=lot_id_col,
                    concentration_col=concentration_col,
                    volume_col=volume_col,
                    ins_by=ins_by,
                )
                results[mixture_alias] = result
                if result.get("INSERTED"):
                    print(
                        f"Inserted mixture lot {mixture_alias}: "
                        f"IC_ID={result['MIXTURE_IC_ID']} "
                        f"ICL_ID={result['MIXTURE_ICL_ID']} "
                        f"LOT_NUM={result['MIXTURE_LOT_NUM']}"
                    )
                else:
                    print(
                        f"Mixture lot already exists; inserted nothing for {mixture_alias}: "
                        f"IC_ID={result['MIXTURE_IC_ID']} "
                        f"ICL_ID={result['MIXTURE_ICL_ID']} "
                        f"LOT_NUM={result['MIXTURE_LOT_NUM']}"
                    )

            mixture_row_ids = {id(member) for group in grouped.values() for member in group}
            output_rows = [
                {field: row.raw.get(field, "") for field in fieldnames}
                for row in rows if id(row) not in mixture_row_ids
            ]
            for mixture_alias, group in grouped.items():
                row = dict(group[0].raw)
                row["LOT_ID"] = results[mixture_alias]["MIXTURE_ICL_ID"]
                row["ICL_IC_ID"] = results[mixture_alias]["MIXTURE_IC_ID"]
                row["ALIAS"] = mixture_alias
                row["ICL_LOT_NUM"] = results[mixture_alias]["MIXTURE_LOT_NUM"]
                row["COMPOUNDINDEX"] = None
                row["MIXTURE_PART_OF"] = None
                row["MIXTURE_PART_COUNT"] = None
                output_rows.append({field: row.get(field, "") for field in fieldnames})

            result_df = pd.DataFrame(output_rows, columns=fieldnames)
            if output_path is not None:
                out = Path(output_path)
                out.parent.mkdir(parents=True, exist_ok=True)
                if out.suffix.lower() == ".csv":
                    result_df.to_csv(out, index=False)
                elif out.suffix.lower() in {".xlsx", ".xlsm"}:
                    result_df.to_excel(out, index=False)
                else:
                    raise ValueError("output_path must end in .csv, .xlsx, or .xlsm")
            session.commit()
        except Exception:
            session.rollback()
            raise

    return result_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", "--excel", dest="source", required=True,
        help="Source CSV or Excel workbook containing component rows"
    )
    parser.add_argument(
        "--output", default=None,
        help="For CSV input, write a rewritten source with mixture_part rows removed and coded mixture rows appended"
    )
    parser.add_argument("--sheet", default=None, help="Excel worksheet name; default is the first sheet")
    parser.add_argument("--host", default="localhost", help="Database host (localhost or localhost)")
    parser.add_argument("--db", default="QB3_HiTS_nightly", help="Database name")
    parser.add_argument("--no-sync-ic-compound-id", action="store_true", help="Leave IC_COMPOUND_ID NULL instead of matching IC_ID")
    parser.add_argument("--mixture-lot-class", default=DEFAULT_MIXTURE_LOT_CLASS, help="Class to use for the mixture lot row")
    parser.add_argument("--mixture-lot-num", type=int, default=None, help="Explicit lot number for every new mixture; default is next available")
    parser.add_argument("--compound-id-col", default=None, help="Component compound ID column; CSV default is ICL_IC_ID, Excel default is SMDC")
    parser.add_argument("--lot-id-col", default=None, help="Component lot ID column; CSV default is lot_id, Excel default is LOT_ID")
    parser.add_argument("--concentration-col", default="Concentration_mM", help="Column for lot-mixture concentration")
    parser.add_argument("--volume-col", default="Volume_uL", help="Column for lot-mixture volume")
    parser.add_argument("--ins-by", type=int, default=569, help="User ID recorded in IC_INS_BY and ICL_INS_BY (default: 569)")
    args = parser.parse_args()

    try:
        insert_mixtures(
            source=args.source,
            output_path=args.output,
            sheet_name=args.sheet,
            host=args.host,
            db=args.db,
            sync_ic_compound_id=not args.no_sync_ic_compound_id,
            mixture_lot_class=args.mixture_lot_class,
            mixture_lot_num=args.mixture_lot_num,
            compound_id_col=args.compound_id_col,
            lot_id_col=args.lot_id_col,
            concentration_col=args.concentration_col,
            volume_col=args.volume_col,
            ins_by=args.ins_by,
        )
    except Exception as exc:
        raise SystemExit(f"ERROR: {exc}") from exc


if __name__ == "__main__":
    main()