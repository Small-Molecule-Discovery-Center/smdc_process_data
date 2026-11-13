"""Reusable helpers for registering plates and wells.

The notebook currently mixes dataframe reshaping, layout mapping, validation,
and SQL insertion. These helpers split those concerns into reusable pieces so
callers can drive the pipeline with a config dict instead of a long parameter
list.

Typical usage:

    plates, wells = prepare_inventory_frames(
        source,
        plate_inputs={"IP_GRID_NUM": "plate", "IP_BARCODE": "plate", "IP_ICO_ID": 267, ...},
        well_inputs={"IW_COORDINATES": "well", "IW_CONTAM": 0, ...},
        size=384,
        engine=engine
    )
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import string
import sys
import re

import pandas as pd
from sqlalchemy import MetaData, Table, and_, or_, select


WELL_ROWS = {
    96: list(string.ascii_uppercase[:8]),
    384: list(string.ascii_uppercase[:16]),
    1536: list(string.ascii_uppercase) + ['AA', 'AB', 'AC', 'AD', 'AE', 'AF'],
}

WELL_COLS = {
    96: list(range(1, 13)),
    384: list(range(1, 25)),
    1536: list(range(1, 49)),
}

MAP_PROFILE_FILES = {
    "interleave": "map_96_to_384_interleaved.csv",
    "quadrant": "map_384_1536_quadrant.csv",
}

# Standard normalized map columns expected in all layout map files
MAP_STANDARD_COLS = ("source_well", "source_quadrant", "dest_well")


def now_string(fmt: str = "%Y-%m-%d %H:%M:%S") -> str:
    """Return the current local timestamp formatted for IP_INS_DATE / IW_INS_DATE.

    Parameters:
        fmt: strftime format string; defaults to the ``YYYY-MM-DD HH:MM:SS`` form
            used by the source notebooks' ``stringnow()`` helper.
    """
    return datetime.now().strftime(fmt)


def validate_required_columns(frame: pd.DataFrame, required: Sequence[str], label: str) -> None:
    """Raise when any required columns are missing from a dataframe.

    Parameters:
        frame: DataFrame to validate.
        required: Column names that must be present.
        label: Friendly name used in error messages.
    """
    missing = [col for col in required if col not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def validate_unique(frame: pd.DataFrame, keys: Sequence[str], label: str) -> None:
    """Raise when duplicate key combinations exist in a dataframe.

    Parameters:
        frame: DataFrame to validate.
        keys: Column names that define a uniqueness key.
        label: Friendly name used in error messages.
    """
    duplicated = frame.duplicated(subset=list(keys), keep=False)
    if duplicated.any():
        dup_rows = frame.loc[duplicated, list(keys)].drop_duplicates()
        raise ValueError(f"{label} has duplicate keys: {dup_rows.to_dict(orient='records')}")


def table_column_names(table: Table) -> list[str]:
    """Return ordered SQLAlchemy table column names.

    Parameters:
        table: Reflected or declared SQLAlchemy table.
    """
    return [column.name for column in table.columns]


def reflect_inventory_tables(engine: Engine) -> tuple[Table, Table]:
    """Reflect INV_PLATE and INV_WELL from the connected database.

    Parameters:
        engine: SQLAlchemy engine connected to the target schema.
    """
    metadata = MetaData()
    plate = Table("INV_PLATE", metadata, autoload_with=engine)
    well = Table("INV_WELL", metadata, autoload_with=engine)
    return plate, well


def reflect_table(engine: Engine, name: str) -> Table:
    """Reflect a single table by name from the connected database.

    Parameters:
        engine: SQLAlchemy engine connected to the target schema.
        name: Table name to reflect (e.g. ``"INV_COLLECTION"``).
    """
    return Table(name, MetaData(), autoload_with=engine)


def to_mapping(
    lookup: Mapping[Any, Any] | pd.Series | pd.DataFrame | None,
    *,
    key_col: str | None = None,
    value_col: str | None = None,
) -> dict[Any, Any]:
    """Normalize mapping-like inputs into a plain dictionary.

    Parameters:
        lookup: Mapping, Series, DataFrame, or None.
        key_col: Key column when lookup is a DataFrame.
        value_col: Value column when lookup is a DataFrame.
    """
    if lookup is None:
        return {}
    if isinstance(lookup, Mapping):
        return dict(lookup)
    if isinstance(lookup, pd.Series):
        return lookup.to_dict()
    if isinstance(lookup, pd.DataFrame):
        if key_col is None or value_col is None:
            raise ValueError("key_col and value_col are required when lookup is a DataFrame")
        validate_required_columns(lookup, [key_col, value_col], "lookup dataframe")
        return dict(zip(lookup[key_col], lookup[value_col]))
    raise TypeError(f"Unsupported lookup type: {type(lookup)!r}")


def well_coordinates(size: int) -> list[str]:
    """Generate row-major plate coordinates for a plate size.

    Parameters:
        size: Plate size; supported values are 96, 384, 1536.
    """
    if size not in WELL_ROWS or size not in WELL_COLS:
        raise ValueError(f"Unsupported plate size: {size}")
    return [f"{row}{col:02d}" for row in WELL_ROWS[size] for col in WELL_COLS[size]]


def well_index_map(size: int) -> dict[str, int]:
    """Generate coordinate-to-index mapping for a plate size.

    Parameters:
        size: Plate size; supported values are 96, 384, 1536.
    """
    return {coord: idx + 1 for idx, coord in enumerate(well_coordinates(size))}

def add_well_column_to_source(source: pd.DataFrame, row: str, col: str, well_col: str) -> pd.DataFrame:
    """Generate a dataframe column with zero-padded well coordinates from row and column values.

    Parameters:
        row: Row letter(s) (e.g., 'A', 'B', ..., 'P').
        col: Column number (1-based).
        well_col: Name of the well column in the source DataFrame.
    """
    source[well_col] = [f"{r}{c:02d}" for r, c in zip(source[row], source[col].astype(int))]
    return source

def load_layout_map(source: str | Path | pd.DataFrame) -> pd.DataFrame:
    """Load a layout map from DataFrame or CSV path.

    Parameters:
        source: DataFrame or CSV path.
    """
    if isinstance(source, pd.DataFrame):
        return source.copy()
    path = Path(source)
    # utf-8-sig strips a leading BOM (present on map_96_to_384_interleaved.csv)
    return pd.read_csv(path, encoding="utf-8-sig")


def resolve_layout_map(
    *,
    profile: str,
    custom_path: str | Path | None = None,
    base_dir: str | Path = Path(__file__).resolve().parent,
) -> pd.DataFrame:
    """Load a map artifact by profile name or custom path.

    Parameters:
        profile: One of `interleave`, `quadrant`, `custom`.
        custom_path: Required when profile is `custom`.
        base_dir: Directory containing map artifacts for named profiles.
    """
    profile_key = profile.strip().lower()
    if profile_key == "custom":
        if custom_path is None:
            raise ValueError("custom_path is required when profile='custom'")
        return load_layout_map(Path(custom_path))
    if profile_key not in MAP_PROFILE_FILES:
        raise ValueError(f"Unknown layout profile: {profile}")
    return load_layout_map(Path(base_dir) / MAP_PROFILE_FILES[profile_key])


def _natural_sort_key(value: Any) -> tuple[int, int, str]:
    """Natural sort so Q1 < Q2 < Q10 rather than lexicographic."""
    text = str(value)
    m = re.fullmatch(r"([A-Za-z]+)(\d+)", text)
    if m:
        prefix, num = m.groups()
        return (0, int(num), prefix)
    return (1, 0, text)


def _build_quadrant_assignment(
    source: pd.DataFrame,
    source_plate_col: str,
    destination_plate_names: Sequence[str],
    source_plate_order: Sequence[str] | Mapping[str, Sequence[str]] | None,
    plates_per_destination: int,
    quadrants: Sequence[str],
) -> pd.DataFrame:
    """Build a plate→(dest_plate, quadrant) assignment table.

    Returns a DataFrame with columns [source_plate_col, '_DEST_PLATE_', '_SRC_QUADRANT_'].

    The final destination may be partially filled: a group of fewer than
    ``plates_per_destination`` source plates simply occupies quadrants
    ``Q1..Qk`` and leaves the rest empty (those wells are padded later).
    """
    if len(quadrants) != plates_per_destination:
        raise ValueError(
            f"Map has {len(quadrants)} quadrant labels {list(quadrants)}, "
            f"but plates_per_destination={plates_per_destination}"
        )

    if isinstance(source_plate_order, Mapping):
        order_map = {str(k): list(v) for k, v in source_plate_order.items()}
        rows: list[tuple[str, str, str]] = []
        for dest in destination_plate_names:
            srcs = order_map.get(str(dest), [])
            if not (1 <= len(srcs) <= plates_per_destination):
                raise ValueError(
                    f"Destination '{dest}' requires 1..{plates_per_destination} source plates, got {len(srcs)}"
                )
            rows.extend((str(src), str(dest), str(q)) for src, q in zip(srcs, quadrants[:len(srcs)]))
    else:
        all_sources = sorted(
            source[source_plate_col].dropna().astype(str).unique().tolist(),
            key=_natural_sort_key,
        )
        if source_plate_order is not None:
            all_sources = [str(s) for s in source_plate_order]
        n_dest = len(destination_plate_names)
        capacity = n_dest * plates_per_destination
        min_needed = (n_dest - 1) * plates_per_destination + 1
        if not (min_needed <= len(all_sources) <= capacity):
            raise ValueError(
                f"{len(all_sources)} source plates do not fit {n_dest} destinations at "
                f"{plates_per_destination} per destination (expected {min_needed}..{capacity}); "
                f"destination count should be ceil(n_sources / {plates_per_destination})"
            )
        rows = []
        for idx, dest in enumerate(destination_plate_names):
            block = all_sources[idx * plates_per_destination:(idx + 1) * plates_per_destination]
            rows.extend((str(src), str(dest), str(q)) for src, q in zip(block, quadrants[:len(block)]))

    return pd.DataFrame(rows, columns=[source_plate_col, "_DEST_PLATE_", "_SRC_QUADRANT_"])


def apply_layout_map(
    source: pd.DataFrame,
    layout_map: pd.DataFrame,
    *,
    source_well_col: str,
    target_well_col: str,
    source_plate_col: str | None = None,
    drop_unmapped: bool = True,
) -> pd.DataFrame:
    """Merge source data with a normalized layout map and populate mapped well column.

    The layout map must have columns ``source_well``, ``source_quadrant``, and
    ``dest_well`` (the standard normalized format).  This helper performs a simple
    well-only remap (matching on ``source_well``); it does **not** reshape multiple
    source plates into a destination plate.  Multi-parent quadrant/interleave
    reshaping now lives in :func:`remap_wells_by_quadrant`.

    Parameters:
        source: Input source dataframe.
        layout_map: Normalized layout map with ``source_well``, ``source_quadrant``,
            and ``dest_well`` columns.
        source_well_col: Well column name in the user dataframe.
        target_well_col: Output column written with mapped destination well.
        source_plate_col: Optional plate identifier column, validated when supplied.
        drop_unmapped: When True (default), source rows with no matching map entry
            are silently dropped (e.g. border wells not included in an interleave
            transfer).  When False, unmapped wells raise a ``ValueError``.
    """
    validate_required_columns(source, [source_well_col], "source dataframe")
    validate_required_columns(layout_map, list(MAP_STANDARD_COLS), "layout map")

    left = source.copy()
    right = layout_map.copy()

    if source_plate_col is not None:
        validate_required_columns(left, [source_plate_col], "source dataframe")

    merged = left.merge(
        right,
        how="left",
        left_on=source_well_col,
        right_on="source_well",
    )

    merged = merged.drop_duplicates()
    mapped_missing = merged["dest_well"].isna()
    if mapped_missing.any():
        if not drop_unmapped:
            diag_cols = [c for c in [source_plate_col, source_well_col] if c]
            missing = merged.loc[mapped_missing, diag_cols].drop_duplicates()
            raise ValueError(f"Layout map did not resolve all wells: {missing.to_dict(orient='records')}")
        merged = merged.loc[~mapped_missing].copy()

    merged[target_well_col] = merged["dest_well"]
    drop_cols = [c for c in ["source_well", "source_quadrant", "dest_well"] if c in merged.columns and c != target_well_col]
    merged = merged.drop(columns=drop_cols)
    return merged


COLLECTION_TABLE_NAME = "INV_COLLECTION"


def _resolve_collection(
    plates: pd.DataFrame,
    *,
    engine: Engine | None,
) -> pd.DataFrame:
    """Fill whichever of IP_ICO_ID / IP_COLLECTION is missing from INV_COLLECTION.

    ``IP_ICO_ID`` (``ICO_ID``) and ``IP_COLLECTION`` (``ICO_NAME``) are 1:1 on
    INV_COLLECTION.  Supply exactly one and the other is looked up (which reflects
    INV_COLLECTION from ``engine``); supply both and nothing is looked up; supply
    neither and this raises.
    """
    has_ico = "IP_ICO_ID" in plates.columns
    has_collection = "IP_COLLECTION" in plates.columns
    if not has_ico and not has_collection:
        raise ValueError("Supply IP_ICO_ID or IP_COLLECTION in plate_inputs")
    if has_ico and has_collection:
        return plates
    if engine is None:
        raise ValueError(
            "engine is required to resolve IP_ICO_ID <-> IP_COLLECTION from INV_COLLECTION"
        )
    cmap = fetch_collection_map(engine)
    if has_ico:
        id_to_name = dict(zip(cmap["ICO_ID"], cmap["ICO_NAME"]))
        plates["IP_COLLECTION"] = plates["IP_ICO_ID"].map(id_to_name)
        unresolved = plates.loc[plates["IP_COLLECTION"].isna(), "IP_ICO_ID"].dropna().drop_duplicates().tolist()
        if unresolved:
            raise ValueError(f"No INV_COLLECTION name for ICO_ID: {unresolved}")
    else:
        name_to_id = dict(zip(cmap["ICO_NAME"], cmap["ICO_ID"]))
        plates["IP_ICO_ID"] = plates["IP_COLLECTION"].map(name_to_id)
        unresolved = plates.loc[plates["IP_ICO_ID"].isna(), "IP_COLLECTION"].dropna().drop_duplicates().tolist()
        if unresolved:
            raise ValueError(f"No INV_COLLECTION id for name: {unresolved}")
    return plates


def _normalize_column_name(value: Any) -> str:
    """Normalize a source/database column name for tolerant name matching."""
    return re.sub(r"[^A-Za-z0-9]+", "", str(value)).upper()


def _find_matching_source_column(source: pd.DataFrame, target: str) -> str | None:
    """Find a source column by exact name first, then normalized name matching."""
    if target in source.columns:
        return target
    norm = _normalize_column_name(target)
    matches = [c for c in source.columns if _normalize_column_name(c) == norm]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous source columns for {target!r}: {matches}")
    return None


def _resolve_value_from_source(
    source: pd.DataFrame,
    value: Any,
    *,
    group_key: pd.Series | None = None,
    label: str = "input",
) -> pd.Series | Any:
    """Resolve a mapping value as a source column, keyed mapping, or literal."""
    if isinstance(value, Mapping):
        if group_key is None:
            raise ValueError(f"{label} is a mapping but no grouping key was supplied")
        keyed = {str(k): v for k, v in value.items()}
        out = group_key.astype(str).map(keyed)
        if out.isna().any():
            missing = group_key.loc[out.isna()].drop_duplicates().tolist()
            raise ValueError(f"{label} is missing entries for: {missing}")
        return out
    if isinstance(value, str):
        col = _find_matching_source_column(source, value)
        if col is not None:
            return source[col]
    return value


# Well fields that hold numbers. Each may be given as a literal (broadcast to every
# row), a source column name (per-row values, e.g. a dose-response concentration
# column), or - for daughter wells - a keyed mapping.
# IW_FINAL_CONC_NM (final assay concentration, nM) can be used in place of, or
# alongside, IW_CONC_UM; no unit conversion is done between them (uM vs nM is the
# caller's responsibility). Rule for compound wells: at least one concentration must be
# filled, and IW_VOL_UL may be blank only where IW_FINAL_CONC_NM is filled.
NUMERIC_WELL_FIELDS = ("IW_VOL_UL", "IW_CONC_UM", "IW_FINAL_CONC_NM")


def _coerce_numeric_well_input(source: pd.DataFrame, key: str, value: Any) -> Any:
    """Typo guard for numeric well inputs given as strings.

    A string is first tried as a source column name (exact, then normalized). If it is
    not a column it must parse as a number; otherwise raise instead of silently
    broadcasting a misspelled column name as a literal into ``IW_VOL_UL``/``IW_CONC_UM``.
    """
    if isinstance(value, str) and _find_matching_source_column(source, value) is None:
        try:
            return float(value)
        except ValueError:
            raise ValueError(
                f"{key} was given {value!r}, which is neither a source column nor a number. "
                f"Source columns: {list(source.columns)}"
            ) from None
    return value


def _validate_numeric_well_fields(
    wells: pd.DataFrame,
    keys: Iterable[str],
    *,
    enforce_concentration_rule: bool = False,
    compound_col: str = "IW_ICL_ID",
) -> pd.DataFrame:
    """Coerce numeric well fields to numbers and sanity check them.

    Every field in ``keys`` must be numeric and non-negative (0 is allowed, e.g. a
    vehicle-only well); nulls are allowed here.

    With ``enforce_concentration_rule`` the following is checked on wells that carry
    a compound (``compound_col`` not null); controls / padding wells are exempt:

    * at least one of ``IW_CONC_UM`` / ``IW_FINAL_CONC_NM`` is filled (both may be);
    * ``IW_VOL_UL`` is filled, unless ``IW_FINAL_CONC_NM`` is filled for that well.
    """
    out = wells
    for key in keys:
        if key not in out.columns:
            continue
        raw = out[key]
        numeric = pd.to_numeric(raw, errors="coerce")
        bad = raw.notna() & numeric.isna()
        if bad.any():
            raise ValueError(f"{key} has non-numeric values: {raw[bad].drop_duplicates().head(5).tolist()}")
        if (numeric < 0).any():
            raise ValueError(f"{key} has negative values: {numeric[numeric < 0].drop_duplicates().head(5).tolist()}")
        out[key] = numeric

    if enforce_concentration_rule:
        def filled(col: str) -> pd.Series:
            if col not in out.columns:
                return pd.Series(False, index=out.index)
            return pd.to_numeric(out[col], errors="coerce").notna()

        has_compound = out[compound_col].notna() if compound_col in out.columns else pd.Series(True, index=out.index)
        where = "IW_COORDINATES" if "IW_COORDINATES" in out.columns else None

        def describe(mask: pd.Series) -> str:
            sample = out.loc[mask, where].head(5).tolist() if where else out.index[mask][:5].tolist()
            return f"{int(mask.sum())} wells with a compound ({compound_col}), e.g. {sample}"

        final_nm = filled("IW_FINAL_CONC_NM")
        no_conc = has_compound & ~(filled("IW_CONC_UM") | final_nm)
        if no_conc.any():
            raise ValueError(
                "Neither IW_CONC_UM nor IW_FINAL_CONC_NM is filled for " + describe(no_conc)
                + ". At least one is required (both are allowed)."
            )
        no_vol = has_compound & ~filled("IW_VOL_UL") & ~final_nm
        if no_vol.any():
            raise ValueError(
                "IW_VOL_UL is blank for " + describe(no_vol)
                + ". It may only be blank where IW_FINAL_CONC_NM is filled."
            )
    return out


def _apply_schema_source_columns(
    frame: pd.DataFrame,
    source: pd.DataFrame,
    schema_columns: Sequence[str],
    *,
    group_col: str | None = None,
    explicit_inputs: Mapping[str, Any] | None = None,
    ignore_columns: Sequence[str] = (),
) -> pd.DataFrame:
    """Populate schema columns from explicit inputs, then exact/name-matched source columns.

    Explicit inputs always win. Source-name matching is only used when a destination
    column was not explicitly specified. For plate rows, the first non-null value per
    group is used; for well rows, values are copied row-for-row.
    """
    out = frame.copy()
    explicit = dict(explicit_inputs or {})
    ignored = set(ignore_columns)

    grouped_cache: dict[str, pd.Series] = {}
    if group_col is not None:
        validate_required_columns(source, [group_col], "source dataframe")

    def source_values(col: str) -> pd.Series:
        if group_col is None:
            return source[col].reset_index(drop=True)
        if col not in grouped_cache:
            grouped_cache[col] = (
                source[[group_col, col]]
                .dropna(subset=[col])
                .groupby(group_col, as_index=True)[col]
                .first()
            )
        return out[group_col].map(grouped_cache[col])

    for target in schema_columns:
        if target in ignored or target in explicit:
            continue
        source_col = _find_matching_source_column(source, target)
        if source_col is not None:
            out[target] = source_values(source_col)

    for target, value in explicit.items():
        if target in ignored:
            continue
        resolved = _resolve_value_from_source(
            source,
            value,
            group_key=out[group_col] if group_col is not None else None,
            label=target,
        )
        if isinstance(resolved, pd.Series):
            if group_col is not None and len(resolved) != len(out):
                # Source-column values may be indexed by source rows; map them per group
                # when the caller used a source column.
                source_col = _find_matching_source_column(source, value) if isinstance(value, str) else None
                if source_col is None:
                    raise ValueError(f"Resolved {target!r} has wrong length: {len(resolved)} != {len(out)}")
                out[target] = source_values(source_col).to_numpy()
            elif group_col is not None and resolved.index.equals(out.index):
                out[target] = resolved.to_numpy()
            else:
                out[target] = resolved.reset_index(drop=True).to_numpy()
        else:
            out[target] = resolved
    return out


def build_plate_rows(
    source: pd.DataFrame,
    *,
    plate_inputs: Mapping[str, Any],
    plate_barcode_suffix: str | None = None,
    engine: Engine | None = None,
    schema_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Build plate rows using one shared dynamic column-resolution path.

    Every exact/name-matched ``IP_*`` column already present in ``source`` is carried
    automatically. Explicit ``plate_inputs`` entries override that auto-match.
    ``IP_BARCODE`` is special: when not explicitly supplied (or supplied as ``None``),
    it defaults to the resolved ``IP_GRID_NUM`` and then receives the optional suffix.
    ``IP_INS_DATE`` is generated when not supplied. ``IP_ICO_ID``/``IP_COLLECTION``
    are resolved as before.
    """
    grid_src = _resolve_source_column(plate_inputs, "IP_GRID_NUM", source)
    source_work = source.copy()
    source_work["__RESOLVE_GRID__"] = source[grid_src]
    base = (
        source[[grid_src]]
        .drop_duplicates()
        .rename(columns={grid_src: "IP_GRID_NUM"})
        .reset_index(drop=True)
    )
    base["__RESOLVE_GRID__"] = base["IP_GRID_NUM"]
    schema = list(schema_columns or [c for c in source.columns if str(c).startswith("IP_")])
    if "IP_GRID_NUM" not in schema:
        schema.insert(0, "IP_GRID_NUM")

    # All automatically matching IP_* columns plus explicit inputs, excluding barcode
    # controls which have custom behavior below.
    controls = {"IP_GRID_NUM", "IP_BARCODE", "IP_BARCODE_SUFFIX", "IP_INS_DATE"}
    plates = _apply_schema_source_columns(
        base,
        source_work,
        schema,
        group_col="__RESOLVE_GRID__",
        explicit_inputs=plate_inputs,
        ignore_columns=controls,
    )
    plates["IP_GRID_NUM"] = base["IP_GRID_NUM"]

    barcode_input = plate_inputs.get("IP_BARCODE")
    if barcode_input is None:
        plates["IP_BARCODE"] = plates["IP_GRID_NUM"]
    else:
        resolved = _resolve_value_from_source(
            source_work, barcode_input, group_key=plates["IP_GRID_NUM"], label="IP_BARCODE"
        )
        if isinstance(resolved, pd.Series):
            source_col = _find_matching_source_column(source, barcode_input) if isinstance(barcode_input, str) else None
            if source_col is None:
                plates["IP_BARCODE"] = resolved.reset_index(drop=True).to_numpy()
            else:
                grouped = (
                    source_work[["__RESOLVE_GRID__", source_col]].dropna(subset=[source_col])
                    .groupby("__RESOLVE_GRID__", as_index=True)[source_col].first()
                )
                plates["IP_BARCODE"] = plates["IP_GRID_NUM"].map(grouped)
        else:
            plates["IP_BARCODE"] = resolved

    suffix = plate_inputs.get("IP_BARCODE_SUFFIX", plate_barcode_suffix)
    if suffix is not None and suffix != "":
        resolved_suffix = _resolve_value_from_source(
            source_work, suffix, group_key=plates["IP_GRID_NUM"], label="IP_BARCODE_SUFFIX"
        )
        if isinstance(resolved_suffix, pd.Series):
            source_col = _find_matching_source_column(source, suffix) if isinstance(suffix, str) else None
            if source_col is not None:
                grouped = (
                    source_work[["__RESOLVE_GRID__", source_col]].dropna(subset=[source_col])
                    .groupby("__RESOLVE_GRID__", as_index=True)[source_col].first()
                )
                plates["IP_BARCODE"] = plates["IP_BARCODE"].astype(str) + plates["IP_GRID_NUM"].map(grouped).astype(str)
            else:
                plates["IP_BARCODE"] = plates["IP_BARCODE"].astype(str) + resolved_suffix.astype(str).to_numpy()
        else:
            plates["IP_BARCODE"] = plates["IP_BARCODE"].astype(str) + str(resolved_suffix)

    if "IP_INS_DATE" in plate_inputs and plate_inputs["IP_INS_DATE"] is not None:
        resolved_date = _resolve_value_from_source(
            source_work, plate_inputs["IP_INS_DATE"], group_key=plates["IP_GRID_NUM"], label="IP_INS_DATE"
        )
        if isinstance(resolved_date, pd.Series):
            source_col = _find_matching_source_column(source, plate_inputs["IP_INS_DATE"]) if isinstance(plate_inputs["IP_INS_DATE"], str) else None
            if source_col is not None:
                grouped = (
                    source[[grid_src, source_col]].dropna(subset=[source_col])
                    .groupby(grid_src, as_index=True)[source_col].first()
                )
                plates["IP_INS_DATE"] = plates["IP_GRID_NUM"].map(grouped)
            else:
                plates["IP_INS_DATE"] = resolved_date.reset_index(drop=True).to_numpy()
        else:
            plates["IP_INS_DATE"] = resolved_date
    elif "IP_INS_DATE" not in plates.columns or plates["IP_INS_DATE"].isna().all():
        plates["IP_INS_DATE"] = now_string()

    plates = _resolve_collection(plates, engine=engine)
    return plates.drop(columns=["__RESOLVE_GRID__"], errors="ignore")

def build_well_rows(
    source: pd.DataFrame,
    *,
    plate_source_col: str,
    well_source_col: str,
    plate_lookup: Mapping[Any, Any] | pd.Series | pd.DataFrame | None,
    well_defaults: Mapping[str, Any],
    size: int,
    require_plate_ids: bool = True,
    well_barcode_source_col: str | None = None,
    plate_lookup_key_col: str | None = None,
    plate_lookup_value_col: str | None = None,
) -> pd.DataFrame:
    """Build well rows through the same dynamic source/name resolver as all other paths."""
    validate_required_columns(source, [plate_source_col, well_source_col], "source dataframe")
    frame = source.copy()
    schema = [c for c in source.columns if str(c).startswith("IW_")]
    frame = _apply_schema_source_columns(
        frame,
        source,
        schema,
        explicit_inputs=well_defaults,
        ignore_columns={"IW_IP_ID", "IW_COORDINATES", "IW_WELL_INDEX", "IW_BARCODE"},
    )
    lookup = to_mapping(plate_lookup, key_col=plate_lookup_key_col, value_col=plate_lookup_value_col)
    frame["IW_IP_ID"] = frame[plate_source_col].map(lookup) if lookup else pd.NA
    frame["IW_COORDINATES"] = frame[well_source_col].astype(str)
    frame["IW_WELL_INDEX"] = frame["IW_COORDINATES"].map(well_index_map(size))
    if require_plate_ids and frame["IW_IP_ID"].isna().any():
        missing = frame.loc[frame["IW_IP_ID"].isna(), plate_source_col].drop_duplicates().tolist()
        raise ValueError(f"Missing plate lookup values for: {missing}")
    if frame["IW_WELL_INDEX"].isna().any():
        missing = frame.loc[frame["IW_WELL_INDEX"].isna(), well_source_col].drop_duplicates().tolist()
        raise ValueError(f"Unrecognized well coordinates for size {size}: {missing}")
    if well_barcode_source_col is not None:
        src_col = _find_matching_source_column(source, well_barcode_source_col)
        if src_col is None:
            raise ValueError(f"IW_BARCODE source column {well_barcode_source_col!r} was not found")
        frame["IW_BARCODE"] = source[src_col].to_numpy()
    elif "IW_BARCODE" in well_defaults and well_defaults["IW_BARCODE"] is not None:
        value = well_defaults["IW_BARCODE"]
        resolved = _resolve_value_from_source(source, value, label="IW_BARCODE")
        frame["IW_BARCODE"] = resolved.to_numpy() if isinstance(resolved, pd.Series) else resolved
    return frame


def drop_empty_wells(wells: pd.DataFrame, compound_col: str = "IW_ICL_ID") -> pd.DataFrame:
    """Drop wells with no compound (``compound_col`` null): empty wells aren't inserted.

    A frame without ``compound_col`` is returned unchanged, since it can't say
    which wells are empty.
    """
    if compound_col not in wells.columns:
        return wells
    return wells[wells[compound_col].notna()].reset_index(drop=True)


def fill_missing_wells(
    wells: pd.DataFrame,
    *,
    plate_col: str,
    well_col: str,
    size: int,
    carry_columns: Sequence[str] = ("IW_IP_ID", "IW_INS_BY", "IW_INS_DATE", "IW_CONTAM"),
) -> pd.DataFrame:
    """Pad each plate to full geometry by adding missing coordinate rows.

    Note:
        ``IW_BARCODE`` is not carried: it holds a scanned matrix-tube barcode, which
        an empty padded well doesn't have. Compound/lot fields such as ``IW_ICL_ID``
        are excluded too, so empty wells do not inherit compound identity.

    Parameters:
        wells: Well dataframe containing at least plate and well columns.
        plate_col: Column identifying destination plate groups.
        well_col: Column containing destination well coordinates.
        size: Plate size used to determine full coordinate set.
        carry_columns: Columns copied from the first row of each plate into new rows.
    """
    validate_required_columns(wells, [plate_col, well_col], "wells dataframe")
    complete = well_coordinates(size)
    missing_rows: list[pd.DataFrame] = []
    for plate_value, group in wells.groupby(plate_col, dropna=False):
        missing = sorted(set(complete) - set(group[well_col].dropna().astype(str)))
        if not missing:
            continue
        template = {}
        if not group.empty:
            first = group.iloc[0]
            for column in carry_columns:
                if column in group.columns:
                    template[column] = first[column]
        rows = []
        for well_value in missing:
            row = {column: pd.NA for column in wells.columns}
            row[plate_col] = plate_value
            row[well_col] = well_value
            if "IW_COORDINATES" in row:
                row["IW_COORDINATES"] = well_value
            row["IW_WELL_INDEX"] = well_index_map(size)[well_value]
            for column, value in template.items():
                row[column] = value
            rows.append(row)
        missing_rows.append(pd.DataFrame(rows))
    if not missing_rows:
        return wells.copy()
    return pd.concat([wells, *missing_rows], ignore_index=True)


def align_to_schema(
    frame: pd.DataFrame,
    schema_columns: Sequence[str],
    *,
    fill_value: Any = pd.NA,
) -> pd.DataFrame:
    """Add missing schema columns and return schema-first ordered dataframe.

    Parameters:
        frame: Input dataframe.
        schema_columns: Ordered schema column names.
        fill_value: Value used for missing columns.
    """
    aligned = frame.copy()
    for column in schema_columns:
        if column not in aligned.columns:
            aligned[column] = fill_value
    extra = [column for column in aligned.columns if column not in schema_columns]
    ordered = aligned[[column for column in schema_columns if column in aligned.columns]]
    if extra:
        ordered = pd.concat([ordered, aligned[extra]], axis=1)
    return ordered


def coerce_to_schema(frame: pd.DataFrame, schema_columns: Sequence[str], *, fill_value: Any = pd.NA) -> pd.DataFrame:
    """Return dataframe restricted to schema columns in schema order.

    Parameters:
        frame: Input dataframe.
        schema_columns: Ordered schema column names.
        fill_value: Value used for missing columns.
    """
    aligned = frame.copy()
    for column in schema_columns:
        if column not in aligned.columns:
            aligned[column] = fill_value
    return aligned[list(schema_columns)].copy()


def prepare_insert_frames(
    plates: pd.DataFrame,
    wells: pd.DataFrame,
    *,
    plate_table: Table,
    well_table: Table,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate and coerce plate/well frames to reflected table schemas.

    Parameters:
        plates: Candidate INV_PLATE rows.
        wells: Candidate INV_WELL rows.
        plate_table: Reflected INV_PLATE table.
        well_table: Reflected INV_WELL table.

    Notes:
        This helper expects wells to already contain resolved IW_IP_ID values.
        If IP_ID values are auto-generated on insert, use the staged helpers:
        prepare_plate_insert_frame -> insert plates -> fetch_inserted_plate_ids
        -> attach_plate_ids_to_wells -> prepare_well_insert_frame.
    """
    plate_columns = table_column_names(plate_table)
    well_columns = table_column_names(well_table)

    validate_required_columns(
        plates,
        ["IP_GRID_NUM", "IP_BARCODE", "IP_IPT_ID", "IP_STATUS_ID", "IP_SCREENING_CENTER"],
        "plate insert frame",
    )
    validate_required_columns(
        wells,
        ["IW_IP_ID", "IW_COORDINATES", "IW_WELL_INDEX"],
        "well insert frame",
    )

    validate_unique(plates, ["IP_GRID_NUM", "IP_BARCODE"], "plate insert frame")
    validate_unique(wells, ["IW_IP_ID", "IW_COORDINATES"], "well insert frame")

    return (
        coerce_to_schema(plates, plate_columns),
        coerce_to_schema(wells, well_columns),
    )


def prepare_plate_insert_frame(plates: pd.DataFrame, *, plate_table: Table) -> pd.DataFrame:
    """Validate and coerce only plate rows for first-phase insert.

    Parameters:
        plates: Candidate INV_PLATE rows.
        plate_table: Reflected INV_PLATE table.
    """
    plate_columns = table_column_names(plate_table)
    validate_required_columns(
        plates,
        ["IP_GRID_NUM", "IP_BARCODE", "IP_IPT_ID", "IP_STATUS_ID", "IP_SCREENING_CENTER"],
        "plate insert frame",
    )
    validate_unique(plates, ["IP_GRID_NUM", "IP_BARCODE"], "plate insert frame")
    return coerce_to_schema(plates, plate_columns)


def fetch_inserted_plate_ids(
    engine: Engine,
    *,
    plate_table: Table,
    key_values: Sequence[Any],
    key_col: str = "IP_GRID_NUM",
    id_col: str = "IP_ID",
    include_cols: Sequence[str] = ("IP_BARCODE",),
) -> pd.DataFrame:
    """Fetch inserted plate IDs using plate keys after plate insert.

    Parameters:
        engine: SQLAlchemy engine.
        plate_table: Reflected INV_PLATE table.
        key_values: Plate key values used to query inserted rows.
        key_col: Plate table key column used for lookup.
        id_col: Plate ID column to return.
        include_cols: Additional plate columns to include in output.
    """
    if not key_values:
        raise ValueError("key_values is empty; nothing to fetch")
    requested_cols = [key_col, id_col, *include_cols]
    stmt_cols = [plate_table.c[col] for col in requested_cols]
    stmt = select(*stmt_cols).where(plate_table.c[key_col].in_(list(key_values)))
    with engine.connect() as conn:
        rows = conn.execute(stmt).fetchall()
    out = pd.DataFrame(rows, columns=requested_cols)
    missing = sorted(set(key_values) - set(out[key_col].tolist()))
    if missing:
        raise ValueError(f"Missing inserted plates for {key_col}: {missing}")
    return out


def attach_plate_ids_to_wells(
    wells: pd.DataFrame,
    plate_ids: pd.DataFrame,
    *,
    well_key_col: str = "plate",
    plate_key_col: str = "IP_GRID_NUM",
    plate_id_col: str = "IP_ID",
    target_col: str = "IW_IP_ID",
    strict: bool = True,
) -> pd.DataFrame:
    """Attach generated plate IDs to well rows in the second phase.

    Parameters:
        wells: Well dataframe containing a plate label column.
        plate_ids: Dataframe containing plate key and generated plate IDs.
        well_key_col: Plate label column in wells.
        plate_key_col: Key column in plate_ids.
        plate_id_col: ID column in plate_ids.
        target_col: Output ID column written on wells.
        strict: If True, raise when any wells remain unresolved.
    """
    validate_required_columns(wells, [well_key_col], "wells dataframe")
    validate_required_columns(plate_ids, [plate_key_col, plate_id_col], "plate id dataframe")
    mapping = dict(zip(plate_ids[plate_key_col], plate_ids[plate_id_col]))
    out = wells.copy()
    out[target_col] = out[well_key_col].map(mapping)
    if strict and out[target_col].isna().any():
        missing = out.loc[out[target_col].isna(), well_key_col].drop_duplicates().tolist()
        raise ValueError(f"Missing plate IDs for well rows: {missing}")
    return out


def prepare_well_insert_frame(wells: pd.DataFrame, *, well_table: Table) -> pd.DataFrame:
    """Validate and coerce only well rows for second-phase insert.

    Parameters:
        wells: Candidate INV_WELL rows with resolved IW_IP_ID.
        well_table: Reflected INV_WELL table.
    """
    well_columns = table_column_names(well_table)
    validate_required_columns(
        wells,
        ["IW_IP_ID", "IW_COORDINATES", "IW_WELL_INDEX"],
        "well insert frame",
    )
    validate_unique(wells, ["IW_IP_ID", "IW_COORDINATES"], "well insert frame")
    return coerce_to_schema(wells, well_columns)


def _assert_unique_columns(source: pd.DataFrame) -> None:
    """Raise a clear error if ``source`` has duplicate column labels.

    Concatenating several spreadsheets and normalizing headers (lowercasing,
    replacing spaces) can collapse distinct headers onto the same name.  Duplicate
    columns make ``source[name]`` return a DataFrame instead of a Series, which
    silently corrupts column-or-literal resolution and later breaks ``pd.concat`` in
    :func:`fill_missing_wells`, so we fail early with the offending names.
    """
    dupes = source.columns[source.columns.duplicated()].unique().tolist()
    if dupes:
        raise ValueError(
            f"source has duplicate column names: {dupes}. Deduplicate before "
            "registering, e.g. source = source.loc[:, ~source.columns.duplicated()]"
        )


def _is_source_column(source: pd.DataFrame, value: Any) -> bool:
    """Return True when ``value`` names an existing column in ``source``.

    Used to decide, for each ``IP_``/``IW_`` input, whether to pull a column from the
    source frame or broadcast the value as a literal.
    """
    return isinstance(value, str) and value in source.columns


def _resolve_source_column(
    inputs: Mapping[str, Any],
    key: str,
    source: pd.DataFrame,
    *,
    required: bool = True,
) -> str | None:
    """Resolve the source column name that ``inputs[key]`` must point to.

    ``IP_GRID_NUM`` and ``IW_COORDINATES`` are structural: they identify the plate
    group and the well coordinate, so their input values must name real source
    columns rather than literals.
    """
    if key not in inputs:
        if required:
            raise ValueError(f"{key!r} is required in inputs and must name a source column")
        return None
    value = inputs[key]
    if not _is_source_column(source, value):
        raise ValueError(
            f"{key!r} must map to an existing source column; got {value!r} "
            f"(source columns: {list(source.columns)})"
        )
    return value


def prepare_inventory_frames(
    source: pd.DataFrame,
    *,
    plate_inputs: Mapping[str, Any],
    well_inputs: Mapping[str, Any],
    size: int,
    plate_lookup: Mapping[Any, Any] | pd.Series | pd.DataFrame | None = None,
    engine: Engine | None = None,
    plate_barcode_suffix: str | None = None,
    require_plate_ids: bool = True,
    fill_missing: bool = False,
    plate_schema_columns: Sequence[str] | None = None,
    well_schema_columns: Sequence[str] | None = None,
    require_concentration: bool | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare NEW-plate inventory rows using the shared dynamic resolver.

    ``require_concentration`` controls the concentration/volume rule for wells that
    carry a compound: at least one of ``IW_CONC_UM`` / ``IW_FINAL_CONC_NM`` must be
    filled (both are allowed), and ``IW_VOL_UL`` may be blank only where
    ``IW_FINAL_CONC_NM`` is filled. ``None`` (default) enforces it whenever
    ``well_inputs`` supplies any of ``IW_VOL_UL`` / ``IW_CONC_UM`` /
    ``IW_FINAL_CONC_NM``; ``True`` always enforces it; ``False`` never does.

    Exact/name-matched ``IP_*``/``IW_*`` source columns are carried automatically;
    explicit inputs override those values. ``IW_BARCODE`` is the scanned
    matrix-tube barcode: it comes only from an ``IW_BARCODE`` source column (or an
    explicit ``well_inputs`` entry), is optional, and is never filled from
    ``IP_BARCODE``. Padded empty wells get no barcode.

    Empty wells are not inserted: rows with no ``IW_ICL_ID`` (no compound) are
    dropped, and plates are not padded to full geometry unless
    ``fill_missing=True`` is passed explicitly (Amanda, Nov 2026).
    """
    _assert_unique_columns(source)
    grid_src = _resolve_source_column(plate_inputs, "IP_GRID_NUM", source)
    well_src = _resolve_source_column(well_inputs, "IW_COORDINATES", source)
    insert_time = now_string()

    plate_inputs = dict(plate_inputs)
    well_inputs = dict(well_inputs)
    plate_inputs.setdefault("IP_INS_DATE", insert_time)
    well_inputs.setdefault("IW_INS_DATE", insert_time)
    numeric_keys = [k for k in NUMERIC_WELL_FIELDS if well_inputs.get(k) is not None]
    for key in numeric_keys:
        well_inputs[key] = _coerce_numeric_well_input(source, key, well_inputs[key])

    plates = build_plate_rows(
        source,
        plate_inputs=plate_inputs,
        plate_barcode_suffix=plate_barcode_suffix,
        engine=engine,
        schema_columns=plate_schema_columns,
    )

    well_base = source.copy()
    well_schema = list(well_schema_columns or [c for c in source.columns if str(c).startswith("IW_")])
    if "IW_COORDINATES" not in well_schema:
        well_schema.append("IW_COORDINATES")
    # IW_BARCODE (scanned tube barcode, optional) is resolved like any other field.
    wells = _apply_schema_source_columns(
        well_base,
        source,
        well_schema,
        explicit_inputs=well_inputs,
        ignore_columns={"IW_COORDINATES", "IW_WELL_INDEX"},
    )
    wells["IW_COORDINATES"] = source[well_src].astype(str).to_numpy()
    # Validate before fill_missing_wells so only real source rows are checked (and after
    # IW_COORDINATES is set so error messages can name the offending wells).
    wells = _validate_numeric_well_fields(
        wells,
        numeric_keys,
        enforce_concentration_rule=bool(numeric_keys) if require_concentration is None else require_concentration,
    )
    if "IW_WELL_INDEX" in well_inputs and well_inputs["IW_WELL_INDEX"] is not None:
        resolved_index = _resolve_value_from_source(source, well_inputs["IW_WELL_INDEX"], label="IW_WELL_INDEX")
        wells["IW_WELL_INDEX"] = resolved_index.to_numpy() if isinstance(resolved_index, pd.Series) else resolved_index
    else:
        wells["IW_WELL_INDEX"] = wells["IW_COORDINATES"].map(well_index_map(size))

    lookup = to_mapping(plate_lookup)
    if lookup:
        wells["IW_IP_ID"] = source[grid_src].map(lookup).to_numpy()
    elif "IW_IP_ID" in wells.columns and not wells["IW_IP_ID"].isna().all():
        pass
    else:
        wells["IW_IP_ID"] = pd.NA

    if wells["IW_WELL_INDEX"].isna().any():
        missing = wells.loc[wells["IW_WELL_INDEX"].isna(), "IW_COORDINATES"].drop_duplicates().tolist()
        raise ValueError(f"Unrecognized well coordinates for size {size}: {missing}")
    if require_plate_ids and wells["IW_IP_ID"].isna().any():
        missing = wells.loc[wells["IW_IP_ID"].isna(), grid_src].drop_duplicates().tolist()
        raise ValueError(f"Missing plate lookup values for: {missing}")

    wells = drop_empty_wells(wells)
    if fill_missing:
        wells = fill_missing_wells(wells, plate_col=grid_src, well_col="IW_COORDINATES", size=size)

    if "IW_BARCODE" not in wells.columns:
        wells["IW_BARCODE"] = None  # no scanned tube barcodes in this input

    return plates, wells

def remap_wells_by_quadrant(
    mother_wells: pd.DataFrame,
    *,
    dest_size: int,
    destination_plate_names: Sequence[str],
    source_plate_col: str,
    source_well_col: str = "IW_COORDINATES",
    layout_map: pd.DataFrame | str | Path | None = None,
    profile: str | None = None,
    custom_path: str | Path | None = None,
    source_plate_order: Sequence[Any] | Mapping[str, Sequence[Any]] | None = None,
    plates_per_destination: int = 4,
    dest_plate_col: str = "_DEST_PLATE_",
    keep_source_in_notes: bool = True,
    overrides: Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Reshape DB-sourced mother wells into daughter plate coordinates by quadrant.

    Unified engine for the ``interleave`` (e.g. 4×96→384), ``quadrant``
    (e.g. 4×384→1536) and ``custom`` layouts.  Ordered parent plates are assigned to
    quadrant slots ``Q1..Qn`` in chunks of ``plates_per_destination``; the wells are
    merged with the normalized map on ``(source_well, source_quadrant)`` and rewritten
    to destination coordinates.  The **final destination may be partially filled**
    (e.g. 13 parents → ``[4, 4, 4, 1]``); unfilled quadrants produce no rows and are
    padded as empty wells downstream via :func:`fill_missing_wells`.

    Rows whose source well has no map entry (e.g. interleave border wells) are dropped.

    Parameters:
        mother_wells: Source well rows fetched from the DB (carry ``IW_ID``, ``IW_IP_ID``, etc.).
        dest_size: Destination plate size (384, 1536, ...) used for ``IW_WELL_INDEX``.
        destination_plate_names: Ordered daughter plate labels; written into ``dest_plate_col``.
            There should be ``ceil(n_parents / plates_per_destination)`` of them.
        source_plate_col: Column in ``mother_wells`` identifying the parent plate
            (e.g. ``IP_GRID_NUM`` or ``IW_IP_ID``).
        source_well_col: Column holding the source well coordinate.
        layout_map: Normalized map dataframe/path (mutually exclusive with ``profile``).
        profile: Named profile (``interleave``/``quadrant``/``custom``) resolved via
            :func:`resolve_layout_map` when ``layout_map`` is not given.
        custom_path: Custom map path when ``profile='custom'``.
        source_plate_order: Explicit parent order — a flat list (chunked sequentially
            into destinations) or a ``{dest_name: [parents]}`` mapping.  Defaults to
            natural-sorted unique parents.
        plates_per_destination: Parents merged into each destination (must equal the
            number of quadrant labels in the map).
        dest_plate_col: Output column naming the destination daughter.
        keep_source_in_notes: Write ``IW_NOTES = "<parent>-<well>"`` for provenance.
        overrides: Column→value overrides applied last to the remapped well rows.
    """
    validate_required_columns(mother_wells, [source_plate_col, source_well_col], "mother wells")
    if profile is not None:
        mapping = resolve_layout_map(profile=profile, custom_path=custom_path)
    elif layout_map is not None:
        mapping = load_layout_map(layout_map)
    else:
        raise ValueError("Provide either profile= or layout_map=")
    validate_required_columns(mapping, list(MAP_STANDARD_COLS), "layout map")

    quadrants = sorted(
        mapping["source_quadrant"].dropna().astype(str).unique().tolist(),
        key=_natural_sort_key,
    )

    assignment = _build_quadrant_assignment(
        mother_wells,
        source_plate_col=source_plate_col,
        destination_plate_names=list(destination_plate_names),
        source_plate_order=source_plate_order,
        plates_per_destination=plates_per_destination,
        quadrants=quadrants,
    )
    assignment = assignment.rename(columns={source_plate_col: "__SRC_PLATE_KEY__"})

    left = mother_wells.copy()
    left["__SRC_PLATE_KEY__"] = left[source_plate_col].astype(str)
    merged = left.merge(assignment, how="inner", on="__SRC_PLATE_KEY__")
    merged = merged.merge(
        mapping,
        how="left",
        left_on=[source_well_col, "_SRC_QUADRANT_"],
        right_on=["source_well", "source_quadrant"],
    )

    # Drop source wells with no map entry (e.g. interleave border columns).
    merged = merged.loc[merged["dest_well"].notna()].copy()

    if keep_source_in_notes:
        merged["IW_NOTES"] = merged[source_plate_col].astype(str) + "-" + merged[source_well_col].astype(str)

    merged["IW_COORDINATES"] = merged["dest_well"]
    merged["IW_WELL_INDEX"] = merged["dest_well"].map(well_index_map(dest_size))
    merged[dest_plate_col] = merged["_DEST_PLATE_"].astype(str)

    drop_cols = [
        "source_well", "source_quadrant", "dest_well",
        "_SRC_QUADRANT_", "_DEST_PLATE_", "__SRC_PLATE_KEY__",
    ]
    drop_cols = [c for c in drop_cols if c in merged.columns and c != dest_plate_col]
    result = merged.drop(columns=drop_cols).reset_index(drop=True)

    # Explicit overrides always win over remapped/source-derived values. Mapping values
    # are keyed by destination plate, while source-column names remain available for
    # callers that want to carry a dynamically-fetched parent field through the remap.
    for column, value in (overrides or {}).items():
        if isinstance(value, Mapping):
            keyed = {str(k): v for k, v in value.items()}
            result[column] = result[dest_plate_col].astype(str).map(keyed)
        elif _is_source_column(mother_wells, value):
            # Align by original row order only when lengths match; otherwise fail loudly
            # instead of silently assigning incorrect data after a many-to-one remap.
            if len(mother_wells) != len(result):
                raise ValueError(
                    f"Override {column!r} references source column {value!r}, but remapping changed row count "
                    f"from {len(mother_wells)} to {len(result)}"
                )
            result[column] = mother_wells[value].values
        else:
            result[column] = value
    return result


def remap_384_to_1536(
    wells: pd.DataFrame,
    *,
    layout_map: str | Path | pd.DataFrame | None = None,
    layout_profile: str = "quadrant",
    custom_layout_path: str | Path | None = None,
    dest_plate_col: str = "IW_IP_ID",
    source_plate_col: str,
    source_well_col: str = "IW_COORDINATES",
    source_plate_order: Mapping[Any, Sequence[str]] | Sequence[str] | None = None,
    keep_original_in_notes: bool = True,
) -> pd.DataFrame:
    """Back-compat wrapper: remap grouped 384-well rows into 1536 coordinates.

    ``source_plate_col`` names the column that identifies each parent plate (for
    example ``IP_GRID_NUM`` or ``IP_BARCODE`` joined from INV_PLATE). It has no
    default: it used to default to ``IW_BARCODE``, which is a per-well tube
    barcode, not a plate identifier, and is empty on most plates.

    Thin adapter over :func:`remap_wells_by_quadrant`.  Input is one dataframe of
    source 384-well rows where ``dest_plate_col`` already identifies the destination
    1536 plate grouping; parents within each group are assigned to quadrant slots in
    ``source_plate_order`` (or natural-sorted).  Unlike the historical implementation,
    unmapped wells are dropped rather than raising (the quadrant map is complete, so
    this is a no-op for 384→1536), and a partially-filled destination is allowed.

    Parameters mirror the original function; see :func:`remap_wells_by_quadrant`.
    """
    validate_required_columns(
        wells, [dest_plate_col, source_plate_col, source_well_col], "wells dataframe"
    )
    if layout_map is not None:
        mapping = load_layout_map(layout_map)
    else:
        mapping = resolve_layout_map(profile=layout_profile, custom_path=custom_layout_path)
    validate_required_columns(mapping, list(MAP_STANDARD_COLS), "1536 layout map")
    n_quadrants = mapping["source_quadrant"].dropna().astype(str).nunique()

    dests = list(dict.fromkeys(wells[dest_plate_col].dropna().tolist()))
    order_map: dict[str, list[str]] = {}
    for dest in dests:
        group = wells[wells[dest_plate_col] == dest]
        if isinstance(source_plate_order, Mapping):
            plates = list(source_plate_order.get(dest, []))
        elif source_plate_order is not None:
            plates = list(source_plate_order)
        else:
            plates = []
        if not plates:
            plates = sorted(
                group[source_plate_col].dropna().astype(str).unique().tolist(),
                key=_natural_sort_key,
            )
        order_map[str(dest)] = [str(p) for p in plates]

    return remap_wells_by_quadrant(
        wells,
        dest_size=1536,
        destination_plate_names=[str(d) for d in dests],
        source_plate_col=source_plate_col,
        source_well_col=source_well_col,
        layout_map=mapping,
        profile=None,
        source_plate_order=order_map,
        plates_per_destination=n_quadrants,
        dest_plate_col=dest_plate_col,
        keep_source_in_notes=keep_original_in_notes,
    )


def stage_remapped_wells(
    remapped_wells: pd.DataFrame,
    *,
    destination_plate_ids: pd.DataFrame,
    parent_to_daughter: Mapping[Any, Any],
    dest_plate_col: str = "_DEST_PLATE_",
    well_inputs: Mapping[str, Any] | None = None,
    well_overrides: Mapping[str, Any] | None = None,
    ins_by: int | None = None,
    ins_date: str | None = None,
) -> pd.DataFrame:
    """Finalize an interleave/quadrant/custom remap through the shared daughter-well core.

    The remap itself only changes coordinates/destination grouping. This helper then
    stages the real daughter ``IW_IP_ID`` while preserving every other source
    ``IW_*`` field, especially ``IW_ICL_ID``. ``IW_BARCODE`` is left empty (see
    stage_daughter_wells).
    """
    validate_required_columns(remapped_wells, ["IW_ID", "IW_IP_ID", dest_plate_col], "remapped wells")
    validate_required_columns(destination_plate_ids, ["IP_GRID_NUM", "IP_ID"], "destination plate dataframe")
    destination_to_id = dict(zip(destination_plate_ids["IP_GRID_NUM"].astype(str), destination_plate_ids["IP_ID"]))
    translated_overrides: dict[str, Any] = {}
    for column, value in (well_overrides or {}).items():
        if isinstance(value, Mapping):
            translated_overrides[column] = {destination_to_id.get(str(k), k): v for k, v in value.items()}
        else:
            translated_overrides[column] = value
    result = stage_daughter_wells(
        remapped_wells,
        new_ip_id_map=parent_to_daughter,
        well_inputs=well_inputs,
        well_overrides=translated_overrides,
        ins_by=ins_by,
        ins_date=ins_date,
    )
    result = result.drop(columns=[dest_plate_col], errors="ignore")
    return result


def attach_source_well_ids(
    wells: pd.DataFrame,
    *,
    source_id_col: str = "IW_ID",
    target_col: str = "IW_SOURCE_ID",
    drop_source_id: bool = False,
) -> pd.DataFrame:
    """Assign source well IDs directly from an existing ID column.

    Parameters:
        wells: Input well dataframe.
        source_id_col: Column containing existing source well IDs.
        target_col: Output column to receive source well IDs.
        drop_source_id: Whether to remove source_id_col after assignment.
    """
    validate_required_columns(wells, [source_id_col], "wells dataframe")
    frame = wells.copy()
    if frame[source_id_col].isna().any():
        raise ValueError(f"{source_id_col} contains null values; cannot assign {target_col}")
    frame[target_col] = frame[source_id_col]
    if drop_source_id:
        frame = frame.drop(columns=[source_id_col])
    return frame


def resolve_source_well_ids_by_lot(
    wells: pd.DataFrame,
    source_wells: pd.DataFrame,
    *,
    lot_col: str = "IW_ICL_ID",
    source_lot_col: str = "IW_ICL_ID",
    source_id_col: str = "IW_ID",
    source_plate_col: str = "IW_IP_ID",
    allowed_source_plate_ids: Iterable[Any] | None = None,
    target_col: str = "IW_SOURCE_ID",
) -> pd.DataFrame:
    """Resolve daughter-plate source well IDs from lot IDs using strict matching.

    Parameters:
        wells: Daughter well dataframe that needs IW_SOURCE_ID values.
        source_wells: Source well pool used for lot-to-well lookup.
        lot_col: Lot ID column in daughter wells.
        source_lot_col: Lot ID column in source well pool.
        source_id_col: Source well ID column in source well pool.
        source_plate_col: Source plate column used for optional filtering.
        allowed_source_plate_ids: Optional source plate IDs to restrict lookup scope.
        target_col: Output column to receive resolved source well IDs.
    """
    validate_required_columns(wells, [lot_col], "wells dataframe")
    validate_required_columns(source_wells, [source_lot_col, source_id_col], "source wells dataframe")
    pool = source_wells.copy()
    if allowed_source_plate_ids is not None:
        validate_required_columns(pool, [source_plate_col], "source wells dataframe")
        pool = pool[pool[source_plate_col].isin(list(allowed_source_plate_ids))]

    candidates = pool[[source_lot_col, source_id_col]].dropna().drop_duplicates()
    dup_lots = candidates[candidates.duplicated(subset=[source_lot_col], keep=False)]
    if not dup_lots.empty:
        ambiguous = dup_lots.sort_values([source_lot_col, source_id_col]).to_dict(orient="records")
        raise ValueError(f"Ambiguous source well IDs per lot in source pool: {ambiguous}")

    lot_to_source = dict(zip(candidates[source_lot_col], candidates[source_id_col]))
    frame = wells.copy()
    frame[target_col] = frame[lot_col].map(lot_to_source)
    unresolved = frame.loc[frame[lot_col].notna() & frame[target_col].isna(), [lot_col]].drop_duplicates()
    if not unresolved.empty:
        raise ValueError(f"Could not resolve source wells for lots: {unresolved[lot_col].tolist()}")
    return frame


def fetch_collection_map(
    engine: Engine,
    collection_table: Table | None = None,
    *,
    id_col: str = "ICO_ID",
    name_col: str = "ICO_NAME",
) -> pd.DataFrame:
    """Fetch ``(ICO_ID, ICO_NAME)`` rows from INV_COLLECTION.

    Used internally by :func:`_resolve_collection` so that supplying only one of
    ``IP_ICO_ID`` (``ICO_ID``) / ``IP_COLLECTION`` (``ICO_NAME``) resolves the other.

    Parameters:
        engine: SQLAlchemy engine.
        collection_table: Reflected INV_COLLECTION table; reflected from ``engine`` by
            name (``INV_COLLECTION``) when omitted.
        id_col: Collection id column (mapped to ``IP_ICO_ID``).
        name_col: Collection name column (mapped to ``IP_COLLECTION``).
    """
    if collection_table is None:
        collection_table = reflect_table(engine, COLLECTION_TABLE_NAME)
    stmt = select(collection_table.c[id_col], collection_table.c[name_col])
    with engine.connect() as conn:
        rows = conn.execute(stmt).fetchall()
    return pd.DataFrame(rows, columns=[id_col, name_col])


def fetch_plates(
    engine: Engine,
    plate_table: Table,
    keys: Sequence[Any],
    *,
    key_col: str = "IP_GRID_NUM",
    ipt_id: int | None = None,
    barcode_contains: str | None = None,
    use_like: bool | None = None,
) -> pd.DataFrame:
    """Fetch source plate rows (including ``IP_ID``) from the DB for daughter creation.

    Mirrors the notebooks' grid-num / IPT / barcode-substring filtering, and also
    satisfies the "look up parent IP_IDs from a barcode list on a later day" case
    (pass ``key_col='IP_BARCODE'`` with ``use_like=False``).

    Parameters:
        engine: SQLAlchemy engine.
        plate_table: Reflected INV_PLATE table.
        keys: Values to match against ``key_col``.
        key_col: Column matched against ``keys`` (``IP_GRID_NUM`` or ``IP_BARCODE``).
        ipt_id: Optional ``IP_IPT_ID`` filter (e.g. 1=master, 9=acoustic, 11=...).
        barcode_contains: Optional post-filter keeping rows whose ``IP_BARCODE``
            contains this substring (e.g. ``'MA01'``, ``'Eco1'``).
        use_like: Use SQL ``LIKE`` per key (default True for ``IP_GRID_NUM``) instead
            of ``IN`` (default for other columns).
    """
    keys = list(keys)
    if not keys:
        raise ValueError("keys is empty; nothing to fetch")
    if use_like is None:
        use_like = key_col == "IP_GRID_NUM"
    col = plate_table.c[key_col]
    key_filter = or_(*[col.like(k) for k in keys]) if use_like else col.in_(keys)
    if ipt_id is not None:
        stmt = select(plate_table).where(and_(plate_table.c.IP_IPT_ID == ipt_id, key_filter))
    else:
        stmt = select(plate_table).where(key_filter)
    df = pd.read_sql(stmt, con=engine)
    if barcode_contains:
        df = df[df["IP_BARCODE"].str.contains(barcode_contains, na=False)].copy()
    if df.empty:
        raise ValueError(
            f"No plates matched {key_col} in {keys} "
            f"(ipt_id={ipt_id}, barcode_contains={barcode_contains})"
        )
    return df.reset_index(drop=True)


def fetch_mother_wells(
    engine: Engine,
    well_table: Table,
    parent_ip_ids: Sequence[Any],
    *,
    plate_table: Table | None = None,
    plate_cols: Sequence[str] = ("IP_GRID_NUM",),
) -> pd.DataFrame:
    """Fetch full source well rows (including ``IW_ID``) for the given parent plates.

    Parameters:
        engine: SQLAlchemy engine.
        well_table: Reflected INV_WELL table.
        parent_ip_ids: Parent ``IP_ID`` values whose wells to fetch.
        plate_table: Optional reflected INV_PLATE table; when given, ``plate_cols``
            are joined in (e.g. ``IP_GRID_NUM``) so callers can order/label parents.
        plate_cols: Plate columns to attach when ``plate_table`` is provided.
    """
    parent_ip_ids = list(parent_ip_ids)
    if not parent_ip_ids:
        raise ValueError("parent_ip_ids is empty; nothing to fetch")
    stmt = select(well_table).where(well_table.c.IW_IP_ID.in_(parent_ip_ids))
    wells = pd.read_sql(stmt, con=engine)
    if wells.empty:
        raise ValueError(f"No wells found for parent plate ids: {parent_ip_ids}")
    if plate_table is not None:
        cols = ["IP_ID", *plate_cols]
        pstmt = select(*[plate_table.c[c] for c in cols]).where(plate_table.c.IP_ID.in_(parent_ip_ids))
        plates = pd.read_sql(pstmt, con=engine)
        wells = wells.merge(plates, how="left", left_on="IW_IP_ID", right_on="IP_ID").drop(
            columns=["IP_ID"], errors="ignore"
        )
    return wells.reset_index(drop=True)


def assert_new_barcodes(engine: Engine, plate_table: Table, barcodes: Sequence[Any]) -> None:
    """Raise if any of the given plate barcodes already exist in the DB.

    Parameters:
        engine: SQLAlchemy engine.
        plate_table: Reflected INV_PLATE table.
        barcodes: Candidate daughter plate barcodes.
    """
    barcodes = list(barcodes)
    if not barcodes:
        raise ValueError("barcodes is empty; nothing to check")
    stmt = select(plate_table.c.IP_BARCODE).where(plate_table.c.IP_BARCODE.in_(barcodes))
    with engine.connect() as conn:
        existing = [row[0] for row in conn.execute(stmt)]
    if existing:
        raise ValueError(f"Barcodes already exist in {plate_table.name}: {sorted(set(existing))}")


def assert_no_existing_wells(engine: Engine, well_table: Table, daughter_ip_ids: Sequence[Any]) -> None:
    """Raise if any wells already exist for the given daughter plate IDs.

    Parameters:
        engine: SQLAlchemy engine.
        well_table: Reflected INV_WELL table.
        daughter_ip_ids: Newly-inserted daughter plate IP_IDs.
    """
    ids = list(daughter_ip_ids)
    if not ids:
        raise ValueError("daughter_ip_ids is empty; nothing to check")
    stmt = select(well_table.c.IW_IP_ID).where(well_table.c.IW_IP_ID.in_(ids))
    with engine.connect() as conn:
        existing = sorted({row[0] for row in conn.execute(stmt)})
    if existing:
        raise ValueError(f"Wells already exist for plate ids: {existing}")


def copy_daughter_plates(
    source_plates: pd.DataFrame,
    *,
    barcode_suffix: str | None = None,
    grid_to_barcode: Mapping[str, str] | None = None,
    plate_inputs: Mapping[str, Any] | None = None,
    overrides: Mapping[str, Any] | None = None,
    ins_by: int | None = None,
    ins_date: str | None = None,
    id_col: str = "IP_ID",
    source_id_col: str = "IP_SOURCE_ID",
    grid_col: str = "IP_GRID_NUM",
    barcode_col: str = "IP_BARCODE",
) -> pd.DataFrame:
    """Build 1:1 daughter plates using the same dynamic resolver as new plates.

    Every source ``IP_*`` column is auto-carried by name match unless explicitly
    overridden. ``IP_ID`` becomes ``IP_SOURCE_ID``. ``IP_INS_DATE`` is refreshed for
    the daughter unless explicitly supplied. ``IP_BARCODE`` follows the shared special
    precedence: explicit input/map/suffix, otherwise daughter grid number. Overrides win.
    """
    validate_required_columns(source_plates, [id_col, grid_col], "source plates")
    inputs = dict(plate_inputs or {})
    if barcode_col not in inputs:
        # Keep copy's explicit barcode controls as inputs for the shared resolver.
        if grid_to_barcode is not None:
            inputs[barcode_col] = {str(k): v for k, v in grid_to_barcode.items()}
        else:
            inputs[barcode_col] = None
    if "IP_BARCODE_SUFFIX" not in inputs:
        inputs["IP_BARCODE_SUFFIX"] = barcode_suffix
    if "IP_INS_DATE" not in inputs:
        inputs["IP_INS_DATE"] = ins_date

    # Build from the parent source without the auto-carry of parent identity/insert date.
    schema = [c for c in source_plates.columns if str(c).startswith("IP_")]
    base = source_plates[[grid_col]].drop_duplicates().rename(columns={grid_col: "IP_GRID_NUM"}).reset_index(drop=True)
    plates = _apply_schema_source_columns(
        base,
        source_plates,
        schema,
        group_col="IP_GRID_NUM",
        explicit_inputs=inputs,
        ignore_columns={"IP_ID", "IP_GRID_NUM", "IP_BARCODE", "IP_BARCODE_SUFFIX", "IP_INS_DATE", source_id_col},
    )
    plates[grid_col] = base["IP_GRID_NUM"]
    plates[source_id_col] = source_plates.groupby(grid_col)[id_col].first().reindex(plates[grid_col]).to_numpy()

    # Barcode special logic. A mapping keyed by grid is accepted explicitly.
    barcode_input = inputs.get(barcode_col)
    if barcode_input is None:
        plates[barcode_col] = plates[grid_col]
    elif isinstance(barcode_input, Mapping):
        keyed = {str(k): v for k, v in barcode_input.items()}
        plates[barcode_col] = plates[grid_col].astype(str).map(keyed)
        if plates[barcode_col].isna().any():
            missing = plates.loc[plates[barcode_col].isna(), grid_col].drop_duplicates().tolist()
            raise ValueError(f"IP_BARCODE mapping is missing entries for: {missing}")
    else:
        resolved = _resolve_value_from_source(source_plates, barcode_input, group_key=plates[grid_col], label="IP_BARCODE")
        if isinstance(resolved, pd.Series):
            src_col = _find_matching_source_column(source_plates, barcode_input) if isinstance(barcode_input, str) else None
            if src_col is not None:
                grouped = source_plates[[grid_col, src_col]].dropna(subset=[src_col]).groupby(grid_col)[src_col].first()
                plates[barcode_col] = plates[grid_col].map(grouped)
            else:
                plates[barcode_col] = resolved.to_numpy()
        else:
            plates[barcode_col] = resolved

    suffix = inputs.get("IP_BARCODE_SUFFIX")
    if suffix is not None and suffix != "":
        resolved_suffix = _resolve_value_from_source(source_plates, suffix, group_key=plates[grid_col], label="IP_BARCODE_SUFFIX")
        if isinstance(resolved_suffix, pd.Series):
            src_col = _find_matching_source_column(source_plates, suffix) if isinstance(suffix, str) else None
            if src_col is not None:
                grouped = source_plates[[grid_col, src_col]].dropna(subset=[src_col]).groupby(grid_col)[src_col].first()
                plates[barcode_col] = plates[barcode_col].astype(str) + plates[grid_col].map(grouped).astype(str)
            else:
                plates[barcode_col] = plates[barcode_col].astype(str) + resolved_suffix.astype(str).to_numpy()
        else:
            plates[barcode_col] = plates[barcode_col].astype(str) + str(resolved_suffix)

    # Refresh insert date unless explicitly supplied through an input/override mapping.
    explicit_date = inputs.get("IP_INS_DATE")
    if explicit_date is None:
        plates["IP_INS_DATE"] = now_string()
    else:
        resolved_date = _resolve_value_from_source(source_plates, explicit_date, group_key=plates[grid_col], label="IP_INS_DATE")
        if isinstance(resolved_date, pd.Series):
            src_col = _find_matching_source_column(source_plates, explicit_date) if isinstance(explicit_date, str) else None
            if src_col is not None:
                grouped = source_plates[[grid_col, src_col]].dropna(subset=[src_col]).groupby(grid_col)[src_col].first()
                plates["IP_INS_DATE"] = plates[grid_col].map(grouped)
            else:
                plates["IP_INS_DATE"] = resolved_date.to_numpy()
        else:
            plates["IP_INS_DATE"] = resolved_date

    if ins_by is not None and "IP_INS_BY" not in inputs:
        plates["IP_INS_BY"] = ins_by

    # Explicit overrides are always final and can be literals, source columns, or grid-keyed mappings.
    for column, value in (overrides or {}).items():
        if isinstance(value, Mapping):
            keyed = {str(k): v for k, v in value.items()}
            plates[column] = plates[grid_col].astype(str).map(keyed)
            if plates[column].isna().any():
                missing = plates.loc[plates[column].isna(), grid_col].drop_duplicates().tolist()
                raise ValueError(f"Override for {column!r} is missing entries for: {missing}")
        elif isinstance(value, str) and _find_matching_source_column(source_plates, value) is not None:
            src_col = _find_matching_source_column(source_plates, value)
            grouped = source_plates[[grid_col, src_col]].dropna(subset=[src_col]).groupby(grid_col)[src_col].first()
            plates[column] = plates[grid_col].map(grouped)
        else:
            plates[column] = value
    if ins_by is not None and (overrides or {}).get("IP_INS_BY") is None and "IP_INS_BY" not in plates.columns:
        plates["IP_INS_BY"] = ins_by
    return plates.reset_index(drop=True)


def stage_daughter_wells(
    source_wells: pd.DataFrame,
    *,
    new_ip_id_map: Mapping[Any, Any],
    well_inputs: Mapping[str, Any] | None = None,
    well_overrides: Mapping[str, Any] | None = None,
    ins_by: int | None = None,
    ins_date: str | None = None,
    id_col: str = "IW_ID",
    source_id_col: str = "IW_SOURCE_ID",
    ip_id_col: str = "IW_IP_ID",
) -> pd.DataFrame:
    """Shared daughter-well core used by copy/interleave/quadrant/custom paths.

    All source ``IW_*`` columns are preserved by default. Explicit ``well_inputs`` and
    ``well_overrides`` win. ``IW_IP_ID`` is the staged daughter ID. ``IW_BARCODE`` is
    the scanned matrix-tube barcode, so it is not copied from the parent: a daughter
    well plate has no tube barcodes. It is empty unless ``well_inputs`` /
    ``well_overrides`` supply one. It is never filled from ``IP_BARCODE``.
    ``IW_ICL_ID`` is never synthesized or replaced, and parent wells without one
    (empty wells) are dropped rather than copied.
    """
    validate_required_columns(source_wells, [id_col, ip_id_col], "source wells")
    # Empty parent wells (no compound) are not copied to the daughter plate.
    source_wells = drop_empty_wells(source_wells)
    wells = source_wells.copy()
    parent_ip_ids = wells[ip_id_col].copy()
    wells[source_id_col] = wells[id_col]
    wells = wells.drop(columns=[id_col])
    wells[ip_id_col] = parent_ip_ids.map(dict(new_ip_id_map))
    if wells[ip_id_col].isna().any():
        missing = parent_ip_ids.loc[wells[ip_id_col].isna()].drop_duplicates().tolist()
        raise ValueError(f"Could not remap {ip_id_col} for parent plate ids: {missing}")

    # Explicit well inputs/overrides are the only things allowed to change normal IW_* data.
    # IW_ICL_ID is deliberately not special-cased: it remains whatever the source well
    # carried (normally the compound/lot ID). Only IW_BARCODE gets daughter-specific logic.
    for mapping_name, mapping_values in (("well_inputs", well_inputs or {}), ("well_overrides", well_overrides or {})):
        for column, value in mapping_values.items():
            if column in {"IW_BARCODE", "IW_INS_DATE", "IW_INS_BY"}:
                continue
            if isinstance(value, Mapping):
                keyed = {str(k): v for k, v in value.items()}
                wells[column] = wells[ip_id_col].astype(str).map(keyed)
                if wells[column].isna().any():
                    missing = wells.loc[wells[column].isna(), ip_id_col].drop_duplicates().tolist()
                    raise ValueError(f"{mapping_name} override for {column!r} is missing entries for: {missing}")
            elif isinstance(value, str):
                source_col = _find_matching_source_column(source_wells, value)
                if source_col is not None:
                    wells[column] = source_wells[source_col].to_numpy()
                else:
                    if column in NUMERIC_WELL_FIELDS:
                        value = _coerce_numeric_well_input(source_wells, column, value)
                    wells[column] = value
            else:
                wells[column] = value

    explicit_numeric = [
        k for k in NUMERIC_WELL_FIELDS
        if k in (well_inputs or {}) or k in (well_overrides or {})
    ]
    # Only enforce the conc/volume rule when the caller is overriding these fields;
    # a plain copy inherits whatever the parent wells already carry.
    wells = _validate_numeric_well_fields(wells, explicit_numeric, enforce_concentration_rule=bool(explicit_numeric))

    barcode_input = (well_overrides or {}).get("IW_BARCODE", (well_inputs or {}).get("IW_BARCODE"))
    if barcode_input is not None:
        if isinstance(barcode_input, Mapping):
            keyed = {str(k): v for k, v in barcode_input.items()}
            wells["IW_BARCODE"] = wells[ip_id_col].astype(str).map(keyed)
            if wells["IW_BARCODE"].isna().any():
                missing = wells.loc[wells["IW_BARCODE"].isna(), ip_id_col].drop_duplicates().tolist()
                raise ValueError(f"IW_BARCODE mapping is missing entries for: {missing}")
        elif isinstance(barcode_input, str) and _find_matching_source_column(source_wells, barcode_input) is not None:
            wells["IW_BARCODE"] = source_wells[_find_matching_source_column(source_wells, barcode_input)].to_numpy()
        else:
            wells["IW_BARCODE"] = barcode_input
    else:
        wells["IW_BARCODE"] = None  # tube barcodes don't transfer to a daughter well plate

    # Insert metadata is refreshed for every daughter registration unless explicitly
    # supplied through well_inputs/well_overrides.
    ins_date_value = (well_overrides or {}).get("IW_INS_DATE", (well_inputs or {}).get("IW_INS_DATE"))
    if ins_date_value is None:
        wells["IW_INS_DATE"] = ins_date if ins_date is not None else now_string()
    else:
        if isinstance(ins_date_value, str):
            source_col = _find_matching_source_column(source_wells, ins_date_value)
            if source_col is not None:
                wells["IW_INS_DATE"] = source_wells[source_col].to_numpy()
            else:
                wells["IW_INS_DATE"] = ins_date_value
        else:
            wells["IW_INS_DATE"] = ins_date_value
    ins_by_value = (well_overrides or {}).get("IW_INS_BY", (well_inputs or {}).get("IW_INS_BY"))
    if ins_by_value is not None:
        if isinstance(ins_by_value, str):
            source_col = _find_matching_source_column(source_wells, ins_by_value)
            if source_col is not None:
                wells["IW_INS_BY"] = source_wells[source_col].to_numpy()
            else:
                wells["IW_INS_BY"] = ins_by_value
        else:
            wells["IW_INS_BY"] = ins_by_value
    elif ins_by is not None:
        wells["IW_INS_BY"] = ins_by
    return wells.reset_index(drop=True)


def copy_daughter_wells(
    source_wells: pd.DataFrame,
    *,
    new_ip_id_map: Mapping[Any, Any],
    well_inputs: Mapping[str, Any] | None = None,
    well_overrides: Mapping[str, Any] | None = None,
    vol_ul: Any | None = None,
    conc_um: Any | None = None,
    final_conc_nm: Any | None = None,
    ins_by: int | None = None,
    ins_date: str | None = None,
    id_col: str = "IW_ID",
    source_id_col: str = "IW_SOURCE_ID",
    ip_id_col: str = "IW_IP_ID",
) -> pd.DataFrame:
    """Build exact-copy daughter wells through the shared staged-well core.

    ``IW_ICL_ID`` is carried from the parent well unchanged. ``IW_BARCODE`` is left
    empty unless supplied (see stage_daughter_wells).
    """
    combined = dict(well_inputs or {})
    if vol_ul is not None:
        combined["IW_VOL_UL"] = vol_ul
    if conc_um is not None:
        combined["IW_CONC_UM"] = conc_um
    if final_conc_nm is not None:
        combined["IW_FINAL_CONC_NM"] = final_conc_nm
    return stage_daughter_wells(
        source_wells,
        new_ip_id_map=new_ip_id_map,
        well_inputs=combined,
        well_overrides=well_overrides,
        ins_by=ins_by,
        ins_date=ins_date,
        id_col=id_col,
        source_id_col=source_id_col,
        ip_id_col=ip_id_col,
    )
def append_frame(engine: Engine, table_name: str, frame: pd.DataFrame) -> None:
    """Append a dataframe to a database table using pandas.to_sql.

    Parameters:
        engine: SQLAlchemy engine connected to destination database.
        table_name: Destination table name.
        frame: Dataframe to append.
    """
    frame.to_sql(name=table_name, con=engine, if_exists="append", index=False)