"""Clean up duplicate INV_COMP_LOT rows from partial registration runs.

DRY RUN is the default. Use --execute to commit.

Selection logic
---------------
1. Load all *_reg.xlsx files under SOURCE_ROOT.
2. LOT_ID in the registration files is INV_COMP_LOT.ICL_ID.
3. Compute keep_num_lots from true duplicates in the successful registration set:
   - a row qualifies for the true-duplicate pool when its platewell is unique
     in source but its smdc occurs more than once.
   - keep_num_lots is the number of qualifying rows for the same (smdc, alias).
   - everything else gets keep_num_lots = 1.
4. For every source LOT_ID, look up its ICL_IC_ID and ICL_ALIAS in the DB.
5. For each (ICL_IC_ID, ICL_ALIAS), keep the newest N collection lots by
   ICL_ID descending, where N = keep_num_lots.
6. Delete all older collection lots for those compounds.
7. For retained collection lots, reset ICL_LOT_NUM to consecutive values starting
   at MAX(existing non-collection ICL_LOT_NUM) + 1, or 1 if none exist.
8. Apply red_flags.csv actions before deleting INV_COMP_LOT:
   - "do not delete" => abort if linked rows exist.
   - "delete from here first if exists" => delete those direct child rows,
     but abort if those child rows themselves have downstream FK data.
   - "ignore - separate on purpose" => report only.

No write occurs unless --execute is supplied. All writes occur in one transaction.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
import numpy as np
from sqlalchemy import MetaData, bindparam, text
from sqlalchemy import Table, delete
from sqlalchemy.engine import Engine



DEFAULT_RED_FLAGS = Path(__file__).with_name("red_flags.csv")
ACTION_DELETE = "delete from here first if exists"
ACTION_BLOCK = "do not delete"
ACTION_IGNORE = "ignore - separate on purpose"


def chunked(items, size: int):
    """Yield successive chunks of `size` items."""
    for start in range(0, len(items), size):
        yield items[start:start + size]


def validate_source(source: pd.DataFrame) -> pd.DataFrame:

    required = {
        "lot_id",
        "plate",
        "row",
        "col",
        "well",
        "comments",
        "smdc",
        "alias",
        "platewell"
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise RuntimeError(f"Source data is missing required columns: {missing}")

    source["lot_id"] = pd.to_numeric(source["lot_id"], errors="raise").astype("int64")
    source["comments"] = source["comments"].fillna("").astype(str)
    source["smdc"] = source["smdc"].astype(str)
    source["alias"] = source["alias"].astype(str)

    return source.sort_values(["platewell", "lot_id"]).reset_index(drop=True)


def add_keep_num_lots(source: pd.DataFrame) -> pd.DataFrame:
    """Add keep_num_lots using the explicit definition."""

    source = source.copy()

    # Compute the number of unique platewells and unique LOT_IDs for each (smdc, alias) group.
    group_stats = (
        source.groupby(["smdc", "alias"])
        .agg(
            n_rows=("lot_id", "nunique"),
            n_platewells=("platewell", "nunique"),
        )
    )

    # Determine how many lots to keep based on the number of unique platewells.
    group_stats["keep_num_lots"] = np.where(
        group_stats["n_platewells"] > 1,
        group_stats["n_platewells"],
        1,
    )

    # Map the keep_num_lots back to the source DataFrame based on smdc and alias.
    source["keep_num_lots"] = (
        source.set_index(["smdc", "alias"]).index
        .map(group_stats["keep_num_lots"])
        .astype(int)
    )

    return source


def load_red_flags(path: Path) -> pd.DataFrame:
    flags = pd.read_csv(path)
    required = {
        "object_a",
        "parent_columns",
        "object_b",
        "child_columns",
        "action",
    }
    missing = sorted(required - set(flags.columns))
    if missing:
        raise RuntimeError(f"red_flags.csv is missing required columns: {missing}")

    flags = flags.copy()
    flags["action_norm"] = flags["action"].fillna("").astype(str).str.strip().str.lower()
    return flags


def _resolve_source_identities_by_smdc(
    conn,
    source: pd.DataFrame,
    batch_size: int = 5000,
) -> pd.DataFrame:
    """Resolve source SMDC values to INV_COMP_LOT identity IDs.

    The registration LOT_ID is a useful exact key while the row exists, but it
    is not a durable key after a prior cleanup.  SMDC is the source-level
    compound identifier, so use INV_COMPOUND.IC_COMPOUND_ID as the fallback.
    """
    source_identity = (
        source[["smdc", "alias", "keep_num_lots"]]
        .drop_duplicates()
        .copy()
    )
    source_identity["smdc_key"] = source_identity["smdc"].astype(str).str.strip()

    smdcs = sorted(source_identity["smdc_key"].dropna().unique().tolist())
    found: list[dict] = []

    for batch in chunked(smdcs, batch_size):
        stmt = text(
            """
            SELECT
                IC_ID,
                IC_COMPOUND_ID
            FROM INV_COMPOUND
            WHERE CAST(IC_COMPOUND_ID AS CHAR) IN :smdcs
            """
        ).bindparams(bindparam("smdcs", expanding=True))
        found.extend(conn.execute(stmt, {"smdcs": batch}).mappings().all())

    mapping = pd.DataFrame(found)
    if mapping.empty:
        return pd.DataFrame(
            columns=["smdc", "alias", "keep_num_lots", "ICL_IC_ID"]
        )

    mapping["smdc_key"] = mapping["IC_COMPOUND_ID"].astype(str).str.strip()
    ambiguous = (
        mapping.groupby("smdc_key")["IC_ID"]
        .nunique()
        .loc[lambda x: x > 1]
    )
    if not ambiguous.empty:
        raise RuntimeError(
            "ABORT: source SMDC maps to multiple INV_COMPOUND.IC_ID values:\n"
            + ambiguous.to_string()
        )

    resolved = source_identity.merge(
        mapping[["smdc_key", "IC_ID"]].drop_duplicates(),
        on="smdc_key",
        how="left",
    )
    resolved = resolved.rename(columns={"IC_ID": "ICL_IC_ID"})
    return resolved[["smdc", "alias", "keep_num_lots", "ICL_IC_ID"]]


def check_already_clean(
    conn,
    source: pd.DataFrame,
    collection: str,
    keep_cols: list,
) -> pd.DataFrame | None:
    """Return a reconstructed source when the database already satisfies retention.

    LOT_ID is the authoritative identity for this workflow. A source LOT_ID that
    is no longer present means that lot was already removed/reconciled; report it
    as a warning and let the normal source_to_db_targets path continue with the
    remaining LOT_IDs. Do not attempt to reconstruct missing LOT_IDs via SMDC.
    """
    lot_ids = sorted(source["lot_id"].drop_duplicates().astype(int).tolist())
    if not lot_ids:
        raise RuntimeError("ABORT: source contains no LOT_IDs.")

    found: list[dict] = []
    for batch in chunked(lot_ids, 5000):
        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_ID IN :ids
            """
        ).bindparams(bindparam("ids", expanding=True))
        found.extend(conn.execute(stmt, {"ids": batch}).mappings().all())

    found_df = pd.DataFrame(found)
    found_ids = set(found_df["ICL_ID"].astype(int)) if not found_df.empty else set()
    missing = sorted(set(lot_ids) - found_ids)

    if missing:
        print(
            f"WARNING: {len(missing):,} source lot IDs are not in the database "
            "(already removed/reconciled); continuing without them."
        )
        return None

    exact = source[
        ["lot_id", "smdc", "alias", "keep_num_lots"]
    ].drop_duplicates(subset=["lot_id"]).merge(
        found_df,
        left_on="lot_id",
        right_on="ICL_ID",
        how="inner",
        validate="one_to_one",
    )

    present_wrong_collection = exact[
        exact["ICL_VEN_COLLECTION"].notna()
        & exact["ICL_VEN_COLLECTION"].ne(collection)
    ]
    if not present_wrong_collection.empty:
        sample = present_wrong_collection[["lot_id", "ICL_VEN_COLLECTION"]].head(20)
        raise RuntimeError(
            f"ABORT: {len(present_wrong_collection):,} source LOT_IDs are in a different collection.\n"
            + sample.to_string(index=False)
        )

    alias_mismatch = exact[
        exact["ICL_ALIAS"].notna()
        & exact["alias"].astype(str).ne(exact["ICL_ALIAS"].astype(str))
    ]
    if not alias_mismatch.empty:
        sample = alias_mismatch[["lot_id", "alias", "ICL_ALIAS"]].head(20)
        raise RuntimeError(
            "ABORT: source alias does not match INV_COMP_LOT.ICL_ALIAS for some LOT_IDs:\n"
            + sample.to_string(index=False)
        )

    # With all LOT_IDs present, use the exact LOT_ID -> IC_ID mapping already in
    # INV_COMP_LOT. This avoids any SMDC fallback and its historical ambiguities.
    expected = (
        exact[["ICL_IC_ID", "ICL_ALIAS", "keep_num_lots"]]
        .drop_duplicates()
        .groupby(["ICL_IC_ID", "ICL_ALIAS"], as_index=False)["keep_num_lots"]
        .max()
    )

    ic_ids = sorted(expected["ICL_IC_ID"].astype(int).unique().tolist())
    actual_rows: list[dict] = []
    for batch in chunked(ic_ids, 5000):
        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_IC_ID IN :ic_ids
              AND ICL_VEN_COLLECTION = :collection
            """
        ).bindparams(bindparam("ic_ids", expanding=True))
        actual_rows.extend(
            conn.execute(
                stmt, {"ic_ids": batch, "collection": collection}
            ).mappings().all()
        )

    actual = pd.DataFrame(actual_rows)
    if actual.empty:
        return None

    actual_counts = (
        actual.groupby(["ICL_IC_ID", "ICL_ALIAS"], as_index=False)
        .size()
        .rename(columns={"size": "actual_count"})
    )
    check = expected.merge(
        actual_counts,
        on=["ICL_IC_ID", "ICL_ALIAS"],
        how="left",
    )
    check["actual_count"] = check["actual_count"].fillna(0).astype(int)

    # Extra aliases for the same source compounds mean cleanup is still needed.
    source_pairs = expected[["ICL_IC_ID", "ICL_ALIAS"]].drop_duplicates()
    unexpected = actual.merge(
        source_pairs,
        on=["ICL_IC_ID", "ICL_ALIAS"],
        how="left",
        indicator=True,
    )
    unexpected = unexpected[unexpected["_merge"] == "left_only"]
    if not unexpected.empty:
        print(
            f"Database has {unexpected['ICL_ID'].nunique():,} collection rows for "
            "the source compounds outside the expected aliases; cleanup is required."
        )
        return None

    bad = check[check["actual_count"] != check["keep_num_lots"]]
    if not bad.empty:
        print(
            "Existing database state does not reconcile with source retention counts; "
            "normal cleanup planning will continue."
        )
        print(bad.head(20).to_string(index=False))
        return None

    print("Database already reconciles with the expected retention counts.")
    print("No delete/update transaction will be started.")
    return get_new_source(conn, collection, source, keep_cols)

def source_to_db_targets(
    conn,
    source: pd.DataFrame,
    collection: str,
    batch_size: int = 1000,
) -> pd.DataFrame:
    """Validate source identities and return DB identity/retention metadata.

    LOT_ID is the authoritative identity. Missing LOT_IDs are treated as lots
    that were already removed/reconciled and are excluded from the DB target set.
    They are never reconstructed through SMDC.
    """
    lot_ids = sorted(source["lot_id"].drop_duplicates().astype(int).tolist())
    found: list[dict] = []

    for chunk in chunked(lot_ids, batch_size):
        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_ID IN :ids
            """
        ).bindparams(bindparam("ids", expanding=True))
        found.extend(conn.execute(stmt, {"ids": chunk}).mappings().all())

    db = pd.DataFrame(found)
    if db.empty:
        db = pd.DataFrame(
            columns=[
                "ICL_ID", "ICL_IC_ID", "ICL_ALIAS",
                "ICL_VEN_COLLECTION", "ICL_LOT_NUM",
            ]
        )

    expected_lots = pd.DataFrame({"lot_id": lot_ids})
    merged = expected_lots.merge(
        db,
        left_on="lot_id",
        right_on="ICL_ID",
        how="left",
        indicator=True,
    )
    missing = merged.loc[merged["_merge"] != "both", "lot_id"].astype(int).tolist()
    if missing:
        print(
            f"WARNING: {len(missing):,} source lot IDs are not in the database "
            "(already removed/reconciled); continuing without them."
        )

    present = merged[merged["_merge"] == "both"].copy()
    wrong_collection = present.loc[
        present["ICL_VEN_COLLECTION"].notna()
        & present["ICL_VEN_COLLECTION"].ne(collection),
        "lot_id",
    ].astype(int).tolist()
    if wrong_collection:
        raise RuntimeError(
            f"ABORT: {len(wrong_collection)} source LOT_IDs are in a different collection "
            f"from {collection!r}. First IDs: {wrong_collection[:20]}"
        )

    source_alias = source[
        ["lot_id", "smdc", "alias", "keep_num_lots"]
    ].drop_duplicates(subset=["lot_id"])
    present_source = source_alias.merge(
        db,
        left_on="lot_id",
        right_on="ICL_ID",
        how="inner",
        validate="one_to_one",
    )
    alias_mismatch = present_source[
        present_source["alias"].astype(str).ne(present_source["ICL_ALIAS"].astype(str))
    ]
    if not alias_mismatch.empty:
        sample = alias_mismatch[["lot_id", "alias", "ICL_ALIAS"]].head(20)
        raise RuntimeError(
            "ABORT: source alias does not match INV_COMP_LOT.ICL_ALIAS for some LOT_IDs:\n"
            + sample.to_string(index=False)
        )

    # Preserve exact DB mappings for LOT_IDs still present.
    present_out = present_source[
        [
            "lot_id", "smdc", "alias", "keep_num_lots",
            "ICL_ID", "ICL_IC_ID", "ICL_ALIAS",
            "ICL_VEN_COLLECTION", "ICL_LOT_NUM",
        ]
    ].copy()

    missing_source = source_alias[~source_alias["lot_id"].isin(present_source["lot_id"])]
    if not missing_source.empty:
        # Missing LOT_IDs were already reported above. They are intentionally
        # excluded from the DB target set because they were already
        # removed/reconciled and must not be reconstructed from SMDC.
        return present_out

    return present_out

def get_retention_plan(
    conn,
    source_db: pd.DataFrame,
    collection: str,
    batch_size: int = 5000,
) -> pd.DataFrame:
    """Return all relevant collection ICL rows with KEEP/DELETE decisions.

    Matching is based on the composite identity:
        (ICL_IC_ID, ICL_ALIAS)

    keep_num_lots is the maximum requested number for each identity.
    Within each identity, the newest rows are determined by descending ICL_ID.
    """

    identities = (
        source_db[
            ["ICL_IC_ID", "ICL_ALIAS", "keep_num_lots"]
        ]
        .drop_duplicates()
        .groupby(
            ["ICL_IC_ID", "ICL_ALIAS"],
            as_index=False,
        )["keep_num_lots"]
        .max()
    )

    if identities.empty:
        raise RuntimeError("ABORT: no source identities found.")

    # ------------------------------------------------------------
    # Fetch all collection rows for all relevant ICL_IC_IDs.
    #
    # We only need ICL_IC_ID in the IN clause because ICL_ALIAS is
    # validated in pandas against the exact (ICL_IC_ID, ICL_ALIAS)
    # identities from source_db.
    # ------------------------------------------------------------

    ic_ids = (
        identities["ICL_IC_ID"]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    matches = []

    for batch in chunked(ic_ids, batch_size):

        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_IC_ID IN :ic_ids
              AND ICL_VEN_COLLECTION = :collection
            """
        ).bindparams(
            bindparam("ic_ids", expanding=True)
        )

        result = conn.execute(
            stmt,
            {
                "ic_ids": batch,
                "collection": collection,
            },
        )

        matches.extend(result.mappings().all())

    if not matches:
        raise RuntimeError(
            f"ABORT: no {collection} rows found for any relevant ICL_IC_ID."
        )

    plan = pd.DataFrame(matches)

    # ------------------------------------------------------------
    # Keep ONLY exact (ICL_IC_ID, ICL_ALIAS) identities from source.
    # ------------------------------------------------------------

    identities = identities.copy()
    identities["ICL_IC_ID"] = identities["ICL_IC_ID"].astype(int)

    plan["ICL_IC_ID"] = plan["ICL_IC_ID"].astype(int)

    plan = plan.merge(
        identities,
        on=["ICL_IC_ID", "ICL_ALIAS"],
        how="inner",
        validate="many_to_one",
    )

    # ------------------------------------------------------------
    # Detect source identities that had no DB rows.
    # ------------------------------------------------------------

    found_identities = (
        plan[
            ["ICL_IC_ID", "ICL_ALIAS"]
        ]
        .drop_duplicates()
        .assign(_found=True)
    )

    missing = (
        identities
        .merge(
            found_identities,
            on=["ICL_IC_ID", "ICL_ALIAS"],
            how="left",
        )
        .loc[lambda x: x["_found"].isna()]
        .drop(columns="_found")
    )

    if not missing.empty:
        raise RuntimeError(
            f"ABORT: no matching {collection} rows found for these "
            f"(ICL_IC_ID, ICL_ALIAS) identities:\n"
            + missing.to_string(index=False)
        )

    # ------------------------------------------------------------
    # Newest rows first within each composite identity.
    # ------------------------------------------------------------

    plan = plan.sort_values(
        ["ICL_IC_ID", "ICL_ALIAS", "ICL_ID"],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    plan["rank_newest"] = (
        plan.groupby(
            ["ICL_IC_ID", "ICL_ALIAS"]
        ).cumcount() + 1
    )

    plan["keep_num_lots"] = plan["keep_num_lots"].astype(int)

    plan["action"] = np.where(
        plan["rank_newest"] <= plan["keep_num_lots"],
        "KEEP",
        "DELETE",
    )

    return plan[
        [
            "ICL_ID",
            "ICL_IC_ID",
            "ICL_ALIAS",
            "ICL_VEN_COLLECTION",
            "ICL_LOT_NUM",
            "keep_num_lots",
            "rank_newest",
            "action",
        ]
    ]


def get_reset_plan(
    conn,
    kept: pd.DataFrame,
    collection: str,
    delete_ids: list[int],
    batch_size: int = 5000,
) -> pd.DataFrame:
    """Build ICL_LOT_NUM reset plan.

    For each ICL_IC_ID:

    - Ignore all rows from `collection` when determining the base lot number.
    - Find MAX(ICL_LOT_NUM) among the other collections that will survive.
    - Retained rows in `collection` are assigned consecutive lot numbers
      beginning at max_other_lot_num + 1.
    - Retained TargetMol rows are ordered by ICL_ID ascending so that their
      lot numbers follow insertion order.
    """

    if kept.empty:
        return pd.DataFrame()

    kept = kept.copy()
    kept["ICL_ID"] = kept["ICL_ID"].astype(int)
    kept["ICL_IC_ID"] = kept["ICL_IC_ID"].astype(int)

    ic_ids = (
        kept["ICL_IC_ID"]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    delete_ids = [int(x) for x in delete_ids]

    # ------------------------------------------------------------
    # Helper for batched IN queries.
    # ------------------------------------------------------------
    def fetch_in_batches(
        sql: str,
        ids: list[int],
        id_param: str,
        extra_params: dict | None = None,
    ) -> pd.DataFrame:

        pieces = []

        for batch in chunked(ids, batch_size):
            stmt = text(sql).bindparams(
                bindparam(id_param, expanding=True)
            )

            params = {id_param: batch}

            if extra_params:
                params.update(extra_params)

            pieces.append(
                pd.read_sql_query(
                    stmt,
                    conn,
                    params=params,
                )
            )

        if not pieces:
            return pd.DataFrame()

        return pd.concat(pieces, ignore_index=True)

    # ------------------------------------------------------------
    # 1. MAX lot number among OTHER collections that will survive.
    #
    # This is deliberately NOT based on TargetMol rows.
    # ------------------------------------------------------------
    max_other = fetch_in_batches(
        """
        SELECT ICL_IC_ID,
               MAX(ICL_LOT_NUM) AS max_other_lot_num
        FROM INV_COMP_LOT
        WHERE ICL_IC_ID IN :ic_ids
          AND (
            ICL_VEN_COLLECTION <> :collection
                OR ICL_VEN_COLLECTION IS NULL
            )
          AND ICL_LOT_NUM IS NOT NULL
        GROUP BY ICL_IC_ID
        """,
        ic_ids,
        "ic_ids",
        {"collection": collection},
    )

    max_other_by_ic = (
        max_other
        .set_index("ICL_IC_ID")["max_other_lot_num"]
        .to_dict()
    )

    # ------------------------------------------------------------
    # 2. Determine retained TargetMol rows.
    #
    # `kept` should already represent exactly the TargetMol rows
    # that survive the deletion.
    # ------------------------------------------------------------
    kept = kept.sort_values(
        ["ICL_IC_ID", "ICL_ID"]
    )

    # ------------------------------------------------------------
    # 3. Assign TargetMol lot numbers.
    # ------------------------------------------------------------
    reset_rows = []

    for ic_id, group in kept.groupby("ICL_IC_ID", sort=False):

        ic_id = int(ic_id)

        max_other_lot = max_other_by_ic.get(ic_id)

        next_num = (
            int(max_other_lot) + 1
            if pd.notna(max_other_lot)
            else 1
        )

        for row in group.itertuples(index=False):

            reset_rows.append(
                {
                    "ICL_ID": int(row.ICL_ID),
                    "ICL_IC_ID": ic_id,
                    "ICL_ALIAS": row.ICL_ALIAS,
                    "ICL_LOT_NUM": row.ICL_LOT_NUM,
                    "new_lot_num": next_num,
                }
            )

            next_num += 1

    reset = pd.DataFrame(reset_rows)

    if reset.empty:
        return reset

    reset = reset.sort_values(
        ["ICL_IC_ID", "ICL_ALIAS", "ICL_ID"]
    ).reset_index(drop=True)

    # ------------------------------------------------------------
    # 4. Sanity check: no collisions with surviving OTHER-collection
    # rows.
    #
    # TargetMol rows are excluded because they are exactly the rows
    # being assigned these new numbers.
    # ------------------------------------------------------------
    surviving_other = fetch_in_batches(
        """
        SELECT
            ICL_ID,
            ICL_IC_ID,
            ICL_LOT_NUM,
            ICL_VEN_COLLECTION
        FROM INV_COMP_LOT
        WHERE
            ICL_IC_ID IN :ic_ids
            AND (
                ICL_VEN_COLLECTION <> :collection
                OR ICL_VEN_COLLECTION IS NULL
            )
            AND ICL_LOT_NUM IS NOT NULL
        """,
        ic_ids,
        "ic_ids",
        {"collection": collection},
    )

    collisions = reset.merge(
        surviving_other[
            ["ICL_IC_ID", "ICL_LOT_NUM"]
        ],
        left_on=["ICL_IC_ID", "new_lot_num"],
        right_on=["ICL_IC_ID", "ICL_LOT_NUM"],
        how="inner",
    )

    if not collisions.empty:
        raise RuntimeError(
            "ABORT: proposed ICL_LOT_NUM values collide with "
            "existing non-TargetMol rows:\n"
            + collisions[
                ["ICL_IC_ID", "new_lot_num"]
            ].drop_duplicates().to_string(index=False)
        )

    # ------------------------------------------------------------
    # 5. Don't update rows whose lot number is already correct.
    #
    # This makes the operation idempotent and also avoids relying on
    # MySQL rowcount semantics for no-op UPDATEs.
    # ------------------------------------------------------------
    reset = reset[
        reset["ICL_LOT_NUM"].fillna(-1).astype(int)
        != reset["new_lot_num"].astype(int)
    ].copy()

    return reset


def get_direct_links(flags: pd.DataFrame) -> pd.DataFrame:
    """Configured direct/extra actions only; recursive descendants are not acted on."""
    direct = flags[
        flags["object_a"].eq("INV_COMP_LOT")
        & flags["object_b"].notna()
        & flags["child_columns"].notna()
        & flags["action_norm"].isin([ACTION_BLOCK, ACTION_DELETE, ACTION_IGNORE])
    ].copy()

    extra = flags[
        flags["object_a"].eq("INV_ALIAS")
        & flags["child_columns"].notna()
        & flags["action_norm"].isin([ACTION_BLOCK, ACTION_DELETE, ACTION_IGNORE])
    ].copy()

    return pd.concat([direct, extra], ignore_index=True).drop_duplicates(
        subset=["object_a", "object_b", "child_columns", "action_norm"]
    )


def find_primary_key(conn, table: str) -> str:
    rows = conn.execute(
        text(
            """
            SELECT COLUMN_NAME
            FROM information_schema.KEY_COLUMN_USAGE
            WHERE TABLE_SCHEMA = DATABASE()
              AND TABLE_NAME = :table_name
              AND CONSTRAINT_NAME = 'PRIMARY'
            ORDER BY ORDINAL_POSITION
            """
        ),
        {"table_name": table},
    ).scalars().all()
    if len(rows) != 1:
        raise RuntimeError(f"Cannot safely delete from {table}: expected one-column PK, found {rows}")
    return str(rows[0])


def count_linked_rows(conn, table: str, columns: list[str], ids: list[int]) -> int:
    if not ids:
        return 0
    predicates = " OR ".join(f"`{c}` IN :ids" for c in columns)
    stmt = text(f"SELECT COUNT(*) FROM `{table}` WHERE {predicates}").bindparams(
        bindparam("ids", expanding=True)
    )
    return int(conn.execute(stmt, {"ids": ids}).scalar_one())


def linked_primary_keys(conn, table: str, columns: list[str], ids: list[int]) -> list[int]:
    pk = find_primary_key(conn, table)
    predicates = " OR ".join(f"`{c}` IN :ids" for c in columns)
    stmt = text(f"SELECT `{pk}` FROM `{table}` WHERE {predicates}").bindparams(
        bindparam("ids", expanding=True)
    )
    return sorted({int(x) for x in conn.execute(stmt, {"ids": ids}).scalars().all()})


def downstream_children(conn, table: str, pk_column: str, child_ids: list[int]) -> list[dict]:
    if not child_ids:
        return []

    fks = conn.execute(
        text(
            """
            SELECT
                kcu.TABLE_NAME AS child_table,
                kcu.COLUMN_NAME AS child_column,
                rc.DELETE_RULE AS delete_rule,
                kcu.CONSTRAINT_NAME AS constraint_name
            FROM information_schema.KEY_COLUMN_USAGE kcu
            JOIN information_schema.REFERENTIAL_CONSTRAINTS rc
              ON rc.CONSTRAINT_SCHEMA = kcu.CONSTRAINT_SCHEMA
             AND rc.CONSTRAINT_NAME = kcu.CONSTRAINT_NAME
            WHERE kcu.TABLE_SCHEMA = DATABASE()
              AND kcu.REFERENCED_TABLE_NAME = :parent_table
              AND kcu.REFERENCED_COLUMN_NAME = :parent_column
            """
        ),
        {"parent_table": table, "parent_column": pk_column},
    ).mappings().all()

    blockers: list[dict] = []
    for fk in fks:
        # No recursive cleanup here: if the immediate linked table has children,
        # abort rather than attempting to walk or delete the deeper graph.
        n = count_linked_rows(conn, str(fk["child_table"]), [str(fk["child_column"])], child_ids)
        if n:
            blockers.append({**fk, "rows": n})
    return blockers


def preflight_dependencies(conn, flags: pd.DataFrame, target_ids: list[int]) -> None:
    links = get_direct_links(flags)
    blockers: list[dict] = []

    for (object_a, object_b), group in links.groupby(["object_a", "object_b"], dropna=False):
        table = str(object_b if pd.notna(object_b) else object_a)
        actions = set(group["action_norm"].astype(str))
        if len(actions) != 1:
            blockers.append({"table": table, "reason": f"conflicting actions: {sorted(actions)}"})
            continue
        action = next(iter(actions))

        columns: list[str] = []
        for value in group["child_columns"].dropna():
            for col in str(value).split(","):
                col = col.strip()
                if col and col not in columns:
                    columns.append(col)
        if not columns:
            blockers.append({"table": table, "reason": "no child_columns configured"})
            continue

        n = count_linked_rows(conn, table, columns, target_ids)
        print(f"Dependency check: {table} [{', '.join(columns)}] -> {n:,} matching rows; action={action}")

        if action == ACTION_BLOCK and n:
            blockers.append({
                "table": table,
                "columns": ", ".join(columns),
                "rows": n,
                "reason": "red_flags.csv says DO NOT DELETE",
            })
            continue

        if action == ACTION_IGNORE:
            continue

        if action == ACTION_DELETE and n:
            child_pk = find_primary_key(conn, table)
            child_ids = linked_primary_keys(conn, table, columns, target_ids)
            downstream = downstream_children(conn, table, child_pk, child_ids)
            if downstream:
                for row in downstream:
                    blockers.append({
                        "table": table,
                        "reason": (
                            f"immediate child {row['child_table']}.{row['child_column']} "
                            f"has {row['rows']:,} rows; recursive cleanup is intentionally blocked"
                        ),
                    })

    if blockers:
        print("\n=== DEPENDENCY BLOCKERS ===")
        for b in blockers:
            print(b)
        raise RuntimeError(f"ABORT: {len(blockers)} dependency/configuration blocker(s). No rows were changed.")


def delete_direct_dependencies(conn, flags: pd.DataFrame, target_ids: list[int]) -> dict[str, int]:
    deleted: dict[str, int] = {}
    links = get_direct_links(flags)

    for (object_a, object_b), group in links.groupby(["object_a", "object_b"], dropna=False):
        table = str(object_b if pd.notna(object_b) else object_a)
        actions = set(group["action_norm"].astype(str))
        if actions != {ACTION_DELETE}:
            continue

        columns: list[str] = []
        for value in group["child_columns"].dropna():
            for col in str(value).split(","):
                col = col.strip()
                if col and col not in columns:
                    columns.append(col)
        if not columns:
            raise RuntimeError(f"Cannot delete from {table}: no child_columns configured")

        predicates = " OR ".join(f"`{c}` IN :ids" for c in columns)
        stmt = text(f"DELETE FROM `{table}` WHERE {predicates}").bindparams(
            bindparam("ids", expanding=True)
        )
        result = conn.execute(stmt, {"ids": target_ids})
        deleted[table] = int(result.rowcount or 0)

    return deleted


def update_lot_numbers(
    conn,
    reset_plan: pd.DataFrame,
    batch_size: int = 1000,
) -> int:
    """Safely renumber lots in two set-based UPDATE phases."""

    if reset_plan.empty:
        return 0

    rows = [
        (int(row.ICL_ID), int(row.new_lot_num))
        for row in reset_plan.itertuples(index=False)
    ]

    total_updated = 0

    for start in range(0, len(rows), batch_size):
        batch = rows[start:start + batch_size]
        ids = [icl_id for icl_id, _ in batch]

        # ------------------------------------------------------------
        # Phase 1: move all affected rows to unique temporary values.
        # Negative ICL_ID guarantees uniqueness if ICL_LOT_NUM is signed.
        # ------------------------------------------------------------
        temp_stmt = text(
            """
            UPDATE INV_COMP_LOT
            SET ICL_LOT_NUM = -ICL_ID
            WHERE ICL_ID IN :ids
            """
        ).bindparams(
            bindparam("ids", expanding=True)
        )

        conn.execute(temp_stmt, {"ids": ids})

        # ------------------------------------------------------------
        # Phase 2: assign final values with one set-based UPDATE.
        # ------------------------------------------------------------
        case_parts = []
        params = {}

        for i, (icl_id, new_num) in enumerate(batch):
            id_param = f"id_{i}"
            num_param = f"num_{i}"

            case_parts.append(
                f"WHEN :{id_param} THEN :{num_param}"
            )
            params[id_param] = icl_id
            params[num_param] = new_num

        case_sql = "\n".join(case_parts)

        final_stmt = text(
            f"""
            UPDATE INV_COMP_LOT
            SET ICL_LOT_NUM = CASE ICL_ID
                {case_sql}
            END
            WHERE ICL_ID IN :target_ids
            """
        ).bindparams(
            bindparam("target_ids", expanding=True)
        )

        params["target_ids"] = ids

        result = conn.execute(final_stmt, params)
        total_updated += int(result.rowcount or 0)

        completed = min(start + len(batch), len(rows))
        print(
            f"Updated lot numbers: "
            f"{completed:,}/{len(rows):,}"
        )

    return total_updated


def get_reconciliation_report(
    conn,
    source_db: pd.DataFrame,
    collection: str,
    batch_size: int = 5000,
) -> pd.DataFrame:
    """
    Compare expected retention counts from the registration source
    against the actual surviving INV_COMP_LOT rows in the database.

    Returns one row per (ICL_IC_ID, ICL_ALIAS) with:
        expected_keep_num_lots
        actual_lot_count
        count_difference
        status

    Also includes the actual surviving ICL_IDs and lot numbers.
    """

    expected = (
        source_db[
            ["ICL_IC_ID", "ICL_ALIAS", "keep_num_lots"]
        ]
        .drop_duplicates()
        .groupby(
            ["ICL_IC_ID", "ICL_ALIAS"],
            as_index=False,
        )["keep_num_lots"]
        .max()
        .rename(columns={"keep_num_lots": "expected_keep_num_lots"})
    )

    if expected.empty:
        raise RuntimeError("No expected retention identities found.")

    ic_ids = (
        expected["ICL_IC_ID"]
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    rows: list[dict] = []

    for batch in chunked(ic_ids, batch_size):
        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_IC_ID IN :ic_ids
              AND ICL_VEN_COLLECTION = :collection
            ORDER BY ICL_IC_ID, ICL_ALIAS, ICL_ID
            """
        ).bindparams(
            bindparam("ic_ids", expanding=True)
        )

        rows.extend(
            conn.execute(
                stmt,
                {
                    "ic_ids": batch,
                    "collection": collection,
                },
            ).mappings().all()
        )

    actual = pd.DataFrame(rows)

    if actual.empty:
        actual = pd.DataFrame(
            columns=[
                "ICL_ID",
                "ICL_IC_ID",
                "ICL_ALIAS",
                "ICL_VEN_COLLECTION",
                "ICL_LOT_NUM",
            ]
        )

    if not actual.empty:
        actual["ICL_IC_ID"] = actual["ICL_IC_ID"].astype(int)

    # Exact composite identity match.
    actual_counts = (
        actual.groupby(
            ["ICL_IC_ID", "ICL_ALIAS"],
            as_index=False,
        )
        .agg(
            actual_lot_count=("ICL_ID", "nunique"),
            actual_ICL_IDs=("ICL_ID", lambda x: sorted(x.astype(int).tolist())),
            actual_lot_nums=(
                "ICL_LOT_NUM",
                lambda x: sorted(
                    int(v)
                    for v in x.dropna()
                ),
            ),
        )
    )

    report = expected.merge(
        actual_counts,
        on=["ICL_IC_ID", "ICL_ALIAS"],
        how="left",
    )

    report["actual_lot_count"] = (
        report["actual_lot_count"]
        .fillna(0)
        .astype(int)
    )

    report["count_difference"] = (
        report["actual_lot_count"]
        - report["expected_keep_num_lots"]
    )

    report["status"] = np.select(
        [
            report["actual_lot_count"] == report["expected_keep_num_lots"],
            report["actual_lot_count"] < report["expected_keep_num_lots"],
            report["actual_lot_count"] > report["expected_keep_num_lots"],
        ],
        [
            "OK",
            "MISSING_LOTS",
            "EXTRA_LOTS",
        ],
        default="UNKNOWN",
    )

    return report.sort_values(
        ["status", "ICL_IC_ID", "ICL_ALIAS"]
    ).reset_index(drop=True)



def validate_platewell_duplicate_source(source: pd.DataFrame) -> pd.DataFrame:
    """Validate the minimal source schema needed for post-cleanup platewell reconciliation."""
    required = {"lot_id", "platewell"}
    missing = sorted(required - set(source.columns))
    if missing:
        raise RuntimeError(
            "Source data is missing required columns for platewell reconciliation: "
            f"{missing}"
        )

    source = source[["lot_id", "platewell"]].copy()
    source["lot_id"] = pd.to_numeric(source["lot_id"], errors="raise").astype("int64")
    source["platewell"] = source["platewell"].fillna("").astype(str)

    if source["platewell"].eq("").any():
        raise RuntimeError(
            "ABORT: platewell reconciliation requires a non-empty platewell for every source row."
        )

    # A LOT_ID must map to exactly one platewell in the source used for this stage.
    lot_platewell_counts = source.groupby("lot_id")["platewell"].nunique()
    inconsistent = lot_platewell_counts[lot_platewell_counts > 1]
    if not inconsistent.empty:
        raise RuntimeError(
            "ABORT: source LOT_IDs map to multiple platewells; cannot safely reconcile them:\n"
            + inconsistent.to_string()
        )

    return source.drop_duplicates().reset_index(drop=True)


def get_platewell_duplicate_plan(
    conn,
    source: pd.DataFrame,
    collection: str,
    batch_size: int = 5000,
) -> pd.DataFrame:
    """Plan only the *remaining* platewell duplicates after normal lot cleanup.

    A duplicate platewell is actionable only when all surviving lots in that
    platewell map through INV_COMPOUND to one IC_COMPOUND_ID and exactly one
    surviving lot has ICL_IC_ID == IC_COMPOUND_ID.

    Such a group gets one KEEP row and DELETE rows for the other lots. Groups
    that cannot be reconciled this way are returned as MANUAL rows and are
    never deleted automatically.
    """
    source = validate_platewell_duplicate_source(source)

    lot_ids = source["lot_id"].drop_duplicates().astype(int).tolist()
    found: list[dict] = []

    for batch in chunked(lot_ids, batch_size):
        stmt = text(
            """
            SELECT
                ICL_ID,
                ICL_IC_ID,
                ICL_ALIAS,
                ICL_VEN_COLLECTION,
                ICL_LOT_NUM
            FROM INV_COMP_LOT
            WHERE ICL_ID IN :ids
            """
        ).bindparams(bindparam("ids", expanding=True))
        found.extend(conn.execute(stmt, {"ids": batch}).mappings().all())

    if not found:
        return pd.DataFrame(
            columns=[
                "platewell",
                "ICL_ID",
                "ICL_IC_ID",
                "ICL_ALIAS",
                "ICL_LOT_NUM",
                "IC_COMPOUND_ID",
                "action",
                "reason",
                "keep_ICL_ID",
            ]
        )

    db = pd.DataFrame(found)

    wrong_collection = db.loc[
        db["ICL_VEN_COLLECTION"] != collection, "ICL_ID"
    ].astype(int).tolist()
    if wrong_collection:
        raise RuntimeError(
            f"ABORT: source LOT_IDs present in INV_COMP_LOT but not in {collection!r}. "
            f"First IDs: {wrong_collection[:20]}"
        )

    current = source.merge(
        db,
        left_on="lot_id",
        right_on="ICL_ID",
        how="inner",
        validate="one_to_one",
    )

    # Lots already deleted by the retention stage are intentionally absent here.
    duplicate_wells = (
        current.groupby("platewell")["ICL_ID"]
        .nunique()
        .loc[lambda x: x > 1]
        .index
        .tolist()
    )

    if not duplicate_wells:
        return pd.DataFrame(
            columns=[
                "platewell",
                "ICL_ID",
                "ICL_IC_ID",
                "ICL_ALIAS",
                "ICL_LOT_NUM",
                "IC_COMPOUND_ID",
                "action",
                "reason",
                "keep_ICL_ID",
            ]
        )

    dupes = current[current["platewell"].isin(duplicate_wells)].copy()
    ic_ids = sorted(dupes["ICL_IC_ID"].dropna().astype(int).unique().tolist())

    compound_rows: list[dict] = []
    for batch in chunked(ic_ids, batch_size):
        stmt = text(
            """
            SELECT
                IC_ID,
                IC_COMPOUND_ID
            FROM INV_COMPOUND
            WHERE IC_ID IN :ic_ids
            """
        ).bindparams(bindparam("ic_ids", expanding=True))
        compound_rows.extend(conn.execute(stmt, {"ic_ids": batch}).mappings().all())

    compound = pd.DataFrame(compound_rows, columns=["IC_ID", "IC_COMPOUND_ID"])
    if not compound.empty:
        compound["IC_ID"] = pd.to_numeric(compound["IC_ID"], errors="raise").astype(int)

    dupes = dupes.merge(
        compound,
        left_on="ICL_IC_ID",
        right_on="IC_ID",
        how="left",
        validate="many_to_one",
    ).drop(columns=["IC_ID"])

    plan_rows: list[dict] = []

    for platewell, group in dupes.groupby("platewell", sort=True):
        group = group.sort_values("ICL_ID").copy()
        mapped = group["IC_COMPOUND_ID"].dropna().unique().tolist()
        missing_mapping = group["IC_COMPOUND_ID"].isna().any()
        distinct_ic_ids = group["ICL_IC_ID"].nunique()

        if missing_mapping:
            reason = "One or more duplicate lots have no INV_COMPOUND.IC_COMPOUND_ID mapping."
            for row in group.itertuples(index=False):
                plan_rows.append({
                    "platewell": platewell,
                    "ICL_ID": int(row.ICL_ID),
                    "ICL_IC_ID": int(row.ICL_IC_ID),
                    "ICL_ALIAS": row.ICL_ALIAS,
                    "ICL_LOT_NUM": row.ICL_LOT_NUM,
                    "IC_COMPOUND_ID": row.IC_COMPOUND_ID,
                    "action": "MANUAL",
                    "reason": reason,
                    "keep_ICL_ID": pd.NA,
                })
            continue

        if len(mapped) != 1:
            reason = (
                "Duplicate platewell lots do not reconcile to one IC_COMPOUND_ID: "
                f"{mapped}"
            )
            for row in group.itertuples(index=False):
                plan_rows.append({
                    "platewell": platewell,
                    "ICL_ID": int(row.ICL_ID),
                    "ICL_IC_ID": int(row.ICL_IC_ID),
                    "ICL_ALIAS": row.ICL_ALIAS,
                    "ICL_LOT_NUM": row.ICL_LOT_NUM,
                    "IC_COMPOUND_ID": row.IC_COMPOUND_ID,
                    "action": "MANUAL",
                    "reason": reason,
                    "keep_ICL_ID": pd.NA,
                })
            continue

        ic_compound_id = int(mapped[0])
        canonical = group[group["ICL_IC_ID"].astype(int) == ic_compound_id]

        if distinct_ic_ids <= 1 or len(canonical) != 1:
            reason = (
                "Platewell duplicate does not have exactly one lot whose ICL_IC_ID "
                "equals the reconciled IC_COMPOUND_ID; manual reconciliation required."
            )
            for row in group.itertuples(index=False):
                plan_rows.append({
                    "platewell": platewell,
                    "ICL_ID": int(row.ICL_ID),
                    "ICL_IC_ID": int(row.ICL_IC_ID),
                    "ICL_ALIAS": row.ICL_ALIAS,
                    "ICL_LOT_NUM": row.ICL_LOT_NUM,
                    "IC_COMPOUND_ID": ic_compound_id,
                    "action": "MANUAL",
                    "reason": reason,
                    "keep_ICL_ID": pd.NA,
                })
            continue

        keep_icl_id = int(canonical.iloc[0]["ICL_ID"])
        for row in group.itertuples(index=False):
            row_icl_id = int(row.ICL_ID)
            action = "KEEP" if row_icl_id == keep_icl_id else "DELETE"
            reason = (
                "ICL_IC_ID matches IC_COMPOUND_ID; this is the canonical lot."
                if action == "KEEP"
                else "Same platewell and reconciled IC_COMPOUND_ID; remove non-canonical IC_ID lot."
            )
            plan_rows.append({
                "platewell": platewell,
                "ICL_ID": row_icl_id,
                "ICL_IC_ID": int(row.ICL_IC_ID),
                "ICL_ALIAS": row.ICL_ALIAS,
                "ICL_LOT_NUM": row.ICL_LOT_NUM,
                "IC_COMPOUND_ID": ic_compound_id,
                "action": action,
                "reason": reason,
                "keep_ICL_ID": keep_icl_id,
            })

    return pd.DataFrame(plan_rows).sort_values(
        ["platewell", "action", "ICL_ID"],
        key=lambda s: s.map({"KEEP": 0, "DELETE": 1, "MANUAL": 2}) if s.name == "action" else s,
    ).reset_index(drop=True)


def reconcile_platewell_duplicates(
    engine: Engine,
    source: pd.DataFrame,
    collection: str = "CSC TargetMol Bioactive",
    red_flags_path: Path = DEFAULT_RED_FLAGS,
    execute: bool = False,
) -> pd.DataFrame:
    """Reconcile remaining platewell duplicates using INV_COMPOUND.IC_COMPOUND_ID.

    This is deliberately a separate post-cleanup stage:
      - it only sees LOT_IDs still present in INV_COMP_LOT;
      - it does not alter ICL_LOT_NUM;
      - rerunning after a successful execution produces no DELETE plan for
        those groups because the duplicate has already been resolved;
      - unresolved groups are reported as MANUAL and are never deleted.
    """
    source = validate_platewell_duplicate_source(source)
    flags = load_red_flags(red_flags_path)

    with engine.connect() as conn:
        plan = get_platewell_duplicate_plan(conn, source, collection)

    if plan.empty:
        print("\n=== PLATEWELL DUPLICATE RECONCILIATION PLAN ===")
        print("No remaining platewell duplicates detected.")
        return plan

    manual = plan[plan["action"] == "MANUAL"].copy()
    deletes = plan[plan["action"] == "DELETE"].copy()
    keeps = plan[plan["action"] == "KEEP"].copy()

    print("\n=== PLATEWELL DUPLICATE RECONCILIATION PLAN ===")
    print(f"Duplicate platewell groups: {plan['platewell'].nunique():,}")
    print(f"Lots to keep as canonical IDs: {len(keeps):,}")
    print(f"Lots to delete: {len(deletes):,}")
    print(f"Lots requiring manual reconciliation: {len(manual):,}")

    if not manual.empty:
        print("\n=== MANUAL PLATEWELL RECONCILIATION WARNINGS ===")
        print(
            manual[
                [
                    "platewell",
                    "ICL_ID",
                    "ICL_IC_ID",
                    "IC_COMPOUND_ID",
                    "reason",
                ]
            ].drop_duplicates().to_string(index=False)
        )

    if deletes.empty:
        print("No automatically reconciled platewell duplicates need deletion.")
        return plan

    if not execute:
        print("\nDRY RUN: no platewell-duplicate rows changed.")
        print("Proposed platewell duplicate deletions:")
        print(
            deletes[
                [
                    "platewell",
                    "ICL_ID",
                    "ICL_IC_ID",
                    "IC_COMPOUND_ID",
                    "keep_ICL_ID",
                ]
            ].to_string(index=False)
        )
        return plan

    metadata = MetaData()
    inv_comp_lot = Table(
        "INV_COMP_LOT",
        metadata,
        autoload_with=engine,
    )

    with engine.begin() as conn:
        # Re-plan inside the transaction so a rerun or concurrent cleanup cannot
        # make this stage act on a stale duplicate set.
        plan = get_platewell_duplicate_plan(conn, source, collection)
        deletes = plan[plan["action"] == "DELETE"].copy()

        if deletes.empty:
            print("No platewell duplicate deletions remain when transaction began; skipping.")
            return plan

        delete_ids = sorted(deletes["ICL_ID"].astype(int).unique().tolist())
        preflight_dependencies(conn, flags, delete_ids)
        delete_direct_dependencies(conn, flags, delete_ids)

        before = int(
            conn.execute(
                text(
                    "SELECT COUNT(*) FROM INV_COMP_LOT WHERE ICL_ID IN :ids"
                ).bindparams(bindparam("ids", expanding=True)),
                {"ids": delete_ids},
            ).scalar_one()
        )

        if before != len(delete_ids):
            # The transaction re-plan should normally make this impossible.
            # Treat any race as an explicit abort rather than partially deleting.
            raise RuntimeError(
                "ABORT: platewell duplicate target set changed before delete; "
                f"expected {len(delete_ids)} rows, found {before}."
            )

        result = conn.execute(
            delete(inv_comp_lot).where(inv_comp_lot.c.ICL_ID.in_(delete_ids))
        )
        deleted = int(result.rowcount or 0)
        if deleted != len(delete_ids):
            raise RuntimeError(
                "ABORT inside transaction: expected to delete "
                f"{len(delete_ids)} platewell duplicate lots, deleted {deleted}."
            )

        remaining = int(
            conn.execute(
                text(
                    "SELECT COUNT(*) FROM INV_COMP_LOT WHERE ICL_ID IN :ids"
                ).bindparams(bindparam("ids", expanding=True)),
                {"ids": delete_ids},
            ).scalar_one()
        )
        if remaining:
            raise RuntimeError(
                f"ABORT inside transaction: {remaining} platewell duplicate rows remain."
            )

        final_plan = get_platewell_duplicate_plan(conn, source, collection)
        remaining_auto = final_plan[final_plan["action"] == "DELETE"]
        if not remaining_auto.empty:
            raise RuntimeError(
                "ABORT inside transaction: automatically reconciled platewell duplicates "
                "remain after deletion:\n"
                + remaining_auto.to_string(index=False)
            )

        print(f"Platewell duplicate INV_COMP_LOT rows deleted: {deleted:,}")

    print("Platewell duplicate reconciliation transaction successful.")
    return final_plan

def get_new_source(
    conn,
    collection: str,
    source: pd.DataFrame | None = None,
    keep_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Pull the authoritative surviving INV_COMP_LOT rows for a collection.

    If source + keep_cols are supplied, merge registration metadata back
    using the authoritative identity:

        source.lot_id == INV_COMP_LOT.ICL_ID

    Database columns remain authoritative. Source metadata is retained
    only for columns requested in keep_cols.
    """

    stmt = text(
        """
        SELECT
            ICL_ID,
            ICL_IC_ID,
            ICL_ALIAS,
            ICL_VEN_COLLECTION,
            ICL_LOT_NUM
        FROM INV_COMP_LOT
        WHERE ICL_VEN_COLLECTION = :collection
        ORDER BY
            ICL_IC_ID,
            ICL_ALIAS,
            ICL_LOT_NUM,
            ICL_ID
        """
    )

    result = conn.execute(
        stmt,
        {"collection": collection},
    )

    new_source = pd.DataFrame(result.mappings().all())

    if new_source.empty:
        raise RuntimeError(
            f"No rows found in INV_COMP_LOT for {collection!r}."
        )

    new_source = new_source.rename(
        columns={
            "ICL_ID": "lot_id",
            "ICL_ALIAS": "alias",
        }
    )

    # ------------------------------------------------------------
    # Optionally restore source-file metadata by exact LOT_ID.
    # ------------------------------------------------------------

    if source is not None and keep_cols:

        missing_cols = sorted(
            set(keep_cols) - set(source.columns)
        )

        if missing_cols:
            raise RuntimeError(
                "Source is missing requested keep_cols: "
                f"{missing_cols}"
            )

        if "lot_id" not in keep_cols:
            raise RuntimeError(
                "keep_cols must include 'lot_id' so source metadata "
                "can be merged to INV_COMP_LOT."
            )

        metadata = source[keep_cols].copy()

        if metadata["lot_id"].duplicated().any():
            dup_ids = (
                metadata.loc[
                    metadata["lot_id"].duplicated(keep=False),
                    "lot_id",
                ]
                .drop_duplicates()
                .tolist()
            )

            raise RuntimeError(
                "ABORT: source contains duplicate lot_id values, so "
                "source metadata cannot be safely merged by LOT_ID. "
                f"First IDs: {dup_ids[:20]}"
            )

        # Do NOT allow source metadata to overwrite database columns.
        metadata_cols = [
            c
            for c in metadata.columns
            if c != "lot_id"
            and c not in new_source.columns
        ]

        new_source = new_source.merge(
            metadata[["lot_id"] + metadata_cols],
            on="lot_id",
            how="left",
            validate="one_to_one",
        )

        # --------------------------------------------------------
        # Any surviving DB row without source metadata is suspicious.
        # Do not guess plate/well information.
        # --------------------------------------------------------

        missing_metadata = new_source[
            new_source[metadata_cols].isna().all(axis=1)
        ]

        if not missing_metadata.empty:
            print(
                f"WARNING: {len(missing_metadata):,} surviving "
                f"INV_COMP_LOT rows have no matching source LOT_ID "
                f"metadata."
            )
            print(
                missing_metadata[
                    ["lot_id", "ICL_IC_ID", "alias"]
                ]
                .head(20)
                .to_string(index=False)
            )

    return new_source.reset_index(drop=True)


def delete_extra_lots(
    engine: Engine,
    source: pd.DataFrame,
    collection: str = "CSC TargetMol Bioactive",
    red_flags_path: Path = DEFAULT_RED_FLAGS,
    execute: bool = False,
    keep_cols: list = "platewell",
) -> pd.DataFrame:
    """Plan or execute cleanup of failed collection registration lots.

    The current logic uses keep_num_lots derived from true duplicates;
    keep_first is not a deletion selector. Set keep_first=True for the
    current recommended behavior.

    source: DataFrame with columns lot_id, smdc, alias, keep_num_lots
    red_flags_path: Path to red_flags.csv for dependency actions
    """
    source = validate_source(source)
    source = add_keep_num_lots(source)
    flags = load_red_flags(red_flags_path)

    print("\n=== SOURCE RETENTION PLAN ===")
    print(f"Source rows: {len(source):,}")
    print(f"Source unique LOT_IDs: {source['lot_id'].nunique():,}")
    print(f"Rows with keep_num_lots > 1: {(source['keep_num_lots'] > 1).sum():,}")
    print(f"Distinct (smdc, alias) groups with keep_num_lots > 1: "
          f"{source.loc[source['keep_num_lots'] > 1, ['smdc', 'alias']].drop_duplicates().shape[0]:,}")

    with engine.connect() as conn:
        already_clean = check_already_clean(conn, source, collection, keep_cols)

        if already_clean is not None:
            print(
                "\nAlready-clean database detected for the retention/reset stage. "
                "No retention delete or lot-number reset will be repeated."
            )
            new_source = already_clean
            plan = None
        else:
            source_db = source_to_db_targets(conn, source, collection)

            if source_db.empty:
                print(
                    "No source LOT_IDs remain in the database; treating the "
                    "registration as already reconciled and continuing."
                )
                plan = pd.DataFrame()
                delete_plan = pd.DataFrame()
                keep_plan = pd.DataFrame()
                new_source = get_new_source(conn, collection, source, keep_cols)
            else:
                plan = get_retention_plan(conn, source_db, collection)
                delete_plan = plan[plan["action"] == "DELETE"].copy()
                keep_plan = plan[plan["action"] == "KEEP"].copy()

            if delete_plan.empty:
                print(f"No {collection} rows need deletion.")
                new_source = get_new_source(conn, collection, source, keep_cols)
                # Do not rerun retention/reset work, but continue to the distinct
                # post-cleanup platewell duplicate stage below.
                plan = plan
                delete_plan = pd.DataFrame()
                keep_plan = pd.DataFrame()
            else:
                delete_ids = sorted(delete_plan["ICL_ID"].astype(int).unique().tolist())
                keep_ids = sorted(keep_plan["ICL_ID"].astype(int).unique().tolist())

                print("\n=== DATABASE RETENTION PLAN ===")
                print(f"collection rows considered: {len(plan):,}")
                print(f"{collection} rows to keep: {len(keep_ids):,}")
                print(f"{collection} rows to delete: {len(delete_ids):,}")
                print(f"True-duplicate source LOT_ID rows: {(source['keep_num_lots'] > 1).sum():,}")

                preflight_dependencies(conn, flags, delete_ids)

                reset_plan = get_reset_plan(conn, keep_plan, collection, delete_ids)

                print("\n=== PROPOSED DELETIONS (first 10) ===")
                print(
                    delete_plan[[
                        "ICL_ID", "ICL_IC_ID", "ICL_ALIAS", "ICL_LOT_NUM", "keep_num_lots", "rank_newest"
                    ]].head(10).to_string(index=False)
                )

                print("\n=== PROPOSED LOT NUMBER RESETS (first 10) ===")
                print(
                    reset_plan[["ICL_ID", "ICL_IC_ID", "ICL_ALIAS",  "ICL_LOT_NUM", "new_lot_num"
                    ]].head(10).to_string(index=False)
                )

                if not execute:
                    print("\nDRY RUN: no database rows changed.")
                    # The post-cleanup stage must still be planned independently.
                    platewell_plan = reconcile_platewell_duplicates(
                        engine, source, collection, red_flags_path, execute=False
                    )
                    return plan

    # When the retention/reset stage was already clean, or has nothing to do,
    # jump directly to the distinct post-cleanup platewell reconciliation stage.
    if plan is None or delete_plan.empty:
        platewell_plan = reconcile_platewell_duplicates(
            engine, source, collection, red_flags_path, execute=execute
        )
        if execute:
            with engine.connect() as conn:
                new_source = get_new_source(conn, collection, source, keep_cols)
        return new_source

    metadata = MetaData()
    inv_comp_lot = Table(
        "INV_COMP_LOT",
        metadata,
        autoload_with=engine,
    )

    batch_size=500

    with engine.begin() as conn:
        deleted_children = delete_direct_dependencies(conn, flags, delete_ids)
        for table, n in deleted_children.items():
            print(f"Deleted {n:,} rows from {table}")

        before = conn.execute(
            text(
                "SELECT COUNT(*) FROM INV_COMP_LOT "
                "WHERE ICL_ID IN :ids"
            ).bindparams(bindparam("ids", expanding=True)),
            {"ids": delete_ids},
        ).scalar_one()
        if int(before) != len(delete_ids):
            raise RuntimeError(
                f"ABORT inside transaction: expected {len(delete_ids)} target rows, found {before}."
            )

        connection_info = conn.execute(
            text("""
                        SELECT
                            CONNECTION_ID() AS connection_id,
                            DATABASE() AS db_name,
                            USER() AS db_user
                    """)
        ).mappings().one()

        print("SQLAlchemy connection:", dict(connection_info))
        print(f"About to DELETE {len(delete_ids):,} INV_COMP_LOT rows...")

        total_deleted = 0

        for start in range(0, len(delete_ids), batch_size):
            batch = delete_ids[start:start + batch_size]

            stmt = delete(inv_comp_lot).where(
                inv_comp_lot.c.ICL_ID.in_(batch)
            )

            result = conn.execute(stmt)

            total_deleted += result.rowcount or 0

            print(
                f"Deleted batch {start:,}-{start + len(batch):,}: "
                f"{result.rowcount:,} rows"
            )
        if total_deleted != len(delete_ids):
            raise RuntimeError(
                f"ABORT inside transaction: expected to delete {len(delete_ids)}, deleted {total_deleted}."
            )

        remaining = conn.execute(
            text(
                "SELECT COUNT(*) FROM INV_COMP_LOT "
                "WHERE ICL_ID IN :ids"
            ).bindparams(bindparam("ids", expanding=True)),
            {"ids": delete_ids},
        ).scalar_one()
        if int(remaining) != 0:
            raise RuntimeError(f"ABORT inside transaction: {remaining} deleted rows remain.")

        updated = update_lot_numbers(conn, reset_plan)
        if updated != len(reset_plan):
            raise RuntimeError(
                f"ABORT inside transaction: expected {len(reset_plan)} lot-number updates, got {updated}."
            )

        print(f"INV_COMP_LOT rows deleted: {total_deleted:,}")
        print(f"ICL_LOT_NUM values updated: {updated:,}")

        new_source = get_new_source(
            conn,
            collection,
            source,
            keep_cols
        )

        reconciliation = get_reconciliation_report(
            conn,
            source_db,
            collection,
        )

        bad = reconciliation[
            reconciliation["status"] != "OK"
            ]

        if not bad.empty:
            raise RuntimeError(
                "Cleanup completed but reconciliation failed:\n"
                + bad.to_string(index=False)
            )

        print("Reconciliation status: ", reconciliation["status"].value_counts(dropna=False))
        print("All checks passed.")

    print("Transaction successful.")

    # IMPORTANT: this stage is intentionally after retention and lot-number reset
    # work and is independently idempotent. It never recalculates ICL_LOT_NUM.
    reconcile_platewell_duplicates(
        engine, source, collection, red_flags_path, execute=execute
    )

    # Refresh the authoritative source after the post-cleanup deletion.
    if execute:
        with engine.connect() as conn:
            new_source = get_new_source(conn, collection, source, keep_cols)

    return new_source
# VERSION_MARKER: delete_extra_lots_warn_v2 - missing LOT_IDs warn and continue; no SMDC fallback