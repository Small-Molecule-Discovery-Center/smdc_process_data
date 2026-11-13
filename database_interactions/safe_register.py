"""Safe front door for registering plates and wells in HiTS.

Use this instead of calling ``append_frame`` / ``to_sql`` yourself. It adds
three guard rails on top of ``inventory_registration``:

1. **Test database by default.** ``connect()`` points at ``QB3_HiTS_nightly``
   unless you name ``QB3_HiTS`` (production) yourself.
2. **Dry run by default.** ``register()`` inserts the plates and wells inside
   one transaction, shows you what it did, and then rolls everything back.
   Nothing is saved until you run it again with ``save=True``.
3. **All or nothing.** Plates and wells go in together. If any well fails, the
   plates are rolled back too, so you never get plates with no wells.

Typical notebook use (one settings cell, then one run cell)::

    import safe_register as sr
    import inventory_registration as ir

    DATABASE = "QB3_HiTS_nightly"   # change to "QB3_HiTS" only when ready
    SAVE = False                     # False = dry run, True = really save

    engine = sr.connect(DATABASE)
    plates, wells = ir.prepare_inventory_frames(source, plate_inputs=..., well_inputs=..., size=384)
    sr.register(engine, plates, sr.wells_for_new_plates(wells), save=SAVE)

For daughter plates, build the plates with ``ir.copy_daughter_plates`` and pass
``sr.wells_for_daughter_plates(mother_wells)`` as the second argument.

Connection settings come from environment variables, never from code; see
``mysql_cnx.create_mysql_cnx`` for the list.
"""

from __future__ import annotations

from typing import Any, Callable

import pandas as pd
from sqlalchemy import MetaData, Table, insert, select
from sqlalchemy.engine import Connection, Engine

try:
    import inventory_registration as ir      # run from this folder (and in tests)
except ImportError:
    from database_interactions import inventory_registration as ir   # installed package
from mysql_cnx import create_mysql_cnx

TEST_DB = "QB3_HiTS_nightly"
PRODUCTION_DB = "QB3_HiTS"

WellBuilder = Callable[[pd.DataFrame], pd.DataFrame]


def connect(db: str = TEST_DB) -> Engine:
    """Open a connection to HiTS. Defaults to the nightly test database."""
    return create_mysql_cnx(db=db)


def wells_for_new_plates(wells: pd.DataFrame, *, well_key_col: str = "plate") -> WellBuilder:
    """Well builder for brand-new plates: matches wells to plates by grid number."""

    def build(inserted: pd.DataFrame) -> pd.DataFrame:
        return ir.attach_plate_ids_to_wells(wells, inserted, well_key_col=well_key_col)

    return build


def wells_for_daughter_plates(mother_wells: pd.DataFrame, **copy_kwargs: Any) -> WellBuilder:
    """Well builder for 1:1 daughter plates: copies mother wells onto each daughter.

    Extra keyword arguments (vol_ul, conc_um, ins_by, ...) go to
    ``ir.copy_daughter_wells``.
    """

    def build(inserted: pd.DataFrame) -> pd.DataFrame:
        new_ip_id_map = dict(zip(inserted["IP_SOURCE_ID"], inserted["IP_ID"]))
        return ir.copy_daughter_wells(
            mother_wells,
            new_ip_id_map=new_ip_id_map,
            **copy_kwargs,
        )

    return build


def _insert(conn: Connection, table: Table, frame: pd.DataFrame) -> None:
    rows = frame.astype(object).where(frame.notna(), None).to_dict(orient="records")
    if rows:
        conn.execute(insert(table), rows)


def register(
    engine: Engine,
    plates: pd.DataFrame,
    build_wells: WellBuilder,
    *,
    save: bool = False,
    plate_table_name: str = "INV_PLATE",
    well_table_name: str = "INV_WELL",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Insert plates, then their wells, in one transaction.

    With ``save=False`` (the default) everything is rolled back at the end, so
    the database is unchanged; use it to check that the database accepts the
    rows. With ``save=True`` the rows are committed.

    Returns the inserted plates (with their new IP_IDs) and the well rows.
    During a dry run the IP_IDs are real-looking but are thrown away.
    """
    metadata = MetaData()
    plate_table = Table(plate_table_name, metadata, autoload_with=engine)
    well_table = Table(well_table_name, metadata, autoload_with=engine)
    db = engine.url.database

    plate_rows = ir.prepare_plate_insert_frame(plates, plate_table=plate_table)
    ir.assert_new_barcodes(engine, plate_table, plate_rows["IP_BARCODE"].tolist())

    conn = engine.connect()
    trans = conn.begin()
    try:
        _insert(conn, plate_table, plate_rows)

        # Look plates up by barcode (checked unique above): daughters share
        # their mother's grid number, so grid number alone is ambiguous.
        barcodes = plate_rows["IP_BARCODE"].tolist()
        cols = ["IP_GRID_NUM", "IP_ID", "IP_BARCODE", "IP_SOURCE_ID"]
        stmt = select(*[plate_table.c[c] for c in cols]).where(plate_table.c.IP_BARCODE.in_(barcodes))
        inserted = pd.DataFrame(conn.execute(stmt).fetchall(), columns=cols)
        if len(inserted) != len(plate_rows):
            raise RuntimeError(
                f"Expected {len(plate_rows)} new plates but found {len(inserted)}; nothing was saved."
            )

        well_rows = ir.prepare_well_insert_frame(build_wells(inserted), well_table=well_table)
        _insert(conn, well_table, well_rows)

        if save:
            trans.commit()
            print(f"SAVED to {db}: {len(plate_rows)} plates, {len(well_rows)} wells.")
        else:
            trans.rollback()
            print(
                f"DRY RUN on {db}: {len(plate_rows)} plates and {len(well_rows)} wells would be saved. "
                "Nothing was saved. Set SAVE = True to save."
            )
    except Exception:
        trans.rollback()
        print(f"ERROR: nothing was saved to {db}. See the message below.")
        raise
    finally:
        conn.close()

    return inserted, well_rows
