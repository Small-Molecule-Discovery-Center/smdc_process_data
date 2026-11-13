"""Unit tests for the pure-dataframe logic in inventory_registration.

DB-touching helpers (fetch_plates, fetch_mother_wells, assert_*) are exercised
manually in the demo notebook; these tests cover the reshape/copy/build logic
that needs no database.  Run with:  pytest database_interactions/
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import inventory_registration as ir

MAP_DIR = Path(__file__).resolve().parent


def _mother_wells(parents, size, *, plate_col="IP_GRID_NUM"):
    """Build synthetic mother wells: every parent carries all `size`-well coords."""
    coords = ir.well_coordinates(size)
    rows = []
    iw_id = 1
    for ip_id, name in enumerate(parents, start=1001):
        for coord in coords:
            rows.append(
                {
                    "IW_ID": iw_id,
                    "IW_IP_ID": ip_id,
                    plate_col: name,
                    "IW_COORDINATES": coord,
                    "IW_ICL_ID": iw_id,  # every well lot-bearing
                    "IW_CONC_UM": 1000,
                    "IW_VOL_UL": 20,
                }
            )
            iw_id += 1
    return pd.DataFrame(rows)


def _quadrant_count(profile):
    m = ir.resolve_layout_map(profile=profile)
    return m, sorted(m["source_quadrant"].dropna().astype(str).unique())


# ---------------------------------------------------------------------------
# remap_wells_by_quadrant
# ---------------------------------------------------------------------------

def test_interleave_full_fill():
    parents = ["P1", "P2", "P3", "P4"]
    mothers = _mother_wells(parents, 96)
    mapping, quads = _quadrant_count("interleave")

    out = ir.remap_wells_by_quadrant(
        mothers,
        dest_size=384,
        destination_plate_names=["D1"],
        source_plate_col="IP_GRID_NUM",
        profile="interleave",
        source_plate_order=parents,
        plates_per_destination=4,
    )

    # Every map row (source_well, quadrant) matches exactly one parent well.
    assert len(out) == len(mapping)
    assert set(out["_DEST_PLATE_"]) == {"D1"}
    # IW_WELL_INDEX derives from dest coord via the 384 index map.
    idx = ir.well_index_map(384)
    assert (out["IW_WELL_INDEX"] == out["IW_COORDINATES"].map(idx)).all()
    # Provenance note "<parent>-<well>".
    assert out["IW_NOTES"].str.match(r"P\d-[A-H]\d\d").all()


def test_interleave_partial_final_daughter():
    """13 parents, 4 per destination -> chunks [4,4,4,1]; last daughter partial."""
    parents = [f"P{i}" for i in range(1, 14)]
    mothers = _mother_wells(parents, 96)
    dest_names = ["D1", "D2", "D3", "D4"]
    mapping, quads = _quadrant_count("interleave")
    q1_rows = (mapping["source_quadrant"].astype(str) == quads[0]).sum()

    # This used to raise "Need 16 source plates for 4 destinations, found 13".
    out = ir.remap_wells_by_quadrant(
        mothers,
        dest_size=384,
        destination_plate_names=dest_names,
        source_plate_col="IP_GRID_NUM",
        profile="interleave",
        source_plate_order=parents,
        plates_per_destination=4,
    )

    assert set(out["_DEST_PLATE_"]) == set(dest_names)
    # D1..D3 are fully filled; D4 holds only the single Q1 parent's wells.
    assert (out["_DEST_PLATE_"] == "D1").sum() == len(mapping)
    assert (out["_DEST_PLATE_"] == "D4").sum() == q1_rows


def test_quadrant_384_to_1536_full():
    parents = ["S1", "S2", "S3", "S4"]
    mothers = _mother_wells(parents, 384)
    out = ir.remap_wells_by_quadrant(
        mothers,
        dest_size=1536,
        destination_plate_names=["BIG1"],
        source_plate_col="IP_GRID_NUM",
        profile="quadrant",
        source_plate_order=parents,
        plates_per_destination=4,
    )
    assert len(out) == 1536
    assert out["IW_COORDINATES"].nunique() == 1536


def test_remap_384_to_1536_wrapper_matches_engine():
    parents = ["S1", "S2", "S3", "S4"]
    mothers = _mother_wells(parents, 384, plate_col="IW_BARCODE")
    mothers["IW_IP_ID"] = 9001  # single destination grouping

    wrapped = ir.remap_384_to_1536(
        mothers,
        dest_plate_col="IW_IP_ID",
        source_plate_col="IW_BARCODE",
        source_plate_order=parents,
    )
    direct = ir.remap_wells_by_quadrant(
        mothers,
        dest_size=1536,
        destination_plate_names=["9001"],
        source_plate_col="IW_BARCODE",
        profile="quadrant",
        source_plate_order=parents,
        plates_per_destination=4,
        dest_plate_col="IW_IP_ID",
    )
    assert len(wrapped) == len(direct) == 1536
    assert sorted(wrapped["IW_COORDINATES"]) == sorted(direct["IW_COORDINATES"])


# ---------------------------------------------------------------------------
# fill_missing_wells
# ---------------------------------------------------------------------------

def test_fill_missing_wells_does_not_carry_barcode():
    wells = pd.DataFrame(
        {
            "IW_IP_ID": [7, 7],
            "IW_COORDINATES": ["A01", "A02"],
            "IW_WELL_INDEX": [1, 2],
            "IW_BARCODE": ["BC1", "BC2"],
            "IW_INS_BY": [569, 569],
        }
    )
    padded = ir.fill_missing_wells(wells, plate_col="IW_IP_ID", well_col="IW_COORDINATES", size=96)
    assert len(padded) == 96
    new_rows = padded[~padded["IW_COORDINATES"].isin(["A01", "A02"])]
    assert new_rows["IW_BARCODE"].isna().all()      # padded empties carry no barcode
    assert new_rows["IW_INS_BY"].eq(569).all()      # but do carry insert metadata
    idx = ir.well_index_map(96)
    assert (padded["IW_WELL_INDEX"] == padded["IW_COORDINATES"].map(idx)).all()


# ---------------------------------------------------------------------------
# copy layout
# ---------------------------------------------------------------------------

def test_copy_daughter_plates():
    source = pd.DataFrame(
        {
            "IP_ID": [26367, 26368],
            "IP_GRID_NUM": ["05314TM", "05315TM"],
            "IP_BARCODE": ["05314TM-MA01", "05315TM-MA01"],
            "IP_IPT_ID": [1, 1],
        }
    )
    out = ir.copy_daughter_plates(
        source, barcode_suffix="-Eco2", overrides={"IP_IPT_ID": 9}, ins_by=569
    )
    assert "IP_ID" not in out.columns
    assert out["IP_SOURCE_ID"].tolist() == [26367, 26368]
    assert out["IP_BARCODE"].tolist() == ["05314TM-Eco2", "05315TM-Eco2"]
    assert out["IP_IPT_ID"].eq(9).all()
    assert out["IP_INS_BY"].eq(569).all()
    assert out["IP_INS_DATE"].notna().all()


def test_copy_daughter_wells_overrides_only_lot_bearing():
    source = pd.DataFrame(
        {
            "IW_ID": [1, 2, 3],
            "IW_IP_ID": [1001, 1001, 1001],
            "IW_ICL_ID": [111, np.nan, 333],  # middle well empty
            "IW_VOL_UL": [np.nan, np.nan, np.nan],
            "IW_CONC_UM": [np.nan, np.nan, np.nan],
            "IW_COORDINATES": ["A01", "A02", "A03"],
        }
    )
    out = ir.copy_daughter_wells(
        source, new_ip_id_map={1001: 5001}, vol_ul=5, conc_um=10000, ins_by=569
    )
    assert "IW_ID" not in out.columns
    assert out["IW_SOURCE_ID"].tolist() == [1, 2, 3]
    assert out["IW_IP_ID"].eq(5001).all()
    # Only lot-bearing wells (rows 0 and 2) get vol/conc.
    assert out.loc[[0, 2], "IW_VOL_UL"].eq(5).all()
    assert out.loc[[0, 2], "IW_CONC_UM"].eq(10000).all()
    assert pd.isna(out.loc[1, "IW_VOL_UL"])
    assert pd.isna(out.loc[1, "IW_CONC_UM"])


def test_copy_daughter_wells_raises_on_unmapped_plate():
    source = pd.DataFrame({"IW_ID": [1], "IW_IP_ID": [999], "IW_ICL_ID": [1]})
    with pytest.raises(ValueError):
        ir.copy_daughter_wells(source, new_ip_id_map={1001: 5001})


# ---------------------------------------------------------------------------
# Operation A is now layout-free
# ---------------------------------------------------------------------------

def test_prepare_inventory_frames_rejects_legacy_kwargs():
    src = pd.DataFrame({"Plate": ["P1"], "Well": ["A01"]})
    with pytest.raises(TypeError):
        ir.prepare_inventory_frames(
            src,
            plate_source_col="Plate",  # removed parameter
            well_source_col="Well",  # removed parameter
            plate_defaults={},
            well_defaults={},
            plate_lookup=None,
            size=96,
        )


def test_prepare_inventory_frames_columns_vs_literals():
    src = pd.DataFrame(
        {
            "Plate": ["P1", "P1", "P2"],
            "Well": ["A01", "A02", "A01"],
            "barcode": ["bc1", "bc1", "bc2"],
            "conc": [10.0, 10.0, 5.0],
        }
    )
    plates, wells = ir.prepare_inventory_frames(
        src,
        plate_inputs={
            "IP_GRID_NUM": "Plate",  # source column -> plate id
            "IP_BARCODE": "barcode",  # source column -> per-plate first value
            "IP_IPT_ID": 11,  # literal broadcast
            "IP_SCREENING_CENTER": "SMDC",  # literal (not a source column)
            "IP_ICO_ID": 267,  # both collection keys given -> no lookup
            "IP_COLLECTION": "Coll A",
        },
        well_inputs={
            "IW_COORDINATES": "Well",  # source column
            "IW_CONC_UM": "conc",  # source column pulled through
            "IW_VOL_UL": 20,  # required alongside IW_CONC_UM unless IW_FINAL_CONC_NM is filled
            "IW_CONTAM": 0,  # literal broadcast
        },
        size=96,
        plate_lookup=None,
        require_plate_ids=False,
        fill_missing=False,
    )
    # one plate row per distinct plate identifier
    assert sorted(plates["IP_GRID_NUM"]) == ["P1", "P2"]
    assert dict(zip(plates["IP_GRID_NUM"], plates["IP_BARCODE"])) == {"P1": "bc1", "P2": "bc2"}
    assert set(plates["IP_IPT_ID"]) == {11}
    assert set(plates["IP_SCREENING_CENTER"]) == {"SMDC"}
    # wells pull columns through and broadcast literals
    assert list(wells["IW_COORDINATES"]) == ["A01", "A02", "A01"]
    assert list(wells["IW_CONC_UM"]) == [10.0, 10.0, 5.0]
    assert set(wells["IW_CONTAM"]) == {0}
    assert list(wells["IW_WELL_INDEX"]) == [1, 2, 1]
    # the IP_GRID_NUM source column is retained on wells for attach/padding
    assert "Plate" in wells.columns


def test_prepare_inventory_frames_barcode_suffix_key():
    src = pd.DataFrame(
        {
            "Plate": ["P1", "P2"],
            "Well": ["A01", "A01"],
            "sfx": ["-a", "-b"],
        }
    )
    # literal suffix via the IP_BARCODE_SUFFIX control key
    plates, _ = ir.prepare_inventory_frames(
        src,
        plate_inputs={
            "IP_GRID_NUM": "Plate",
            "IP_BARCODE_SUFFIX": "-384",
            "IP_ICO_ID": 1,
            "IP_COLLECTION": "c",
        },
        well_inputs={"IW_COORDINATES": "Well"},
        size=96,
        require_plate_ids=False,
        fill_missing=False,
    )
    assert dict(zip(plates["IP_GRID_NUM"], plates["IP_BARCODE"])) == {
        "P1": "P1-384",
        "P2": "P2-384",
    }
    # IP_BARCODE_SUFFIX is a control key, not written as a DB column
    assert "IP_BARCODE_SUFFIX" not in plates.columns

    # per-plate suffix pulled from a source column, overriding the parameter
    plates2, _ = ir.prepare_inventory_frames(
        src,
        plate_inputs={
            "IP_GRID_NUM": "Plate",
            "IP_BARCODE_SUFFIX": "sfx",
            "IP_ICO_ID": 1,
            "IP_COLLECTION": "c",
        },
        well_inputs={"IW_COORDINATES": "Well"},
        size=96,
        plate_barcode_suffix="-ignored",
        require_plate_ids=False,
        fill_missing=False,
    )
    assert dict(zip(plates2["IP_GRID_NUM"], plates2["IP_BARCODE"])) == {
        "P1": "P1-a",
        "P2": "P2-b",
    }


def _prepare(src, **well_inputs):
    return ir.prepare_inventory_frames(
        src,
        plate_inputs={"IP_GRID_NUM": "Plate", "IP_BARCODE_SUFFIX": "-x", "IP_ICO_ID": 1, "IP_COLLECTION": "c"},
        well_inputs={"IW_COORDINATES": "Well", **well_inputs},
        size=96,
        require_plate_ids=False,
        fill_missing=True,
    )


def test_prepare_inventory_frames_never_copies_ip_barcode_to_wells():
    # IW_BARCODE is a scanned tube barcode; with none scanned, every well is empty.
    src = pd.DataFrame({"Plate": ["P1", "P1", "P2"], "Well": ["A01", "A02", "A01"]})
    plates, wells = _prepare(src, IW_BARCODE=None)
    assert set(plates["IP_BARCODE"]) == {"P1-x", "P2-x"}
    assert wells["IW_BARCODE"].isna().all()
    plates, wells = _prepare(src)
    assert wells["IW_BARCODE"].isna().all()


def test_prepare_inventory_frames_uses_scanned_iw_barcode_column():
    src = pd.DataFrame({"Plate": ["P1", "P1"], "Well": ["A01", "A02"], "IW_BARCODE": ["T001", "T002"]})
    _, wells = _prepare(src)
    real = wells.set_index("IW_COORDINATES").loc[["A01", "A02"], "IW_BARCODE"]
    assert real.tolist() == ["T001", "T002"]
    padded = wells[~wells["IW_COORDINATES"].isin(["A01", "A02"])]
    assert len(padded) == 94 and padded["IW_BARCODE"].isna().all()


def test_copy_daughter_wells_does_not_carry_tube_barcodes():
    source = pd.DataFrame({"IW_ID": [1], "IW_IP_ID": [1001], "IW_ICL_ID": [111], "IW_BARCODE": ["T001"]})
    out = ir.copy_daughter_wells(source, new_ip_id_map={1001: 5001})
    assert out["IW_BARCODE"].isna().all()


def _collection_engine():
    """In-memory SQLite with an INV_COLLECTION(ICO_ID, ICO_NAME) table for lookups."""
    from sqlalchemy import Column, Integer, MetaData, String, Table, create_engine

    engine = create_engine("sqlite://")
    metadata = MetaData()
    table = Table(
        "INV_COLLECTION",
        metadata,
        Column("ICO_ID", Integer, primary_key=True),
        Column("ICO_NAME", String),
    )
    metadata.create_all(engine)
    with engine.begin() as conn:
        conn.execute(
            table.insert(),
            [{"ICO_ID": 267, "ICO_NAME": "Coll A"}, {"ICO_ID": 268, "ICO_NAME": "Coll B"}],
        )
    return engine


def test_prepare_inventory_frames_collection_lookup():
    src = pd.DataFrame({"Plate": ["P1", "P2"], "Well": ["A01", "A01"]})
    engine = _collection_engine()

    # IP_ICO_ID given, IP_COLLECTION looked up from INV_COLLECTION (reflected internally)
    plates, _ = ir.prepare_inventory_frames(
        src,
        plate_inputs={"IP_GRID_NUM": "Plate", "IP_ICO_ID": 267},
        well_inputs={"IW_COORDINATES": "Well"},
        size=96,
        engine=engine,
        require_plate_ids=False,
        fill_missing=False,
    )
    assert set(plates["IP_COLLECTION"]) == {"Coll A"}

    # IP_COLLECTION given, IP_ICO_ID looked up
    plates2, _ = ir.prepare_inventory_frames(
        src,
        plate_inputs={"IP_GRID_NUM": "Plate", "IP_COLLECTION": "Coll B"},
        well_inputs={"IW_COORDINATES": "Well"},
        size=96,
        engine=engine,
        require_plate_ids=False,
        fill_missing=False,
    )
    assert set(plates2["IP_ICO_ID"]) == {268}

    # both supplied -> no lookup, no engine needed
    plates3, _ = ir.prepare_inventory_frames(
        src,
        plate_inputs={"IP_GRID_NUM": "Plate", "IP_ICO_ID": 999, "IP_COLLECTION": "Whatever"},
        well_inputs={"IW_COORDINATES": "Well"},
        size=96,
        require_plate_ids=False,
        fill_missing=False,
    )
    assert set(plates3["IP_ICO_ID"]) == {999}
    assert set(plates3["IP_COLLECTION"]) == {"Whatever"}

    # neither supplied -> raise
    with pytest.raises(ValueError):
        ir.prepare_inventory_frames(
            src,
            plate_inputs={"IP_GRID_NUM": "Plate"},
            well_inputs={"IW_COORDINATES": "Well"},
            size=96,
            require_plate_ids=False,
            fill_missing=False,
        )

    # only one supplied without engine -> raise
    with pytest.raises(ValueError):
        ir.prepare_inventory_frames(
            src,
            plate_inputs={"IP_GRID_NUM": "Plate", "IP_ICO_ID": 267},
            well_inputs={"IW_COORDINATES": "Well"},
            size=96,
            require_plate_ids=False,
            fill_missing=False,
        )


def test_build_plate_rows_ip_barcode_none_defaults_to_grid_num():
    src = pd.DataFrame({"Plate": ["P1", "P2"], "Well": ["A01", "A01"]})
    # explicit None -> barcode from IP_GRID_NUM, then suffix appended (no "None-...")
    plates = ir.build_plate_rows(
        src,
        plate_inputs={
            "IP_GRID_NUM": "Plate",
            "IP_BARCODE": None,
            "IP_BARCODE_SUFFIX": "-MA01",
            "IP_ICO_ID": 1,
            "IP_COLLECTION": "c",
        },
    )
    assert dict(zip(plates["IP_GRID_NUM"], plates["IP_BARCODE"])) == {
        "P1": "P1-MA01",
        "P2": "P2-MA01",
    }


def test_prepare_inventory_frames_rejects_duplicate_source_columns():
    src = pd.DataFrame([["P1", "A01", "x"]], columns=["Plate", "Well", "Plate"])
    with pytest.raises(ValueError, match="duplicate column"):
        ir.prepare_inventory_frames(
            src,
            plate_inputs={"IP_GRID_NUM": "Plate", "IP_ICO_ID": 1, "IP_COLLECTION": "c"},
            well_inputs={"IW_COORDINATES": "Well"},
            size=96,
            require_plate_ids=False,
            fill_missing=True,
        )


def test_prepare_inventory_frames_requires_structural_columns():
    src = pd.DataFrame({"Plate": ["P1"], "Well": ["A01"]})
    # IP_GRID_NUM value that is not a source column must raise
    with pytest.raises(ValueError):
        ir.prepare_inventory_frames(
            src,
            plate_inputs={"IP_GRID_NUM": "does_not_exist"},
            well_inputs={"IW_COORDINATES": "Well"},
            size=96,
            require_plate_ids=False,
        )
    # IW_COORDINATES missing entirely must raise
    with pytest.raises(ValueError):
        ir.prepare_inventory_frames(
            src,
            plate_inputs={"IP_GRID_NUM": "Plate"},
            well_inputs={"IW_CONTAM": 0},
            size=96,
            require_plate_ids=False,
        )


def test_load_layout_map_strips_bom():
    m = ir.load_layout_map(MAP_DIR / "map_96_to_384_interleaved.csv")
    assert "source_well" in m.columns
    assert "﻿source_well" not in m.columns
