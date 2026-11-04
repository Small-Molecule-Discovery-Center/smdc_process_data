"""Tests for safe_register using an in-memory SQLite stand-in for HiTS."""

import pandas as pd
import pytest
from sqlalchemy import Column, Float, Integer, MetaData, String, Table, create_engine, func, select

import safe_register as sr


@pytest.fixture
def engine():
    eng = create_engine("sqlite://")
    md = MetaData()
    Table(
        "INV_PLATE", md,
        Column("IP_ID", Integer, primary_key=True, autoincrement=True),
        Column("IP_GRID_NUM", String), Column("IP_BARCODE", String, unique=True),
        Column("IP_IPT_ID", Integer), Column("IP_STATUS_ID", Integer),
        Column("IP_SCREENING_CENTER", Integer), Column("IP_SOURCE_ID", Integer),
    )
    Table(
        "INV_WELL", md,
        Column("IW_ID", Integer, primary_key=True, autoincrement=True),
        Column("IW_IP_ID", Integer, nullable=False), Column("IW_COORDINATES", String),
        Column("IW_WELL_INDEX", Integer, nullable=False), Column("IW_CONC_UM", Float),
        Column("IW_SOURCE_ID", Integer), Column("IW_BARCODE", String),
    )
    md.create_all(eng)
    return eng


def _plates():
    return pd.DataFrame({
        "IP_GRID_NUM": ["G1", "G2"], "IP_BARCODE": ["B1", "B2"],
        "IP_IPT_ID": 1, "IP_STATUS_ID": 1, "IP_SCREENING_CENTER": 1,
    })


def _wells(index=1):
    return pd.DataFrame({
        "plate": ["G1", "G1", "G2"], "IW_COORDINATES": ["A01", "A02", "A01"],
        "IW_WELL_INDEX": [1, 2, index], "IW_CONC_UM": [10.0, 10.0, 10.0],
    })


def _count(engine, name):
    t = Table(name, MetaData(), autoload_with=engine)
    with engine.connect() as c:
        return c.execute(select(func.count()).select_from(t)).scalar()


def test_dry_run_saves_nothing(engine):
    inserted, wells = sr.register(engine, _plates(), sr.wells_for_new_plates(_wells()))
    assert len(inserted) == 2 and len(wells) == 3
    assert _count(engine, "INV_PLATE") == 0
    assert _count(engine, "INV_WELL") == 0


def test_save_commits_plates_and_wells(engine):
    sr.register(engine, _plates(), sr.wells_for_new_plates(_wells()), save=True)
    assert _count(engine, "INV_PLATE") == 2
    assert _count(engine, "INV_WELL") == 3


def test_failed_wells_roll_back_plates(engine):
    with pytest.raises(Exception):
        sr.register(engine, _plates(), sr.wells_for_new_plates(_wells(index=None)), save=True)
    assert _count(engine, "INV_PLATE") == 0
    assert _count(engine, "INV_WELL") == 0


def test_existing_barcode_is_refused(engine):
    sr.register(engine, _plates(), sr.wells_for_new_plates(_wells()), save=True)
    with pytest.raises(ValueError, match="already exist"):
        sr.register(engine, _plates(), sr.wells_for_new_plates(_wells()), save=True)


def test_daughter_plates_copy_mother_wells(engine):
    sr.register(engine, _plates(), sr.wells_for_new_plates(_wells()), save=True)
    with engine.connect() as c:
        mothers = pd.read_sql("select * from INV_WELL", c)
        mother_plates = pd.read_sql("select * from INV_PLATE", c)
    daughters = mother_plates.drop(columns=["IP_ID"]).assign(
        IP_BARCODE=["B1-D", "B2-D"], IP_SOURCE_ID=mother_plates["IP_ID"]
    )
    inserted, wells = sr.register(
        engine, daughters, sr.wells_for_daughter_plates(mothers), save=True
    )
    assert set(inserted["IP_SOURCE_ID"]) == set(mother_plates["IP_ID"])
    assert set(wells["IW_IP_ID"]) == set(inserted["IP_ID"])
    assert _count(engine, "INV_WELL") == 6


def test_connect_defaults_to_nightly(monkeypatch):
    for v in ("hits_host", "hits_db"):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("mysql_user", "u")
    monkeypatch.setenv("mysql80_pwd", "p")
    monkeypatch.setenv("mysql80_port", "3306")
    assert sr.connect().url.database == sr.TEST_DB

