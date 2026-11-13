"""Offline tests for plate_maps.lookup_icl_id_by_plate_and_ic_id.

Runs against an in-memory SQLite database that mimics the three INV_* tables
(IP_BARCODE is NOCASE to mimic MySQL's default collation), so no HiTS connection
is needed. Run with:  pytest tests/test_plate_maps_lookup.py
"""
import sys
import types
from pathlib import Path

# Stub only what is missing, so the module imports without a database/extras.
try:
    import mysql_cnx  # noqa: F401
except Exception:
    sys.modules['mysql_cnx'] = types.SimpleNamespace(create_mysql_cnx=lambda *a, **k: None)
try:
    import humanize  # noqa: F401
except ModuleNotFoundError:
    sys.modules['humanize'] = types.SimpleNamespace(precisedelta=lambda *a, **k: '')
try:
    import seaborn  # noqa: F401
except ModuleNotFoundError:
    sys.modules['seaborn'] = types.SimpleNamespace(set_theme=lambda *a, **k: None, heatmap=lambda *a, **k: None)

import numpy as np
import pandas as pd
import pytest
from sqlalchemy import create_engine

import plate_maps as pm

def make_engine(url='sqlite://'):
    e=create_engine(url)
    with e.begin() as c:
        c.exec_driver_sql("create table INV_PLATE(IP_ID integer, IP_BARCODE text collate nocase)")
        c.exec_driver_sql("create table INV_COMP_LOT(ICL_ID integer, ICL_IC_ID integer, ICL_LOT_NUM integer)")
        c.exec_driver_sql("create table INV_WELL(IW_ID integer, IW_IP_ID integer, IW_COORDINATES text, IW_ICL_ID integer)")
        c.exec_driver_sql("insert into INV_PLATE values (1,'2024_NEK'),(2,'NEK_R1_P0-MA01'),(3,'0123')")
        c.exec_driver_sql("insert into INV_COMP_LOT values (500,1191231,1),(501,1191226,1),(502,1191333,1),(503,1191333,2),(504,77,1)")
        c.exec_driver_sql("insert into INV_WELL values (1,1,'N13',500),(2,1,'N18',501),(3,2,'B05',502),(4,2,'B06',503),(5,3,'A01',504)")
    return e

def test_basic_matched_pairs_and_float_ids():
    e=make_engine()
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191231.0,'1191226'],['2024_nek','2024_NEK'],cnx=e)
    assert dict(zip(r.SMDC_ID,r.LOT_ID))=={1191231:500,1191226:501}
    assert set(r['Source plate'])=={'2024_NEK'}

def test_ambiguous_lot_needs_lot_or_well():
    e=make_engine()
    with pytest.raises(ValueError,match='more than one lot'):
        pm.lookup_icl_id_by_plate_and_ic_id([1191333],['NEK_R1_P0-MA01'],cnx=e)
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191333],['NEK_R1_P0-MA01'],lots=[2],cnx=e); assert list(r.LOT_ID)==[503]
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191333],['NEK_R1_P0-MA01'],wells=['b05'],cnx=e); assert list(r.LOT_ID)==[502]
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191333],['NEK_R1_P0-MA01'],strict=False,cnx=e); assert sorted(r.LOT_ID)==[502,503]

def test_missing_raises_or_warns_and_replace_vals():
    e=make_engine()
    with pytest.raises(ValueError,match='no matching well'):
        pm.lookup_icl_id_by_plate_and_ic_id([1191231,999],['2024_NEK','2024_NEK'],cnx=e)
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191231,999],['2024_NEK','2024_NEK'],replace_vals=[999],cnx=e); assert list(r.LOT_ID)==[500]
    r=pm.lookup_icl_id_by_plate_and_ic_id([1191231,999],['2024_NEK','2024_NEK'],strict=False,cnx=e); assert list(r.LOT_ID)==[500]
    with pytest.raises(ValueError,match='no matching well'):   # right compound, wrong plate
        pm.lookup_icl_id_by_plate_and_ic_id([1191231],['NEK_R1_P0-MA01'],cnx=e)
    with pytest.raises(ValueError,match='no matching well'):   # NaN plate
        pm.lookup_icl_id_by_plate_and_ic_id([1191231],[np.nan],cnx=e)

def test_leading_zero_plate_kept_and_length_check_and_empty():
    e=make_engine()
    assert list(pm.lookup_icl_id_by_plate_and_ic_id([77],['0123'],cnx=e).LOT_ID)==[504]
    with pytest.raises(ValueError,match='zipped'):
        pm.lookup_icl_id_by_plate_and_ic_id([1,2],['a'],cnx=e)
    r=pm.lookup_icl_id_by_plate_and_ic_id([5],['x'],replace_vals=[5],cnx=e)
    assert r.empty and list(r.columns)==['SMDC_ID','Lot','Source plate','Source well','LOT_ID']

def test_dr_file_roundtrip_matches_docstring_recipe():
    e=make_engine()
    df=pd.DataFrame({'SMDC_ID':[1191231,1191333,1191333],'Lot':[1,1,2],'Source plate':['2024_NEK','NEK_R1_P0-MA01','NEK_R1_P0-MA01'],'Source well':['N13','B05','B06']})
    res=pm.lookup_icl_id_by_plate_and_ic_id(df['SMDC_ID'],df['Source plate'],lots=df['Lot'],wells=df['Source well'],cnx=e)
    df['_plate']=df['Source plate'].str.upper()
    out=df.merge(res,left_on=['SMDC_ID','Lot','_plate','Source well'],right_on=['SMDC_ID','Lot','Source plate','Source well'],how='left',suffixes=('','_db')).drop(columns=['_plate','Source plate_db'])
    assert list(out.LOT_ID)==[500,502,503] and len(out)==3


def test_lookups_are_parameterised(monkeypatch, tmp_path):
    """A barcode with a quote in it must be treated as data, not SQL."""
    # File-backed, because these helpers dispose their engine after each call.
    e = make_engine(f"sqlite:///{tmp_path / 'hits.db'}")
    with e.begin() as c:
        c.exec_driver_sql("insert into INV_PLATE values (9,'O''BRIEN-01')")
        c.exec_driver_sql("insert into INV_WELL values (9,9,'A01',504)")
    monkeypatch.setattr(pm, '_connect', lambda: e)
    out = pm.lookup_id_and_lot_by_plate_and_well(plates=["O'BRIEN-01"], wells=['A01'])
    assert out[['SMDC_ID', 'Lot', 'Source well']].values.tolist() == [[77, 1, 'A01']]
    out = pm.lookup_plate_by_lot_num(ids=['77'], lot_nums=['1'], like_filter=None)
    assert set(out['Source plate'].str.upper()) == {'0123', "O'BRIEN-01"}
    out = pm.lookup_lot_by_plate_and_ic_id(ids=[77], plates=["O'BRIEN-01"], matched=True)
    assert out['Source well'].tolist() == ['A01']


def test_read_db_defaults_to_production(monkeypatch):
    calls = {}
    sqlalchemy = __import__('sqlalchemy')

    def fake(**k):
        calls.update(k)
        return sqlalchemy.create_engine('sqlite://')
    monkeypatch.setattr(pm, 'create_mysql_cnx', fake)
    pm._connect()
    assert calls == {'db': pm.READ_DB}
    assert pm.READ_DB in ('QB3_HiTS', __import__('os').environ.get('hits_db'))


def test_connect_fails_loudly_when_hits_is_unreachable(monkeypatch, tmp_path):
    # A path in a missing folder makes sqlite fail on first connect, the way a
    # dead SSH tunnel makes MySQL fail: create_engine itself succeeds.
    sqlalchemy = __import__('sqlalchemy')
    bad = sqlalchemy.create_engine(f"sqlite:///{tmp_path}/no/such/dir/x.db")
    monkeypatch.setattr(pm, 'create_mysql_cnx', lambda **k: bad)
    with pytest.raises(ConnectionError, match='SSH tunnel'):
        pm._connect()


def test_db_error_in_lot_lookup_is_not_turned_into_na(monkeypatch):
    from sqlalchemy.exc import OperationalError

    def dead(**k):
        raise OperationalError('select', {}, Exception('Lost connection'))
    monkeypatch.setattr(pm, 'lookup_plate_by_lot_num', dead)
    dids = pd.DataFrame({'SMDC_ID': [5.0], 'Lot': [1.0], 'Destination well': ['A01']})
    row = pd.Series({'Source_plate': 'use lots', 'Grid_num': 'P1'})
    lots = pd.DataFrame({'row': ['A'], 'col': ['01'], 'Lot': [1.0]})
    monkeypatch.setattr(pm, 'read_platemap_block', lambda **k: lots.copy())
    with pytest.raises(OperationalError):
        pm.get_source_plate_info(row, dids.drop(columns='Lot'), None, None, [])


def test_unreachable_hits_in_lot_lookup_is_not_turned_into_na(monkeypatch):
    # _connect() turns a dead tunnel into ConnectionError; that must stop the
    # run too, not fall through to the 'NA' fallback.
    def dead(**k):
        raise ConnectionError('Cannot reach HiTS. Is the SSH tunnel up?')
    monkeypatch.setattr(pm, 'lookup_plate_by_lot_num', dead)
    dids = pd.DataFrame({'SMDC_ID': [5.0], 'Destination well': ['A01']})
    row = pd.Series({'Source_plate': 'use lots', 'Grid_num': 'P1'})
    lots = pd.DataFrame({'row': ['A'], 'col': ['01'], 'Lot': [1.0]})
    monkeypatch.setattr(pm, 'read_platemap_block', lambda **k: lots.copy())
    with pytest.raises(ConnectionError, match='SSH tunnel'):
        pm.get_source_plate_info(row, dids, None, None, [])


def test_lot_lookup_by_user_needs_users():
    with pytest.raises(ValueError, match='users='):
        pm.lookup_lot_by_user(ids=['5'], users=None)
