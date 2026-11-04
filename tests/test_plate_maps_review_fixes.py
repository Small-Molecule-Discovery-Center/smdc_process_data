"""Regression tests for the PR #3 review findings in older plate_maps code."""
import os

import pandas as pd
import pytest

import plate_maps as pm


def test_1536_rows_aa_af_stay_separate(tmp_path, monkeypatch):
    wells = ['A01', 'AA01', 'B01', 'AF02']
    db = pd.DataFrame({
        'IW_COORDINATES': wells, 'Vendor_Alias': ['v1', 'v2', 'v3', 'v4'],
        'SMDC_ID': [1, 2, 3, 4], 'Common_Name': ['n1', 'n2', 'n3', 'n4'],
        'Alias_Type': ['Common_Name'] * 4, 'Lot_Num': [1, 1, 1, 1],
        'Conc_uM': [10] * 4, 'Conc_nM': [None] * 4, 'SMILES': ['C'] * 4,
    })
    monkeypatch.setattr(pm.pd, 'read_sql', lambda *a, **k: db.copy())
    df = pm.make_platemaps_from_barcode('P1', outdir=None, cnx=object())
    assert dict(zip(df.Well, df.Row)) == {'A01': 'A', 'AA01': 'AA', 'B01': 'B', 'AF02': 'AF'}

    pm.make_platemaps_from_barcode('P1', outdir=str(tmp_path), cnx=object())
    sheet = pd.read_excel(tmp_path / 'P1_plate_maps.xlsx', sheet_name='SMDC_ID', index_col=0)
    assert sheet.index.tolist() == ['A', 'B', 'AA', 'AF']      # plate order, not lexical
    assert sheet.loc['A', 1] == 1 and sheet.loc['AA', 1] == 2  # not averaged together


def test_control_flags_added_for_biomek_input():
    out = pm._ensure_control_flags(pd.DataFrame({'SMDC_ID': [1, 2]}))
    assert out.Sci_PosControl.tolist() == [False, False]
    assert out.Sci_NegControl.dtype == bool


def test_control_flags_keep_existing_values():
    out = pm._ensure_control_flags(pd.DataFrame({'Sci_PosControl': [True, None], 'Sci_NegControl': [False, True]}))
    assert out.Sci_PosControl.tolist() == [True, False]
    assert out.Sci_NegControl.tolist() == [False, True]


def test_grid_num_append_keeps_numeric_barcodes_as_text(tmp_path):
    pm._append_grid_num(str(tmp_path), pd.DataFrame({'grid_num': ['000123']}))
    gn = pm._append_grid_num(str(tmp_path), pd.DataFrame({'grid_num': ['000456', '000123']}))
    assert gn.grid_num.tolist() == ['000123', '000456']
    assert (tmp_path / 'grid_num.txt').read_text().split() == ['grid_num', '000123', '000456']


def test_user_lots_same_lot_twice_is_fine():
    lots = pd.DataFrame({'SMDC_ID': [5, 5, 6], 'Lot': [2, 2, 1]})
    assert pm._unique_user_lots(lots).to_dict() == {5: 2, 6: 1}


def test_user_lots_conflict_raises_clearly():
    lots = pd.DataFrame({'SMDC_ID': [5, 5, 6], 'Lot': [2, 3, 1]})
    with pytest.raises(ValueError, match=r'5: \[2, 3\]'):
        pm._unique_user_lots(lots)


def _batch(tmp_path, content):
    """One assay folder with toload_0/20260101/grid_num.txt."""
    toload = tmp_path / '2863_test' / 'toload_0' / '20260101'
    toload.mkdir(parents=True)
    (toload / 'grid_num.txt').write_text(content)
    return tmp_path / '2863_test'


def _info():
    return pd.DataFrame({'Raw_assay_ID': ['2863'], 'Assay_name': ['test'], 'Date_of_screen': ['20260101'],
                         'Raw_data_file': ['x.csv'], 'Grid_num': ['P1']})


def test_clean_up_refuses_to_overwrite_an_earlier_batch(tmp_path, monkeypatch):
    assay = _batch(tmp_path, 'grid_num\nP2\n')
    (assay / 'loaded' / '20260101').mkdir(parents=True)
    (assay / 'loaded' / '20260101' / 'grid_num.txt').write_text('grid_num\nP1\n')
    monkeypatch.setattr(pm, '_read_info_sheet', lambda infile: (None, _info()))
    with pytest.raises(FileExistsError, match='grid_num.txt'):
        pm.clean_up_folders('info.xlsx', basedir=str(tmp_path))
    assert (assay / 'toload_0' / '20260101' / 'grid_num.txt').exists()          # nothing moved
    assert (assay / 'loaded' / '20260101' / 'grid_num.txt').read_text() == 'grid_num\nP1\n'


def test_clean_up_moves_and_allows_identical_rerun(tmp_path, monkeypatch):
    assay = _batch(tmp_path, 'grid_num\nP1\n')
    (assay / 'loaded' / '20260101').mkdir(parents=True)
    (assay / 'loaded' / '20260101' / 'grid_num.txt').write_text('grid_num\nP1\n')  # half-finished earlier run
    monkeypatch.setattr(pm, '_read_info_sheet', lambda infile: (None, _info()))
    pm.clean_up_folders('info.xlsx', basedir=str(tmp_path))
    assert not (assay / 'toload_0' / '20260101').exists()
    assert (assay / 'loaded' / '20260101' / 'grid_num.txt').read_text() == 'grid_num\nP1\n'
