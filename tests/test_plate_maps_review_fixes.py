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


# --- plate-map fixtures -------------------------------------------------------
import io
import string

ROWS = list(string.ascii_uppercase[:16])


def _sheet(blocks, separator='header', blanks=1):
    """A 384-well map sheet as read_excel returns it. blocks: list of fill(r, c)."""
    lines = []
    for i, fill in enumerate(blocks):
        if i:
            lines += [[None] * 25] * blanks
            if separator == 'header':
                lines.append([f'block{i}'] + list(range(1, 25)))
        lines += [[r] + [fill(r, c) for c in range(1, 25)] for r in ROWS]
    buf = io.BytesIO()
    pd.DataFrame([['IDs'] + list(range(1, 25))] + lines).to_excel(buf, header=False, index=False)
    buf.seek(0)
    return pd.read_excel(buf, engine='openpyxl')


def test_blocks_with_headers_and_extra_empty_rows_read_correctly():
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: float(c)], blanks=3)
    conc = pm.read_platemap_block(sheet, block=1, value_name='conc')
    assert conc.row.iloc[0] == 'A' and conc.conc.iloc[0] == 1.0
    assert sorted(conc.row.unique()) == ROWS


def test_blank_only_separator_raises_instead_of_shifting_rows():
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: float(c)], separator='blank')
    with pytest.raises(ValueError, match='Block 1'):
        pm.read_platemap_block(sheet, block=1, value_name='conc')


def test_control_dict_keys_are_canonical():
    assert pm._normalize_control_dict({'sci_poscontrol': 'No Enzyme', 'Sci_NegControl': ' DMSO '}) == \
        {'Sci_PosControl': 'noenzyme', 'Sci_NegControl': 'dmso'}


def test_destination_info_flags_and_nulls_controls_in_any_case():
    def ids(r, c):
        return 'NO ENZYME' if c == 1 else ('dmso ' if c == 2 else 1000 + c)
    dp = {'P1': _sheet([ids])}
    row = pd.Series({'Grid_num': 'P1', 'Raw_assay_ID': 1})
    controls = pm._normalize_control_dict({'Sci_PosControl': 'No Enzyme', 'Sci_NegControl': 'DMSO'})
    replace = {'No Enzyme': float('nan'), 'DMSO': float('nan')}
    dids, _ = pm.get_destination_plate_info(dp, row, replace, controls)
    a01 = dids[dids['Destination well'] == 'A01'].iloc[0]
    a02 = dids[dids['Destination well'] == 'A02'].iloc[0]
    assert a01.Sci_PosControl and not a01.Sci_NegControl and pd.isna(a01.SMDC_ID)
    assert a02.Sci_NegControl and pd.isna(a02.SMDC_ID)
    assert dids.Sci_PosControl.sum() == 16 and dids.Sci_NegControl.sum() == 16


def test_final_control_pass_skips_biomek_worklists():
    dids = pd.DataFrame({'Destination well': ['A01'], 'Sci_PosControl': [False], 'Sci_NegControl': [False]})
    row = pd.Series({'Grid_num': 'W1', 'Source_plate': 'Biomek'})
    out = pm._apply_control_flags_from_map(dids, {'W1': pd.DataFrame({'a': [1]})}, row, {'Sci_PosControl': 'x'})
    assert out.equals(dids)


def test_final_control_pass_reads_grid_maps():
    dp = {'P1': _sheet([lambda r, c: 'Pos' if (r, c) == ('B', 3) else 1000 + c])}
    dids = pd.DataFrame({'Destination well': ['B03', 'A01'], 'Sci_PosControl': [False, True],
                         'Sci_NegControl': [False, False]})
    row = pd.Series({'Grid_num': 'P1', 'Source_plate': 'use lots'})
    out = pm._apply_control_flags_from_map(dids, dp, row, {'Sci_PosControl': ' pos'})
    assert out.Sci_PosControl.tolist() == [True, False]


def test_id_plate_pairs_come_from_the_same_wells():
    dids = pd.DataFrame({'SMDC_ID': [5.0, 6.0, 5.0, 7.0, None],
                         'Source plate': ['PA', 'PB', 'PA', 'NA', 'PC']})
    assert pm._id_plate_pairs(dids) == (['5', '6'], ['PA', 'PB'])
    assert pm._id_plate_pairs(pd.DataFrame({'SMDC_ID': [1.0]})) == ([], [])
