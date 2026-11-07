"""Regression tests for the PR #3 review findings in older plate_maps code."""
import os

import numpy as np
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


@pytest.mark.parametrize('separator,blanks', [('blank', 1), ('blank', 3), ('header', 0), ('header', 2)])
def test_every_block_layout_reads_the_right_rows(separator, blanks):
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: float(c) + ROWS.index(r) / 100,
                    lambda r, c: f'L{r}{c}'], separator=separator, blanks=blanks)
    conc = pm.read_platemap_block(sheet, block=1, value_name='conc')
    lots = pm.read_platemap_block(sheet, block=2, value_name='lot')
    a01 = conc[(conc.row == 'A') & (conc.col == '01')].conc.item()
    p24 = conc[(conc.row == 'P') & (conc.col == '24')].conc.item()
    assert (a01, p24) == (1.0, 24.15)
    assert lots[(lots.row == 'A') & (lots.col == '01')].lot.item() == 'LA1'
    assert len(conc) == len(lots) == 384


def test_missing_block_raises_clearly():
    sheet = _sheet([lambda r, c: 1000 + c])
    with pytest.raises(ValueError, match='block 1 does not exist'):
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


def _grouped_info():
    info = pd.DataFrame({
        'Raw_assay_ID': [2863, 2863, 2863], 'Assay_name': ['A'] * 3, 'Date_of_screen': [20260101] * 3,
        'Raw_data_file': ['run.csv'] * 3, 'Grid_num': ['P1', 'P2', 'P3;P4'],
        'Source_plate': ['use map', 'use lots', 'XYZ-MA01'], 'Lot_info': ['use map', '', '']})
    (group, _), = list(pm._iter_load_groups(info))
    return info, group.iloc[0]


def test_each_plate_gets_its_own_info_settings():
    info, row = _grouped_info()
    assert row.Source_plate == 'use map;use lots;XYZ-MA01'   # what the group row carries
    p1 = pm._plate_info_row(info, row, 'P1')
    p2 = pm._plate_info_row(info, row, 'P2')
    p4 = pm._plate_info_row(info, row, 'P4')
    assert (p1.Grid_num, p1.Source_plate, p1.Lot_info) == ('P1', 'use map', 'use map')
    assert (p2.Source_plate, p2.Lot_info) == ('use lots', '')
    assert (p4.Grid_num, p4.Source_plate) == ('P4', 'XYZ-MA01')   # plate listed with another on one Info row
    assert row.Grid_num == 'P1;P2;P3;P4'                            # group row left alone


def test_plate_not_found_falls_back_to_group_row():
    info, row = _grouped_info()
    out = pm._plate_info_row(info, row, 'P9')
    assert out.Grid_num == 'P9' and out.Source_plate == row.Source_plate


def test_folder_pipeline_passes_each_plate_its_own_row(tmp_path, monkeypatch):
    info, _ = _grouped_info()
    share = tmp_path / 'share'
    share.mkdir()
    (share / 'run.csv').write_text('data')
    seen = {}

    def fake_make(dp, row, users, **k):
        seen[row.Grid_num] = row.Source_plate
        return pd.DataFrame({'Destination well': ['A01'], 'Sci_PosControl': [False], 'Sci_NegControl': [False]})
    monkeypatch.setattr(pm, '_read_info_sheet', lambda infile: ({}, info))
    monkeypatch.setattr(pm, 'make_dr_platemaps', fake_make)
    (tmp_path / 'out' / '1_older_assay' / 'toload_0' / '20250101').mkdir(parents=True)
    result = pm.make_dr_platemaps_and_folders('info.xlsx', basedir=str(tmp_path / 'out'), share_root=str(share))
    assert seen == {'P1': 'use map', 'P2': 'use lots', 'P3': 'XYZ-MA01', 'P4': 'XYZ-MA01'}
    folder = next((tmp_path / 'out').glob('2863*/toload_0/20260101'))
    assert pm.load_folders(result) == [str(folder)]
    assert (folder / 'grid_num.txt').read_text().split() == ['grid_num', 'P1', 'P2', 'P3', 'P4']
    assert (folder / 'run.csv').exists()


def _dids():
    return pd.DataFrame({'Destination well': ['A01', 'P24'], 'SMDC_ID': [1001, 1024]})


def test_cherry_pick_map_without_block_1_leaves_concentration_blank():
    sheet = _sheet([lambda r, c: 1000 + c])
    out = pm.get_concentration_info(_dids(), sheet)
    assert out['[compound] uM'].isna().all() and len(out) == 2


def test_cherry_pick_block_1_with_plate_suffix_leaves_concentration_blank():
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: 'MA01' if (r, c) == ('A', 1) else None])
    out = pm.get_concentration_info(_dids(), sheet)
    assert out['[compound] uM'].isna().all()


def test_block_1_concentrations_are_read():
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: 50000 if (r, c) == ('A', 1) else None])
    out = pm.get_concentration_info(_dids(), sheet).set_index('Destination well')
    assert out.loc['A01', '[compound] uM'] == 50000 and pd.isna(out.loc['P24', '[compound] uM'])


def test_block_1_mixing_numbers_and_text_raises():
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: 'oops' if c == 2 else 10.0])
    with pytest.raises(ValueError, match='mixes concentrations with text'):
        pm.get_concentration_info(_dids(), sheet)


@pytest.mark.parametrize('typed', [
    r'Z:\SMDC\Screens\Jinsu\Copy of DR Upload.xlsx',
    'Z:/SMDC/Screens/Jinsu/Copy of DR Upload.xlsx',
    r'\\montreal\Shared\SMDC\Screens\Jinsu\Copy of DR Upload.xlsx',
    '/Volumes/Shared/SMDC/Screens/Jinsu/Copy of DR Upload.xlsx',
    '/mnt/z/SMDC/Screens/Jinsu/Copy of DR Upload.xlsx',
])
def test_share_paths_typed_anywhere_resolve_on_this_machine(tmp_path, typed):
    target = tmp_path / 'SMDC' / 'Screens' / 'Jinsu' / 'Copy of DR Upload.xlsx'
    target.parent.mkdir(parents=True)
    target.write_text('x')
    assert pm.resolve_share_path(typed, share_root=str(tmp_path)) == str(target)


def test_existing_and_non_share_paths_are_left_alone(tmp_path):
    assert pm.resolve_share_path(str(tmp_path), share_root='/nowhere') == str(tmp_path)
    assert pm.resolve_share_path('./', share_root='/nowhere') == './'
    assert pm.resolve_share_path(r'C:\Users\me\x.xlsx', share_root='/nowhere') == r'C:\Users\me\x.xlsx'


def test_folder_pipeline_reads_a_windows_infile_from_wsl(tmp_path, monkeypatch):
    seen = {}
    def fake_read(infile):
        seen['infile'] = infile
        raise RuntimeError('stop')
    monkeypatch.setattr(pm, '_read_info_sheet', fake_read)
    with pytest.raises(RuntimeError, match='stop'):
        pm.make_dr_platemaps_and_folders(infile=r'Z:\SMDC\a.xlsx', basedir=r'Z:\SMDC\toload',
                                         share_root=str(tmp_path))
    assert seen['infile'] == os.path.join(str(tmp_path), 'SMDC/a.xlsx')


def test_copy_survives_a_mount_that_refuses_timestamps(tmp_path, monkeypatch, caplog):
    src = tmp_path / 'a' / 'run.csv'
    src.parent.mkdir()
    src.write_text('data')
    (src.parent / '.DS_Store').write_text('junk')
    (src.parent / 'sub').mkdir()
    (src.parent / 'sub' / 'b.csv').write_text('more')

    def refuse(*a, **k):
        raise PermissionError(1, 'Operation not permitted')
    monkeypatch.setattr(pm.shutil, 'copystat', refuse)

    out = tmp_path / 'out.csv'
    pm._copy_file(str(src), str(out))
    assert out.read_text() == 'data'
    assert 'could not keep the original timestamps' in caplog.text

    dst = tmp_path / 'loaded'
    pm._copy_tree(str(src.parent), str(dst), ignore=pm.ignore_files)
    assert (dst / 'run.csv').read_text() == 'data'
    assert (dst / 'sub' / 'b.csv').read_text() == 'more'
    assert not (dst / '.DS_Store').exists()


def test_copy_keeps_timestamps_when_allowed(tmp_path):
    src = tmp_path / 'run.csv'
    src.write_text('data')
    os.utime(src, (1_000_000_000, 1_000_000_000))
    out = tmp_path / 'out.csv'
    pm._copy_file(str(src), str(out))
    assert int(os.path.getmtime(out)) == 1_000_000_000


def test_load_folders_needs_the_returned_info():
    with pytest.raises(ValueError, match='load folders'):
        pm.load_folders(pd.DataFrame({'Grid_num': ['P1']}))
    info = pd.DataFrame({'Grid_num': ['P1']})
    info.attrs['load_folders'] = ['/x/1_a/toload_0/20260101']
    assert pm.load_folders(info) == ['/x/1_a/toload_0/20260101']


def _well(dids, w):
    return dids[dids['Destination well'] == w].iloc[0]


def test_id_dash_lot_cells_pin_the_lot_and_leave_other_wells_open():
    def ids(r, c):
        return {1: '12346-1', 2: ' 12347 - 3 ', 3: 'No-Enzyme'}.get(c, 1000 + c)
    dp = {'P1': _sheet([ids])}
    row = pd.Series({'Grid_num': 'P1', 'Raw_assay_ID': 1})
    dids, _ = pm.get_destination_plate_info(dp, row, {'No-Enzyme': float('nan')}, {})
    a01, a02, a03, a04 = (_well(dids, w) for w in ['A01', 'A02', 'A03', 'A04'])
    assert (a01.SMDC_ID, a01.Lot) == (12346, 1)
    assert (a02.SMDC_ID, a02.Lot) == (12347, 3)
    assert pd.isna(a03.SMDC_ID)                    # a label with a dash is not a lot
    assert a04.SMDC_ID == 1004 and pd.isna(a04.Lot)  # left for the lot lookups


def test_malformed_id_dash_lot_fails_loudly():
    dp = {'P1': _sheet([lambda r, c: '12346-1-2' if c == 1 else 1000 + c])}
    row = pd.Series({'Grid_num': 'P1', 'Raw_assay_ID': 1})
    with pytest.raises(ValueError, match='12346-1-2'):
        pm.get_destination_plate_info(dp, row, {}, {})


def test_pinned_lots_win_and_use_map_fills_the_rest():
    sheet = _sheet([lambda r, c: '12346-1' if c == 1 else 1000 + c,
                    lambda r, c: 10.0,
                    lambda r, c: 7])
    dids = pd.DataFrame({'Destination well': ['A01', 'A02'], 'SMDC_ID': [12346.0, 1002.0],
                         'Lot': [1.0, np.nan], 'Source plate': ['SP1', 'SP1']})
    row = pd.Series({'Grid_num': 'P1', 'Source_plate': 'SP1', 'Lot_info': 'use map'})
    out = pm.get_lot_info(row, dids, sheet, users=[], replace_vals=[])
    assert _well(out, 'A01').Lot == 1 and _well(out, 'A02').Lot == 7


def test_lot_source_plate_mode_keeps_pins_and_reads_block_2_for_the_rest(monkeypatch):
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: 10.0, lambda r, c: 7])
    dids = pd.DataFrame({'Destination well': ['A01', 'A02'], 'SMDC_ID': [12346.0, 1002.0],
                         'Lot': [1.0, np.nan]})
    seen = {}

    def fake_lookup(ids, lot_nums, **k):
        seen['pairs'] = list(zip(ids, lot_nums))
        return pd.DataFrame({'SMDC_ID': [12346.0, 1002.0], 'Lot': [1.0, 7.0],
                             'Source plate': ['SP1', 'SP2'], 'Source well': ['A01', 'B01']})
    monkeypatch.setattr(pm, 'lookup_plate_by_lot_num', fake_lookup)
    row = pd.Series({'Grid_num': 'P1', 'Source_plate': 'use lots'})
    out = pm.get_source_plate_info(row, dids, sheet, cnx=None, replace_vals=[])
    assert sorted(seen['pairs']) == [('1002.0', '7'), ('12346.0', '1.0')] or \
        sorted(seen['pairs']) == [('1002.0', '7.0'), ('12346.0', '1.0')]
    assert _well(out, 'A01')['Source plate'] == 'SP1' and _well(out, 'A02')['Source plate'] == 'SP2'


def test_all_pinned_lots_need_no_block_2(monkeypatch):
    sheet = _sheet([lambda r, c: 1000 + c, lambda r, c: 10.0])
    dids = pd.DataFrame({'Destination well': ['A01'], 'SMDC_ID': [12346.0], 'Lot': [1.0]})
    monkeypatch.setattr(pm, 'lookup_plate_by_lot_num', lambda ids, lot_nums, **k: pd.DataFrame(
        {'SMDC_ID': [12346.0], 'Lot': [1.0], 'Source plate': ['SP1'], 'Source well': ['A01']}))
    out = pm.get_source_plate_info(pd.Series({'Grid_num': 'P1', 'Source_plate': 'use lots'}), dids, sheet,
                                   cnx=None, replace_vals=[])
    assert _well(out, 'A01')['Source plate'] == 'SP1'


def test_blank_source_plate_is_an_error_not_use_lots():
    dids = pd.DataFrame({'Destination well': ['A01'], 'SMDC_ID': [1001.0]})
    with pytest.raises(ValueError, match='Source_plate is blank for P1'):
        pm.get_source_plate_info(pd.Series({'Grid_num': 'P1', 'Source_plate': np.nan}), dids,
                                 _sheet([lambda r, c: 1000 + c]), cnx=None, replace_vals=[])


def test_folder_pipeline_stops_before_writing_when_source_plate_is_blank(tmp_path, monkeypatch):
    info = pd.DataFrame({'Grid_num': ['P1', 'P2'], 'Raw_assay_ID': [1, 1], 'Assay_name': ['a', 'a'],
                         'Date_of_screen': [20260101, 20260101], 'Raw_data_file': ['r.csv', 'r.csv']})
    monkeypatch.setattr(pm, '_read_info_sheet', lambda infile: ({}, info))
    with pytest.raises(ValueError, match='P1, P2'):
        pm.make_dr_platemaps_and_folders('x.xlsx', basedir=str(tmp_path / 'out'), share_root=str(tmp_path))
    assert not (tmp_path / 'out').exists()
    assert pm.source_plate_barcodes(info) == []


@pytest.mark.parametrize('value,mode', [
    ('', 'blank'), (np.nan, 'blank'),
    ('NA', 'na'), ('n/a', 'na'), ('N.A.', 'na'), ('none', 'na'), ('n.a. (by hand)', 'na'),
    ('use my name', 'na'), ('Use My Name', 'na'), ('use users', 'na'),
    ('use lots', 'lot'), ('Lots', 'lot'), ('use map', 'map'), ('MAP', 'map'),
    ('use biomek / use echo', 'biomek'),
    ('SMDC-MAP-384-01-Eco1', 'barcodes'), ('PLOT12-MA01; X-MA01', 'barcodes'), ('NA-0042', 'barcodes'),
])
def test_source_plate_keywords_are_not_matched_inside_barcodes(value, mode):
    assert pm._source_plate_mode(value) == mode


def test_unknown_use_keyword_is_an_error():
    with pytest.raises(ValueError, match="Unknown Source_plate keyword"):
        pm._source_plate_mode('use plates')


def test_info_sheet_keeps_na_typed_in_source_plate(tmp_path):
    path = tmp_path / 'upload.xlsx'
    info = pd.DataFrame({'Grid_num': ['P1', 'P2', 'P3'], 'Source_plate': ['NA', 'n/a', None],
                         'Raw_assay_ID': [1, 1, 1], 'Assay_name': ['a'] * 3,
                         'Date_of_screen': [20260101] * 3, 'Raw_data_file': ['r.csv'] * 3})
    with pd.ExcelWriter(path) as xw:
        info.to_excel(xw, sheet_name='Info', index=False)
    _, read = pm._read_info_sheet(str(path))
    assert read.Source_plate.iloc[0] == 'NA' and read.Source_plate.iloc[1] == 'n/a'
    assert pd.isna(read.Source_plate.iloc[2])
    with pytest.raises(ValueError, match='blank for P3'):
        pm._check_source_plates(read)


def test_source_plate_barcodes_skips_keywords_but_keeps_map_barcodes():
    info = pd.DataFrame({'Source_plate': ['use lots', 'NA', 'SMDC-MAP-01; B2', 'B2\nB3']})
    assert pm.source_plate_barcodes(info) == ['SMDC-MAP-01', 'B2', 'B3']


def test_source_plate_help_lists_exact_options_in_quotes():
    for option in ('"use lots"', '"use map"', '"use biomek"', '"NA"', '"n/a"', '"n.a."', '"none"', '"use my name"', '";"'):
        assert option in pm.SOURCE_PLATE_HELP


def test_unknown_source_plate_barcode_says_so(monkeypatch):
    empty = pd.DataFrame(columns=['IW_COORDINATES', 'Vendor_Alias', 'SMDC_ID', 'Common_Name',
                                  'Alias_Type', 'Lot_Num', 'Conc_uM', 'Conc_nM', 'SMILES'])
    monkeypatch.setattr(pm.pd, 'read_sql', lambda *a, **k: empty.copy())
    with pytest.raises(ValueError, match="'NOPE-MA01' has no wells in HiTS"):
        pm.make_platemaps_from_barcode('NOPE-MA01', outdir=None, cnx=object())
