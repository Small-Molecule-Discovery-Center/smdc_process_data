import pytest
"""Smoke tests for the pieces that can run without a database.

Verifies that _normalize_share_path reproduces the ORIGINAL inline logic exactly
(so extracting it changed nothing), that _assay_dir_name is stable across the
int/float/str spellings that used to break clean_up_folders, and that
handle_zero_concentrations drops zeros without prompting.
"""
import sys, types

# Stub only what is missing, so the module loads with no database present
# and without optional extras. None of them are
# touched by the code under test.
try:
    import mysql_cnx  # noqa: F401  (only builds engines lazily; no DB needed)
except Exception:
    sys.modules['mysql_cnx'] = types.SimpleNamespace(create_mysql_cnx=lambda *a, **k: None)
try:
    import humanize  # noqa: F401
except ModuleNotFoundError:
    sys.modules['humanize'] = types.SimpleNamespace(precisedelta=lambda *a, **k: '0 seconds')
try:
    import seaborn  # noqa: F401
except ModuleNotFoundError:
    sys.modules['seaborn'] = types.SimpleNamespace(set_theme=lambda *a, **k: None,
                                                   heatmap=lambda *a, **k: None)

import pandas as pd
import plate_maps as pm


def original_normalize(raw):
    """The inline logic exactly as it was before extraction.

    Kept verbatim, INCLUDING its off-by-one on the 'Montreal/Shared/' branch, so
    the parity check below still proves the extraction was faithful everywhere
    else. That one branch is a deliberate divergence and is asserted separately.
    """
    f = str(raw)
    f = f.replace('\\\\', '/').replace('\\', '/').replace('Z:', '').replace('Y:', '').replace('A:', '').replace('/10.0.0.1/shared/', '')
    if f.startswith('/shared/'):
        f = f[8:]
    elif f.startswith('Montreal/Shared/'):
        f = f[17:]  # the bug: prefix is 16 chars
    elif f.startswith('shared/'):
        f = f[7:]
    elif f.startswith('/Volumes/Shared/'):
        f = f[16:]
    elif f.lower().startswith('smb://montreal.ucsf.edu/shared/'):
        f = f[31:]
    elif f.startswith('smb://montreal/Shared/'):
        f = f[22:]
    elif f.startswith('mb://montreal.ucsf.edu/shared/'):
        f = f[30:]
    elif f.lower().startswith("/montreal.ucsf.edu/shared/"):
        f = f[26:]
    elif f.startswith('/'):
        f = f[1:]
    return f


# The one input whose handling was deliberately changed. Excluded from the
# parity check and asserted on its own below.
FIXED_CASE = 'Montreal/Shared/SMDC/run11.txt'

CASES = [
    r'Z:\SMDC\screens\run1.rpt',
    r'\\montreal.ucsf.edu\shared\SMDC\run2.csv',
    '/shared/SMDC/screens/run3.txt',
    'shared/SMDC/run4.txt',
    '/Volumes/Shared/SMDC/run5.txt',
    'smb://montreal.ucsf.edu/shared/SMDC/run6.txt',
    'SMB://Montreal.ucsf.edu/Shared/SMDC/run7.txt',
    'smb://montreal/Shared/SMDC/run8.txt',
    'mb://montreal.ucsf.edu/shared/SMDC/run9.txt',
    '/montreal.ucsf.edu/shared/SMDC/run10.txt',
    '/10.0.0.1/shared/SMDC/run12.txt',
    'Y:/SMDC/run13.txt',
    '/leading/slash/run14.txt',
]

@pytest.mark.parametrize('case', CASES)
def test_normalize_share_path_matches_original(case):
    assert pm._normalize_share_path(case) == original_normalize(case)


def test_montreal_shared_off_by_one_is_fixed():
    # 'Montreal/Shared/' used to slice [17:] against a 16-character prefix,
    # eating the 'S' of 'SMDC'. If this fails, the off-by-one has come back.
    assert pm._normalize_share_path(FIXED_CASE) == 'SMDC/run11.txt'
    assert original_normalize(FIXED_CASE) == 'MDC/run11.txt'  # documents what was wrong


@pytest.mark.parametrize('prefix', ['Montreal/Shared/', '/shared/', 'shared/', '/Volumes/Shared/',
                                    'smb://montreal.ucsf.edu/shared/', 'smb://montreal/Shared/',
                                    'mb://montreal.ucsf.edu/shared/', '/montreal.ucsf.edu/shared/'])
def test_each_prefix_strips_exactly_itself(prefix):
    assert pm._normalize_share_path(prefix + 'SMDC/x.txt') == 'SMDC/x.txt'


@pytest.mark.parametrize('case', CASES + [FIXED_CASE])
def test_normalized_path_is_relative(case):
    # A leading '/' would make os.path.join discard the share root.
    assert not pm._normalize_share_path(case).startswith('/')


def test_assay_dir_name_agrees_across_spellings():
    assert len({pm._assay_dir_name(r, 'MyAssay') for r in [123, 123.0, '123', '123.0']}) == 1
    assert pm._assay_dir_name('map', 'MyAssay') == 'MyAssay'


def test_handle_zero_concentrations_drops_zeros_without_prompting(monkeypatch):
    monkeypatch.setattr('builtins.input', lambda *a: pytest.fail('handle_zero_concentrations prompted'))
    df = pd.DataFrame({'Destination well': ['A01', 'A02', 'A03', 'A04'],
                       '[compound] uM': [10.0, 1.0, 0.0, 0.0]})
    out = pm.handle_zero_concentrations(df)
    assert out['[compound] uM'].tolist() == [10.0, 1.0]
    assert len(pm.handle_zero_concentrations(pd.DataFrame({'[compound] uM': [10.0, 1.0]}))) == 2


def test_read_platemap_block_rejects_bad_size():
    # A bad size must raise ValueError, not UnboundLocalError or a silent pass.
    with pytest.raises(ValueError, match='1596'):
        pm.read_platemap_block(plate_map=pd.DataFrame(), block=0, size=1596)


def _dr_folder(tmp_path):
    d = tmp_path / '2863_test' / 'toload_0' / '20260101'
    d.mkdir(parents=True)
    (tmp_path / '2863_test' / 'loaded').mkdir()
    pd.DataFrame({'SMDC_ID': [5, 5, 6, 7], 'Lot': [1, 1, 1, None],
                  'Source plate': ['A', 'B', 'A', 'A'],
                  'Source well': ['A01', 'A02', 'A03', 'A04']}).to_excel(d / 'P1_DR_plate.xlsx', index=False)
    pd.DataFrame({'SMDC_ID': [9]}).to_excel(tmp_path / '2863_test' / 'loaded' / 'old_DR_plate.xlsx', index=False)
    return d / 'P1_DR_plate.xlsx'


def test_dr_plate_files_skips_loaded(tmp_path):
    path = _dr_folder(tmp_path)
    assert list(pm.dr_plate_files(str(tmp_path), ['2863'])) == [str(path)]
    assert list(pm.dr_plate_files(str(tmp_path), ['9999'])) == []


def test_dr_plate_files_does_not_match_longer_assay_id(tmp_path):
    path = _dr_folder(tmp_path)
    other = tmp_path / '28630_other' / 'toload_0'
    other.mkdir(parents=True)
    pd.DataFrame({'SMDC_ID': [1]}).to_excel(other / 'X_DR_plate.xlsx', index=False)
    assert list(pm.dr_plate_files(str(tmp_path), ['2863'])) == [str(path)]
    assert list(pm.dr_plate_files(str(tmp_path), [2863])) == [str(path)]


@pytest.mark.parametrize('keep', [[2], [-1], [True], [0.0]])
def test_remove_duplicates_rejects_bad_keep_before_writing(tmp_path, keep):
    path = _dr_folder(tmp_path)
    dups = pm.find_duplicate_source_plates(str(tmp_path), ['2863'])
    with pytest.raises(ValueError, match='keep'):
        pm.remove_duplicate_source_plates(dups, keep=keep, write=True)
    assert len(pd.read_excel(path)) == 4


def test_remove_duplicates_keeps_literal_na_and_numeric_plates(tmp_path):
    d = tmp_path / '2863_test' / 'toload_0'
    d.mkdir(parents=True)
    path = d / 'P1_DR_plate.xlsx'
    pd.DataFrame({'SMDC_ID': [5, 5, 'NA'], 'Lot': [1, 1, 'NA'],
                  'Source plate': [1001, 1002, 'NA'],
                  'Source well': ['A01', 'A02', 'A03']}).to_excel(path, index=False)
    dups = pm.find_duplicate_source_plates(str(tmp_path), ['2863'])
    assert dups.SPlateWell.iloc[0] == ['1001 A01', '1002 A02']
    pm.remove_duplicate_source_plates(dups, keep=[0], write=True)
    left = pd.read_excel(path, keep_default_na=False)
    assert left['Source plate'].tolist() == [1001, 'NA']
    assert left['Lot'].tolist() == [1, 'NA']


def test_find_and_remove_duplicate_source_plates(tmp_path):
    path = _dr_folder(tmp_path)
    dups = pm.find_duplicate_source_plates(str(tmp_path), ['2863'])
    assert dups.SMDC_ID.tolist() == [5]
    assert dups.SPlateWell.iloc[0] == ['A A01', 'B A02']

    assert pm.remove_duplicate_source_plates(dups, keep=[1]) == {str(path): (4, 3)}
    assert len(pd.read_excel(path)) == 4            # dry run left the file alone

    pm.remove_duplicate_source_plates(dups, keep=[1], write=True)
    left = pd.read_excel(path)
    assert left[left.SMDC_ID == 5]['Source plate'].tolist() == ['B']


def test_find_duplicates_empty(tmp_path):
    (tmp_path / '2863_test').mkdir()
    assert pm.find_duplicate_source_plates(str(tmp_path), ['2863']).empty


def test_fill_control_lots(tmp_path):
    path = _dr_folder(tmp_path)
    assert pm.fill_control_lots(str(tmp_path), ['2863'], {7: 3}) == [str(path)]
    assert pd.read_excel(path).loc[3, 'Lot'] != 3   # dry run
    pm.fill_control_lots(str(tmp_path), ['2863'], {7: 3}, write=True)
    assert pd.read_excel(path).loc[3, 'Lot'] == 3


def test_source_plate_barcodes():
    info = pd.DataFrame({'Source_plate': ['X-MA01; Y-Eco1', 'use lots', None, 'X-MA01', 'n.a.']})
    assert pm.source_plate_barcodes(info) == ['X-MA01', 'Y-Eco1']


def test_share_root_env_and_standard_mounts(monkeypatch, tmp_path):
    monkeypatch.setenv('smdc_share_root', str(tmp_path))
    assert pm.default_share_root() == str(tmp_path)
    monkeypatch.delenv('smdc_share_root')
    real_isdir = pm.os.path.isdir
    monkeypatch.setattr(pm.os.path, 'isdir', lambda p: p == '/mnt/z/' or (p != '/Volumes/Shared/' and real_isdir(p)))
    assert pm.default_share_root() == '/mnt/z/'


def test_wsl_and_windows_paths_normalize():
    assert pm._normalize_share_path('/mnt/z/SMDC/run.csv') == 'SMDC/run.csv'
    assert pm._normalize_share_path(r'Z:\SMDC\run.csv') == 'SMDC/run.csv'


def _block_sheet(nrows, ncols, nblocks=2):
    labels = pm._row_labels(nrows)
    rows = []
    for b in range(nblocks):
        rows.append([None] + list(range(1, ncols + 1)))           # column-number header
        rows += [[r] + [f'{b}{r}{c}' for c in range(1, ncols + 1)] for r in labels]
        rows.append([None] * (ncols + 1))                         # blank between blocks
    return pd.DataFrame(rows)


@pytest.mark.parametrize('size, nrows, ncols', [(96, 8, 12), (384, 16, 24), (1536, 32, 48)])
def test_read_platemap_block_infers_plate_size(size, nrows, ncols):
    sheet = _block_sheet(nrows, ncols)
    assert pm._infer_platemap_size(sheet) == size
    got = pm.read_platemap_block(plate_map=sheet, block=1, value_name='v')
    assert len(got) == size
    last = pm._row_labels(nrows)[-1]
    assert got.iloc[-1].tolist() == [last, f'{ncols:02d}', f'1{last}{ncols}']
