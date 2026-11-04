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

fails = 0
for c in CASES:
    got, want = pm._normalize_share_path(c), original_normalize(c)
    status = "ok " if got == want else "DIFF"
    if got != want:
        fails += 1
    print(f"  [{status}] {c!r}\n         -> {got!r}")
print(f"_normalize_share_path matches original on {len(CASES)-fails}/{len(CASES)} unchanged cases")
assert fails == 0

# The deliberate divergence: 'Montreal/Shared/' used to slice [17:] against a
# 16-character prefix, eating the 'S' of 'SMDC' and producing a path that could
# never resolve. This asserts the FIXED behaviour; if it ever fails, the
# off-by-one has come back.
print("\nMontreal/Shared/ off-by-one (deliberate divergence from original):")
fixed = pm._normalize_share_path(FIXED_CASE)
print(f"  was (buggy): {original_normalize(FIXED_CASE)!r}")
print(f"  now (fixed): {fixed!r}")
assert fixed == 'SMDC/run11.txt', f"expected 'SMDC/run11.txt', got {fixed!r}"
assert original_normalize(FIXED_CASE) == 'MDC/run11.txt'  # documents what was wrong

# Every branch must strip exactly its prefix and no more.
for prefix in ['Montreal/Shared/', '/shared/', 'shared/', '/Volumes/Shared/',
               'smb://montreal.ucsf.edu/shared/', 'smb://montreal/Shared/',
               'mb://montreal.ucsf.edu/shared/', '/montreal.ucsf.edu/shared/']:
    got = pm._normalize_share_path(prefix + 'SMDC/x.txt')
    assert got == 'SMDC/x.txt', f"{prefix!r} -> {got!r}, expected 'SMDC/x.txt'"
print("  every prefix branch strips exactly its own length")

# No result may start with '/', or os.path.join would discard the share root.
for c in CASES + [FIXED_CASE]:
    assert not pm._normalize_share_path(c).startswith('/'), c
print("no result starts with '/' (os.path.join stays safe)")

# _assay_dir_name: the spellings that used to make the two functions disagree.
print("\n_assay_dir_name across spellings:")
for rid in [123, 123.0, '123', '123.0']:
    print(f"  {rid!r:>10} -> {pm._assay_dir_name(rid, 'MyAssay')!r}")
assert len({pm._assay_dir_name(r, 'MyAssay') for r in [123, 123.0, '123', '123.0']}) == 1
print("  all four agree")
print(f"  {'map'!r} -> {pm._assay_dir_name('map', 'MyAssay')!r}")
assert pm._assay_dir_name('map', 'MyAssay') == 'MyAssay'

# handle_zero_concentrations: drops zeros, no prompt.
print("\nhandle_zero_concentrations:")
df = pd.DataFrame({
    'Destination well': ['A01', 'A02', 'A03', 'A04'],
    '[compound] uM': [10.0, 1.0, 0.0, 0.0],
})
sys.stdin = None  # any input() call would now raise rather than hang
out = pm.handle_zero_concentrations(df)
print(f"  {len(df)} wells in, {len(out)} out, kept {out['[compound] uM'].tolist()}")
assert len(out) == 2 and 0.0 not in out['[compound] uM'].tolist()

out2 = pm.handle_zero_concentrations(pd.DataFrame({'[compound] uM': [10.0, 1.0]}))
assert len(out2) == 2
print("  no zeros -> unchanged")

# read_platemap_block now rejects a bad size instead of UnboundLocalError.
try:
    pm.read_platemap_block(plate_map=pd.DataFrame(), block=0, size=1596)
except ValueError as e:
    print(f"\nread_platemap_block(size=1596) -> ValueError: {e}")

print("\nALL SMOKE TESTS PASSED")


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
