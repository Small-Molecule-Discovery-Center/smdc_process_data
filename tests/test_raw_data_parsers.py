import csv
import pytest
"""Regression checks for the reusable EnVision reader helpers."""

import unittest
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np
import pandas as pd

import raw_data_parsers


class EnvisionDataParsingTests(unittest.TestCase):
    def test_accepts_plate_and_background_info_before_measurement(self):
        raw = pd.DataFrame(
            [
                ["Plate information", None, None],
                ["Plate", "Repeat", "Barcode"],
                [1, 1, "B123"],
                ["Background information", None, None],
                ["Plate", "Label", "Result"],
                [1, "assay", 0],
                ["Results for measurement", None, None],
            ]
        )
        plate = pd.Index([0])
        background = pd.Index([3])
        result = pd.Index([6])

        raw_data_parsers._validate_envision_section_order(
            raw, "example.csv", plate, background, result
        )

    def test_accepts_background_then_plate_info_before_measurement(self):
        raw = pd.DataFrame(
            [
                ["Background information", None, None],
                ["Plate", "Label", "Result"],
                [1, "assay", 0],
                ["Plate information", None, None],
                ["Plate", "Repeat", "Barcode"],
                [1, 1, "B123"],
                ["Results for measurement", None, None],
            ]
        )

        raw_data_parsers._validate_envision_section_order(
            raw, "example.csv", pd.Index([3]), pd.Index([0]), pd.Index([6])
        )

    def test_rejects_data_before_plate_and_background_info_with_reexport_guidance(self):
        raw = pd.DataFrame(
            [
                ["Results for measurement", None, None],
                ["Plate information", None, None],
                ["Plate", "Repeat", "Barcode"],
                [1, 1, "B123"],
                ["Background information", None, None],
            ]
        )

        with self.assertRaisesRegex(ValueError, "Re-export the data") as error:
            raw_data_parsers._validate_envision_section_order(
                raw, "wrong_order.csv", pd.Index([1]), pd.Index([4]), pd.Index([0])
            )

        self.assertIn(
            "must have exactly one Plate information and one Background information",
            str(error.exception),
        )

    def test_rejects_missing_metadata_section(self):
        raw = pd.DataFrame([["Results for measurement"]])

        with self.assertRaisesRegex(ValueError, "one Plate information"):
            raw_data_parsers._validate_envision_section_order(
                raw, "missing_info.csv", pd.Index([]), pd.Index([]), pd.Index([0])
            )

    def test_file_loader_orders_names_and_skips_loaded_subtrees(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "run10.csv").touch()
            (root / "run2.csv").touch()
            loaded = root / "loaded"
            loaded.mkdir()
            (loaded / "skip.csv").touch()

            with patch.object(raw_data_parsers.pd, "read_csv", return_value=pd.DataFrame()):
                _, filenames = raw_data_parsers._envision_load_csv_files(str(root))

        self.assertEqual(filenames, ["run2.csv", "run10.csv"])

    def test_all_blank_measurement_wells_are_restored(self):
        columns = [0] + list(range(1, 13))
        rows = [["Results for measurement"] + [None] * 12]
        rows.append([None] + list(range(1, 13)))
        for row_label in "ABCDEFGH":
            values = [float(column) for column in range(1, 13)]
            if row_label == "A":
                values[1] = np.nan
            rows.append([row_label] + values)
        raw = pd.DataFrame(rows, columns=columns)
        data_indexes = pd.DataFrame({"idx": [0], 0: ["measurement"]})
        plate_info = pd.DataFrame({"Plate": [1], "Repeat": [1], "Barcode": ["B123"]})

        wells = raw_data_parsers._envision_extract_all_data(
            raw, data_indexes, "example.csv", plate_info, "96"
        )

        self.assertEqual(len(wells), 96)
        missing_well = wells.loc[(wells["row"] == "A") & (wells["col"] == 2)]
        self.assertEqual(len(missing_well), 1)
        self.assertTrue(pd.isna(missing_well["measurement"].iloc[0]))


if __name__ == "__main__":
    unittest.main()


def test_tidy_steps_refuse_a_missing_folder(tmp_path):
    missing = str(tmp_path / 'nope')
    for fn in (raw_data_parsers.ensure_calc_spacing_in_dir, raw_data_parsers.fill_zprime_blocks_in_dir,
               raw_data_parsers._envision_load_csv_files):
        with pytest.raises(FileNotFoundError, match='Folder not found'):
            fn(missing)


def test_tidy_steps_say_how_many_files_they_found(tmp_path, caplog):
    (tmp_path / 'a.csv').write_text('x,y\n1,2\n')
    (tmp_path / 'b.csv').write_text('x,y\n1,2\n')
    with caplog.at_level('INFO', logger='raw_data_parsers'):
        assert raw_data_parsers.ensure_calc_spacing_in_dir(str(tmp_path)) == (2, 0, 0)
        assert raw_data_parsers.fill_zprime_blocks_in_dir(str(tmp_path)) == (0, 0)
    assert 'Found 2 CSV file(s)' in caplog.text and 'nothing changed' in caplog.text


def test_tidy_envision_files_touches_only_the_files_given(tmp_path):
    wide = tmp_path / 'toload.csv'
    wide.write_text('\n'.join([','.join(['A'] + ['1'] * 26), '', ','.join(['B'] + ['2'] * 24)]) + '\n')
    other = tmp_path / 'not_listed.csv'
    other.write_text(','.join(['A'] + ['1'] * 26) + '\n')
    rpt = tmp_path / 'run.rpt'
    rpt.write_text('x')
    out = raw_data_parsers.tidy_envision_files([str(wide), str(rpt)])
    assert out['file'].tolist() == [str(wide)] and out['rows trimmed'].tolist() == [1]
    lines = wide.read_text().split('\n')
    assert len(lines[0].split(',')) == 25 and lines[1] == ''      # trimmed, blank line kept
    assert len(other.read_text().split(',')) == 27                # not listed: untouched
    again = raw_data_parsers.tidy_envision_files(str(wide))
    assert again[['rows trimmed', 'blocks re-arranged', "Z' blocks filled"]].iloc[0].sum() == 0  # tidy


def test_tidy_envision_files_names_missing_files(tmp_path):
    with pytest.raises(FileNotFoundError, match='nope.csv'):
        raw_data_parsers.tidy_envision_files([str(tmp_path / 'nope.csv')])


def test_tidy_flags_reads_with_different_well_counts(tmp_path, monkeypatch):
    f = tmp_path / 'run_4310.csv'
    f.write_text('x\n')
    parsed = pd.DataFrame({'plate': ['p'] * 4, 'row': list('ABCD'), 'col': [1] * 4,
                           'ch1': [1, 1, 1, 9], 'ch2': [1, 1, 1, None], 'well': ['A01', 'B01', 'C01', 'D01']})
    monkeypatch.setattr(raw_data_parsers, 'parse_envision_files', lambda **k: parsed)
    out = raw_data_parsers.tidy_envision_files([str(f)])
    assert out['wells per read'].tolist() == ['4/3'] and bool(out['reads differ'].iloc[0])


def _envision_text(stray=False):
    head = ',' + ','.join(f'{c:02d}' for c in range(1, 13))
    lines = ['Plate information ,' + (',,,,9' if stray else ''), 'Plate,Repeat,Barcode,Measurement date',
             '1,1,,6/11/2026', '', 'Background information', 'Plate,Label,Result', '1,x,0', '1,y,0', '',
             'Results for ch1', head]
    lines += [','.join([r] + ['5'] * 12) for r in 'ABCDEFGH']
    lines += ['', '', 'Basic assay information ,', 'Assay ID: ,,,,4310']
    return '\n'.join(lines) + '\n'


def test_check_envision_csv_finds_pasted_values(tmp_path, monkeypatch):
    monkeypatch.setattr(raw_data_parsers, '_wells_per_read', lambda *a, **k: [96])
    good, bad = tmp_path / 'good.csv', tmp_path / 'bad.csv'
    good.write_text(_envision_text())
    bad.write_text(_envision_text(stray=True))
    assert raw_data_parsers.check_envision_csv(str(good)) == []
    problems = raw_data_parsers.check_envision_csv(str(bad))
    assert len(problems) == 1 and problems[0].startswith('line 1 (Plate information)')


def test_check_envision_csv_flags_reads_with_different_wells(tmp_path, monkeypatch):
    monkeypatch.setattr(raw_data_parsers, '_wells_per_read', lambda *a, **k: [253, 265, 253])
    f = tmp_path / 'f.csv'
    f.write_text(_envision_text())
    assert 'different numbers of wells (253/265/253)' in raw_data_parsers.check_envision_csv(str(f))[0]


def test_block_layout_puts_row_a_on_line_12():
    header = [''] + [f'{c:02d}' for c in range(1, 25)]
    block = [['Plate information'], ['Plate', 'Repeat'], ['1', '1'], [],
             ['Background information'], ['Plate', 'Label'], ['1', 'x'], ['1', 'y'], [],
             ['Results for ch1'], header] + [[r] + ['1'] * 24 for r in 'ABCDEFGHIJKLMNOP'] + [[]]
    assert raw_data_parsers._block_layout_problems(block + block) == []
    no_gap = block[:3] + block[4:]                                     # blank after Plate information missing
    assert 'row A is on line 11' in raw_data_parsers._block_layout_problems(no_gap)[0]
    calc = block[:9] + [['Calculated results: Calc 1']] + block[10:]   # 'Calculated results' title is fine too
    assert raw_data_parsers._block_layout_problems(calc) == []
    results_first = block[9:] + block[:9]
    assert 'comes before any Plate information' in raw_data_parsers._block_layout_problems(results_first)[0]


def test_block_layout_tells_blank_row_p_from_a_missing_row_in_any_block():
    header = [''] + [f'{c:02d}' for c in range(1, 25)]
    top = [['Plate information'], ['Plate', 'Repeat'], ['1', '1'], [],
           ['Background information'], ['Plate', 'Label'], ['1', 'x'], ['1', 'y'], [],
           ['Results for ch1'], header]
    full = top + [[r] + ['1'] * 24 for r in 'ABCDEFGHIJKLMNOP'] + [[]]
    blank_p = top + [[r] + ['1'] * 24 for r in 'ABCDEFGHIJKLMNO'] + [['P'] + [''] * 24] + [[]]
    assert raw_data_parsers._block_layout_problems(full + blank_p) == []     # empty row P is a real row
    missing_c = top + [[r] + ['1'] * 24 for r in 'ABDEFGHIJKLMNOP'] + [[], []]
    problems = raw_data_parsers._block_layout_problems(full + missing_c)  # second block only
    assert len(problems) == 1 and "expected row C" in problems[0] and "found 'D'" in problems[0]
    assert problems[0].startswith(f"line {len(full) + 14}:")


def _block(label='1', bg_rows=('1,x', '1,y'), title='Results for ch1'):
    header = ',' + ','.join(f'{c:02d}' for c in range(1, 25))
    return {'plate': ['Plate information', 'Plate,Repeat,Barcode', f'{label},1,'],
            'background': ['Background information', 'Plate,Label', *bg_rows],
            'results': [title, header] + [','.join([r] + [label] * 24) for r in 'ABCDEFGHIJKLMNOP']}


def _tidy_layout(*blocks):
    lines = []
    for b in blocks:
        lines += b['plate'] + [''] + b['background'] + [''] + b['results'] + ['']
    return lines + ['', 'Basic assay information ,', 'Assay ID: ,,,,4310']


def _arrange(tmp_path, lines):
    f = tmp_path / 'run.csv'
    f.write_text('\n'.join(lines) + '\n')
    result = raw_data_parsers.arrange_envision_blocks_in_file(str(f))
    return result, f.read_text().split('\n')[:-1]


def test_arrange_leaves_a_correct_export_alone(tmp_path):
    good = _tidy_layout(_block('1'), _block('2'))
    assert _arrange(tmp_path, good) == ((False, 0), good)


def test_arrange_pads_a_single_background_row(tmp_path):
    one = _block('1', bg_rows=('1,x',))
    lines = _tidy_layout(one, _block('2'))
    assert 'tidy step adds the missing blank line' in raw_data_parsers._block_layout_problems(
        list(csv.reader(lines)))[0]
    (changed, n), out = _arrange(tmp_path, lines)
    assert changed and n == 1
    assert out[:12] == one['plate'] + [''] + one['background'] + ['', ''] + one['results'][:2] + ['A,' + ','.join(['1'] * 24)]
    assert raw_data_parsers._block_layout_problems(list(csv.reader(out))) == []
    assert _arrange(tmp_path, out)[0] == (False, 0)                      # only once


def test_arrange_moves_plate_and_background_info_above_results_in_every_block(tmp_path):
    b1, b2 = _block('1', title='Calculated results: Calc 1'), _block('2')
    results_first = []
    for b in (b1, b2):
        results_first += b['results'] + [''] + b['plate'] + [''] + b['background'] + ['']
    results_first += ['', 'Basic assay information ,', 'Assay ID: ,,,,4310']
    (changed, n), out = _arrange(tmp_path, results_first)
    assert changed and n == 2
    assert out == _tidy_layout(b1, b2)
    assert raw_data_parsers._block_layout_problems(list(csv.reader(out))) == []


def test_arrange_adds_missing_blank_lines_between_sections(tmp_path):
    b = _block('1')
    squashed = b['plate'] + b['background'] + b['results'] + ['', '', 'Basic assay information ,']
    (changed, n), out = _arrange(tmp_path, squashed)
    assert changed and n == 1
    assert out[11].startswith('A,') and out[:4] == b['plate'] + ['']


def test_arrange_refuses_blocks_it_cannot_pair(tmp_path, caplog):
    b = _block('1')
    lines = b['plate'] + [''] + b['results'] + [''] + b['results'] + [''] + b['background'] + ['']
    with caplog.at_level('WARNING', logger='raw_data_parsers'):
        assert _arrange(tmp_path, lines) == ((False, 0), lines)
    assert 'second results section' in caplog.text and 're-export' in caplog.text
    extra = _tidy_layout(b)[:-3] + b['results'] + ['']               # a results section with no block
    with caplog.at_level('WARNING', logger='raw_data_parsers'):
        assert _arrange(tmp_path, extra)[0] == (False, 0)
    assert 'last block has only results' in caplog.text


@pytest.mark.parametrize('nrows, ncols, last', [(8, 12, 'H'), (16, 24, 'P'), (32, 48, 'AF')])
def test_block_layout_checks_row_letters_for_96_384_and_1536(nrows, ncols, last):
    labels = raw_data_parsers._grid_row_labels(nrows)
    assert labels[0] == 'A' and labels[-1] == last and len(set(labels)) == nrows
    header = [''] + [str(c) for c in range(1, ncols + 1)]
    top = [['Plate information'], ['Plate', 'Repeat'], ['1', '1'], [],
           ['Background information'], ['Plate', 'Label'], ['1', 'x'], ['1', 'y'], [],
           ['Results for ch1'], header]
    ok = top + [[r] + ['1'] * ncols for r in labels[:-1]] + [[labels[-1]] + [''] * ncols] + [[]]
    assert raw_data_parsers._block_layout_problems(ok + ok, nrows) == []
    shifted = top + [[]] + [[r] + ['1'] * ncols for r in labels] + [[]]   # extra blank line: grid starts a line late
    problems = raw_data_parsers._block_layout_problems(ok + shifted, nrows)
    assert len(problems) == 1 and problems[0].startswith(f"line {len(ok) + 1} (block starting here)")


def test_block_layout_flags_grid_rows_pipeline_pilot_would_drop():
    header = [''] + [f'{c:02d}' for c in range(1, 25)]
    top = [['Plate information'], ['Plate', 'Repeat'], ['1', '1'], [],
           ['Background information'], ['Plate', 'Label'], ['1', 'x'], ['1', 'y'], [],
           ['Results for ch1'], header]
    rows = [[r] + ['1.5', '-2', ''] + ['3'] * 21 for r in 'ABCDEFGHIJKLMNOP']
    assert raw_data_parsers._block_layout_problems(top + rows + [[]]) == []
    rows[4][3] = '1.2E+05'
    problems = raw_data_parsers._block_layout_problems(top + rows + [[]])
    assert len(problems) == 1 and 'row E' in problems[0] and 'shifts the rows below' in problems[0]


def test_tidy_keeps_all_48_columns_of_a_1536_grid(tmp_path):
    header = ',' + ','.join(f'{c:02d}' for c in range(1, 49))
    labels = raw_data_parsers._grid_row_labels(32)
    lines = (['Plate information', 'Plate,Repeat,Barcode', '1,1,', '', 'Background information', 'Plate,Label',
              '1,x', '1,y', '', 'Results for ch1', header] + [','.join([r] + ['2'] * 48) for r in labels]
             + ['', '', 'Basic assay information ,'])
    f = tmp_path / 'p1536.csv'
    f.write_text('\n'.join(lines) + '\n')
    out = raw_data_parsers.tidy_envision_files([str(f)])
    assert out['rows trimmed'].tolist() == [0] and out['problems'].tolist() == ['']
    assert f.read_text().split('\n')[:-1] == lines


def test_arrange_always_puts_row_a_on_line_12(tmp_path):
    b = _block('1')
    short_plate = dict(b, plate=b['plate'][:2])                    # plate section one row short
    no_bg_rows = dict(b, background=b['background'][:2])           # header only, no background rows
    two_line_title = dict(b, background=b['background'][:3],      # one background row and an extra
                          results=b['results'][:1] + ['Calc note'] + b['results'][1:])  # line above the grid
    for blk in (short_plate, no_bg_rows, two_line_title):
        (changed, _), out = _arrange(tmp_path, _tidy_layout(blk, _block('2')))
        rows = list(csv.reader(out))
        starts = [i for i, r in enumerate(rows) if r and r[0] == 'Plate information']
        assert [rows[s + 11][0] for s in starts] == ['A', 'A']


def test_arrange_warns_when_row_a_cannot_reach_line_12(tmp_path, caplog):
    b = _block('1', bg_rows=('1,x', '1,y', '1,z'))
    with caplog.at_level('WARNING', logger='raw_data_parsers'):
        _arrange(tmp_path, _tidy_layout(b))
    assert 'too many rows to put row A on line 12' in caplog.text


def test_tidy_says_to_close_a_locked_file(tmp_path, monkeypatch):
    f = tmp_path / 'locked.csv'
    f.write_text('A,1\n')
    def locked(*a, **k):
        raise PermissionError(13, 'Permission denied')
    monkeypatch.setattr(raw_data_parsers, 'trim_csv_columns_in_file', locked)
    with pytest.raises(PermissionError, match='locked.csv: it is open in another program'):
        raw_data_parsers.tidy_envision_files([str(f)])


def test_layout_check_flags_wrong_row_or_column_counts_in_any_block():
    def blk(rows='ABCDEFGHIJKLMNOP', ncols=24):
        header = [''] + [f'{c:02d}' for c in range(1, ncols + 1)]
        return ([['Plate information'], ['Plate', 'Repeat'], ['1', '1'], [], ['Background information'],
                 ['Plate', 'Label'], ['1', 'x'], ['1', 'y'], [], ['Results for ch1'], header]
                + [[r] + ['1'] * ncols for r in rows] + [[]])
    ok = blk()
    assert raw_data_parsers._block_layout_problems(ok + ok) == []
    short = raw_data_parsers._block_layout_problems(ok + blk('ABCDEFGHIJKLMNO'))     # A-O, then a blank line
    assert len(short) == 1 and 'expected row P' in short[0] and 'blank line' in short[0]
    narrow = raw_data_parsers._block_layout_problems(ok + blk(ncols=23))
    assert len(narrow) == 1 and 'should be 1 to 24' in narrow[0] and '23 columns' in narrow[0]
