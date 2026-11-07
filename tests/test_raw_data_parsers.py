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
