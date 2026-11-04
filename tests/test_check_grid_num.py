"""check_upload_folder: the notebook's EnVision plate-order step."""
import pytest

import check_grid_num as cg


def _folder(tmp_path, grid, barcodes):
    d = tmp_path / '1_test' / 'toload_0' / '20260101'
    d.mkdir(parents=True)
    (d / 'grid_num.txt').write_text('\n'.join(grid) + '\n')
    blocks = ''.join(f'Plate information\nPlate,Barcode\n{i},{b}\n\n' for i, b in enumerate(barcodes, 1))
    (d / 'export.csv').write_text(blocks)


def test_swapped_order_is_flagged(tmp_path):
    _folder(tmp_path, ['P1', 'P2'], ['P2', 'P1'])
    problems = cg.check_upload_folder(str(tmp_path))
    assert problems.verdict.tolist() == ['ORDER_MISMATCH']


def test_matching_order_is_ok(tmp_path):
    _folder(tmp_path, ['P1', 'P2'], ['P1', 'P2'])
    assert cg.check_upload_folder(str(tmp_path)).empty


def test_no_grid_files_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        cg.check_upload_folder(str(tmp_path))


def _files(tmp_path, grid, files):
    """files: {csv name: [barcodes in that export]}."""
    d = tmp_path / '1_test' / 'toload_0' / '20260101'
    d.mkdir(parents=True)
    (d / 'grid_num.txt').write_text('\n'.join(grid) + '\n')
    for name, barcodes in files.items():
        blocks = ''.join(f'Plate information\nPlate,Barcode\n{i},{b}\n\n' for i, b in enumerate(barcodes, 1))
        (d / name).write_text(blocks)


def test_single_plate_file_named_for_another_plate_is_flagged(tmp_path):
    # Grid and barcode agree, but the export is named for a different plate.
    _files(tmp_path, ['P1'], {'Q9.csv': ['P1']})
    assert cg.check_upload_folder(str(tmp_path)).verdict.tolist() == ['ORDER_MISMATCH']


def test_multi_plate_exports_in_order_are_ok(tmp_path):
    _files(tmp_path, ['P1', 'P2', 'P3', 'P4'], {'run_a.csv': ['P1', 'P2'], 'run_b.csv': ['P3', 'P4']})
    assert cg.check_upload_folder(str(tmp_path)).empty


def test_files_sorted_like_the_loader(tmp_path):
    # Lexical order would be P1, P10, P2; the loader reads P1, P2, P10.
    _files(tmp_path, ['P1', 'P2', 'P10'], {'P1.csv': ['P1'], 'P2.csv': ['P2'], 'P10.csv': ['P10']})
    assert cg.check_upload_folder(str(tmp_path)).empty


def test_inexact_barcode_match_needs_review(tmp_path):
    # RUN2026-P01 vs RUN2026-P99 share a number but are different plates.
    _folder(tmp_path, ['RUN2026-P01', 'RUN2026-P02'], ['RUN2026-P99', 'RUN2026-P02'])
    assert cg.check_upload_folder(str(tmp_path)).verdict.tolist() == ['REVIEW']
