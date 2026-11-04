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
