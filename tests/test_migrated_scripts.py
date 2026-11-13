"""Code moved from data-proc import cleanly and find their files from the repo."""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_package_imports_work():
    code = ("from database_interactions import plate_maps, check_grid_num, inventory_registration; "
            "from data_viz import plot_curves; print('ok')")
    out = subprocess.run([sys.executable, '-c', code], cwd=REPO, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == 'ok'


def test_tethering_analysis_imports_from_the_repo():
    out = subprocess.run([sys.executable, '-c', 'import process_data.tethering_analysis as t; print(t.pm.__file__, t.pc.__file__)'],
                         cwd=REPO, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert str(REPO / 'database_interactions') in out.stdout and str(REPO / 'data_viz') in out.stdout
