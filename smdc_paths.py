"""Make the repo's script folders importable from any notebook.

    import smdc_paths          # once, at the top of a notebook
    import plate_maps as pm    # database_interactions/
    import plot_curves as pc   # data_viz/

The modules in database_interactions/ and data_viz/ are plain scripts, not a
package, so `pip install -e .` cannot install them under their own names.
This module is installed (see pyproject.toml) and adds those folders, found
relative to this file, to sys.path. Nothing is hard-coded to one machine.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
FOLDERS = [REPO / "database_interactions", REPO / "data_viz"]

for _folder in FOLDERS:
    if str(_folder) not in sys.path:
        sys.path.insert(0, str(_folder))
