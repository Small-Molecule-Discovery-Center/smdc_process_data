"""Make the repo's script folders importable from any notebook.

    import smdc_paths          # once, at the top of a notebook
    import plate_maps as pm    # database_interactions/
    import plot_curves as pc   # data_viz/

Legacy shim for old notebooks that use flat names. database_interactions/ and
data_viz/ are installed packages (see pyproject.toml), so new code should use
`from database_interactions import plate_maps` and `from data_viz import
plot_curves` instead. Importing this module adds those two folders, found
relative to this file, to sys.path so the flat `import plate_maps` /
`import plot_curves` still work. Nothing is hard-coded to one machine.
"""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
FOLDERS = [REPO / "database_interactions", REPO / "data_viz"]

for _folder in FOLDERS:
    if str(_folder) not in sys.path:
        sys.path.insert(0, str(_folder))
