# SMDC Data Processing Tools
A set of functions to preprocess various kinds of raw data for upload or further analysis

### Bruker SPR 9pt DR Plate Formatting: 
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Small-Molecule-Discovery-Center/smdc_preprocess_data/blob/main/Create_Bruker_9pt_DR_Input_Files.ipynb)

### pION PAMPA preprocessing:
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Small-Molecule-Discovery-Center/smdc_preprocess_data/blob/main/preprocess_pampa.ipynb)

### AVIDD Disulfide Tethering Match IDs:
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Small-Molecule-Discovery-Center/smdc_preprocess_data/blob/main/AVIDD_match_disulfides.ipynb)

### Xevo MassLynx Data Extraction:
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Small-Molecule-Discovery-Center/smdc_process_data/blob/main/Extract_Xevo_MassSpec_Data.ipynb)

### Plot LCMS Spectra
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Small-Molecule-Discovery-Center/smdc_process_data/blob/main/Plot_LCMS_Spectra.ipynb)

#### database_interactions
Notebooks and modules for querying HiTS, making plate maps and entering new data.
After `uv sync` or `pip install -e .` (see below), import these from any folder with
`from database_interactions import plate_maps, check_grid_num`. From inside this folder the flat
`import plate_maps` also works.

- **`prepare_hits_upload.ipynb`**: the operator notebook for a HiTS data upload. It follows
  Part C and E of the data upload SOP: build the `toload` folders, run the EnVision barcode /
  plate-order check, fix duplicate source plates, spot-check, and clean up after upload.
- **`plate_maps.py`**: builds plate-map workbooks from HiTS and the `toload` folders. It only
  reads from HiTS (production `QB3_HiTS` unless `hits_db` is set). Raw data is read from
  the Montreal share, found automatically at `/Volumes/Shared/` (Mac) or `Z:` /
  `/mnt/z/` (Windows/WSL). For another drive, pass `share_root=` or set `smdc_share_root`.
- **`check_grid_num.py`**: checks each `grid_num.txt` against the plates and barcodes in the
  EnVision export. Also runs on the command line: `python check_grid_num.py ROOT -o report.csv`.
- **`PLATE_MAPS.md`**: full reference for `plate_maps.py`. It covers the workbook format, the `Source_plate` dispatch words, how lots are resolved, and troubleshooting.
- **`safe_register.py`**: register plates and wells (including daughter plates). It connects to
  the test database `QB3_HiTS_nightly` unless you name `QB3_HiTS`, does a dry run (insert,
  report, roll back) unless you pass `save=True`, and saves plates and wells together, so a
  failure leaves nothing half-saved. The reshaping helpers it uses are in
  `inventory_registration.py`.

#### preprocess_data
Before uploading data to the HiTS db, reorganize or prepare for instrument use

#### process_data
Analyses to replace pipeline pilot functionality

#### postprocess_data
After exporting data from HiTs, further analysis for hit selection and assay comparisons

#### data_viz
Functions for drawing DR curves, etc. (plate map drawing moved to `database_interactions/plate_maps.py`)

#### other_analyses
Random useful requested analyses

#### SMDC_Institution_Info
Data visualizations on SMDC activities over time

#### data
example data files for PAMPA and Bruker SPR

### Install shared utilities

From this checkout, run `uv sync`, which builds `.venv` with this repo, everything the plate
maps, upload notebook and `plot_curves` import, and JupyterLab. Without uv, run
`python -m pip install -e .` (add `pip install jupyterlab` if you need it). The extras add the
optional pieces: `'.[upload]'` adds `rapidfuzz` for fuzzy barcode matching in the EnVision
check, and `'.[tethering]'` adds Pillow for the tethering plots.
This is the canonical source of the `smdc-data-proc` distribution, formerly
installed from data-proc. Reinstall from here in each consuming environment.
The install name `smdc-data-proc` is not an import name. Imports are `import raw_data_parsers`,
`import mysql_cnx`, `from process_data import tethering_analysis` and
`from database_interactions import plate_maps` (also `check_grid_num`, `inventory_registration`,
`safe_register`). Installing into another environment (for example a separate Notebooks venv)
works the same way: activate it, then run `uv pip install -e <path to this repo>`. Rerun that after
`pyproject.toml` changes, and restart the Jupyter kernel.

The scripts in `database_interactions/` and `data_viz/` (`plate_maps`, `check_grid_num`,
`inventory_registration`, `safe_register`, `plot_curves`, ...) are plain modules, not a
package. From a notebook anywhere, including outside this repo, add one line first:

```python
import smdc_paths          # installed by pip install -e .; finds this checkout
import plate_maps as pm
import plot_curves as pc
```

No `sys.path.append('/Users/.../data-proc')` lines are needed any more. data-proc is retired,
and everything it held now lives here.

`raw_data_parsers.extract_rpt_data(infile, samples_only=False, outdir="./")`
extracts Xevo sample metadata and writes spectra to `outdir/samples`.
The existing tethering analysis entry point delegates to this implementation.
Its existing shell/grep-based extraction behavior is preserved.

MySQL defaults to localhost and uses `mysql_user`, `mysql80_pwd`, and
`mysql80_port`. Remote hosts use `mysql_user`, `mysql_pwd`, and `mysql_port`.
The database defaults to `hits_db` or `QB3_HiTS`. Environment variables are
read when creating an engine, so importing the module needs no credentials.
Tethering analysis still uses the separate data-proc plate_maps/plot_curves
helpers for database analysis and plotting.

### Tests

All tests live in `tests/`. Run them from the repo root with `pytest`.
