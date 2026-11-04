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
Notebooks for querying the database, making plate maps and entering new data

**Registering plates and wells (incl. daughter plates):** use `safe_register.py`.
It connects to the test database `QB3_HiTS_nightly` unless you name `QB3_HiTS`,
does a dry run (insert, report, roll back) unless you pass `save=True`, and puts
plates and wells in together so a failure leaves nothing half-saved. The reshaping
helpers it uses live in `inventory_registration.py`. Run the tests with
`pytest database_interactions/`.

**Plate maps:** `plate_maps.py` (also in `database_interactions/`) builds plate-map
workbooks from HiTS and the `toload` folders for data upload. It only reads from
HiTS (production `QB3_HiTS` unless `hits_db` is set). Raw data is read from
`/Volumes/Shared/` unless `smdc_share_root` is set.

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

From this checkout, run `python -m pip install -e .` (or
`python -m pip install -e '.[tethering]'` for plotting dependencies).
This is the canonical source of the `smdc-data-proc` distribution, formerly
installed from data-proc. Reinstall from here in each consuming environment.
Imports remain `import raw_data_parsers`, `import mysql_cnx`, and
`from process_data import tethering_analysis`.

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
