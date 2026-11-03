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

**Splitting salts and solvates off registration SMILES:** use `extract_salts_solvates.py`
before compound registration. It reads the salt/solvate reference tables
(`INV_LOT_SALT`, `INV_LOT_SOLVATE`) from the database when you pass a SQLAlchemy
`engine=` (SELECT only, in a read-only transaction), or from CSV exports of those
tables, which also act as a fallback if the database read fails. Run the tests with
`pytest database_interactions/test_extract_salts_solvates.py` (needs rdkit).

#### preprocess_data
Before uploading data to the HiTS db, reorganize or prepare for instrument use

#### process_data
Analyses to replace pipeline pilot functionality

#### postprocess_data
After exporting data from HiTs, further analysis for hit selection and assay comparisons

#### data_viz
Functions for drawing plate maps, DR curves, etc.

#### other_analyses
Random useful requested analyses

#### SMDC_Institution_Info
Data visualizations on SMDC activities over time

#### data
example data files for PAMPA and Bruker SPR
