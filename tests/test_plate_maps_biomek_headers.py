"""Biomek worklist headers: the dose-response template's column names reach the manifest."""
import pandas as pd

import plate_maps as pm


def test_biomek_template_concentration_header_is_read():
    # The Biomek worklist template's header is 'Concentration (uM)'; it normalises to
    # 'concentrationum', which used to be unmapped, so concentrations never reached '[compound] uM'.
    cols = pm._standardize_dids_columns(pd.DataFrame(columns=['SMDC ID', 'Concentration (uM)'])).columns
    assert list(cols) == ['SMDC_ID', '[compound] uM']
    assert '[compound] uM' in pm._standardize_dids_columns(pd.DataFrame(columns=['Concentration'])).columns


def test_make_dr_platemaps_reads_biomek_template_concentrations(monkeypatch):
    monkeypatch.setattr(pm, '_connect', lambda: None)
    monkeypatch.setattr(pm, 'get_lot_info', lambda row, dids, maps, users, rv: dids)
    # Headers exactly as in TEMPLATE_dose_response_plates.xlsx's long (Biomek) sheet.
    sheet = pd.DataFrame({'SMDC ID': [101, 101], 'Lot': [1, 1],
                          'Source Plate': ['SRC1', 'SRC1'], 'Source Well': ['F08', 'F08'],
                          'Destination Plate': ['DEST1', 'DEST1'], 'Destination Well': ['A10', 'A11'],
                          'Transfer Volume': [125.0, 62.5], 'Concentration (uM)': [32.0, 16.0]})
    row = pd.Series({'Grid_num': 'DEST1', 'Source_plate': 'use biomek', 'Raw_assay_ID': 2741,
                     'Date_of_screen': 20260101})
    out = pm.make_dr_platemaps({'DEST1': sheet}, row, users=[]).set_index('Destination well')
    assert out.loc['A10', '[compound] uM'] == 32.0 and out.loc['A11', '[compound] uM'] == 16.0
