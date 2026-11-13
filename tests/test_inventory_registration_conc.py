"""Tests for the concentration / volume rule in prepare_inventory_frames.

For wells that carry a compound: at least one of IW_CONC_UM / IW_FINAL_CONC_NM must be
filled (both allowed), and IW_VOL_UL may be blank only where IW_FINAL_CONC_NM is filled.
"""
import numpy as np
import pandas as pd
import pytest

import inventory_registration as ir

PLATE = {"IP_GRID_NUM": "Plate", "IP_IPT_ID": 1, "IP_ICO_ID": 1, "IP_COLLECTION": "C"}


def _src(**cols):
    base = {"Plate": ["P1"] * 4, "Well": ["A01", "A02", "A03", "A04"], "lot": [11, 12, 13, 14]}
    base.update(cols)
    return pd.DataFrame(base)


def _run(src, require=None, **wells):
    inputs = {"IW_COORDINATES": "Well", "IW_ICL_ID": "lot", **wells}
    return ir.prepare_inventory_frames(
        src, plate_inputs=dict(PLATE), well_inputs=inputs, size=96,
        require_plate_ids=False, fill_missing=False, require_concentration=require,
    )[1]


def test_conc_um_and_volume_ok():
    w = _run(_src(), IW_CONC_UM=10000, IW_VOL_UL=20)
    assert (w["IW_CONC_UM"] == 10000).all()


def test_final_conc_alone_allows_blank_volume():
    w = _run(_src(fc=[1.0, 2.0, 3.0, 4.0]), IW_FINAL_CONC_NM="fc")
    assert list(w["IW_FINAL_CONC_NM"]) == [1.0, 2.0, 3.0, 4.0]
    assert "IW_VOL_UL" not in w.columns or w["IW_VOL_UL"].isna().all()


def test_both_concentrations_may_be_filled():
    w = _run(_src(fc=[1.0] * 4), IW_CONC_UM=10, IW_FINAL_CONC_NM="fc", IW_VOL_UL=5)
    assert w["IW_CONC_UM"].notna().all() and w["IW_FINAL_CONC_NM"].notna().all()


def test_neither_concentration_raises():
    with pytest.raises(ValueError, match="Neither IW_CONC_UM nor IW_FINAL_CONC_NM"):
        _run(_src(), IW_VOL_UL=20)


def test_one_row_missing_both_concentrations_raises_but_other_row_can_use_final_only():
    src = _src(um=[10.0, np.nan, np.nan, 10.0], fc=[np.nan, 5.0, np.nan, np.nan])
    with pytest.raises(ValueError, match=r"2 wells.*\['A03'\]|1 wells.*\['A03'\]"):
        _run(src, IW_CONC_UM="um", IW_FINAL_CONC_NM="fc", IW_VOL_UL=20)


def test_blank_volume_needs_final_conc_on_that_row():
    src = _src(um=[10.0] * 4, fc=[np.nan, 5.0, np.nan, np.nan], vol=[20, np.nan, np.nan, 20])
    with pytest.raises(ValueError, match=r"IW_VOL_UL is blank for 1 wells.*\['A03'\]"):
        _run(src, IW_CONC_UM="um", IW_FINAL_CONC_NM="fc", IW_VOL_UL="vol")
    src.loc[2, "vol"] = 20  # fix A03; A02 is blank but has IW_FINAL_CONC_NM
    w = _run(src, IW_CONC_UM="um", IW_FINAL_CONC_NM="fc", IW_VOL_UL="vol")
    assert pd.isna(w.loc[1, "IW_VOL_UL"])


def test_wells_without_compound_are_exempt_and_not_inserted():
    src = _src(um=[10.0, 10.0, np.nan, np.nan], vol=[20, 20, np.nan, np.nan])
    src.loc[[2, 3], "lot"] = np.nan  # controls / empty
    w = _run(src, IW_CONC_UM="um", IW_VOL_UL="vol")  # no conc/vol error for them...
    assert len(w) == 2                                 # ...and they aren't inserted


def test_require_flag_modes():
    src = _src()
    _run(src)  # default: nothing supplied -> rule not enforced
    _run(src, require=False, IW_VOL_UL=20)  # explicitly off
    with pytest.raises(ValueError):
        _run(src, require=True)  # explicitly on, nothing supplied


def test_typo_negative_and_text_still_rejected():
    with pytest.raises(ValueError, match="neither a source column nor a number"):
        _run(_src(), IW_CONC_UM="concc", IW_VOL_UL=20)
    with pytest.raises(ValueError, match="negative"):
        _run(_src(c=[-1.0, 1, 1, 1]), IW_CONC_UM="c", IW_VOL_UL=20)
    with pytest.raises(ValueError, match="non-numeric"):
        _run(_src(c=["x", "1", "1", "1"]), IW_CONC_UM="c", IW_VOL_UL=20)
