"""plot_curves fixes from the Nov 13 Copilot review on PR #5 (synthetic data)."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

from data_viz import plot_curves as pc


def _raw(with_aeid):
    logc = np.repeat(np.linspace(-9, -4, 6), 2)
    resp = 100 / (1 + 10 ** ((-6.5 - logc) * 1.0))
    raw = pd.DataFrame({"Compound": "c1", "Target": "T", "logconc": logc, "pct": resp,
                        "RealPlate": "P1", "Plate": "P1", "Well_Type": "compound"})
    neg = pd.DataFrame({"Compound": "dmso", "Target": "T", "logconc": -9.0, "pct": [0.0, 1.0, -1.0, 0.5],
                        "RealPlate": "P1", "Plate": "P1", "Well_Type": "negativeControl"})
    raw = pd.concat([raw, neg], ignore_index=True)
    if with_aeid:
        raw["aeid"] = 7
    return raw


def test_fit_with_no_matching_rows_returns_empty_frame():
    df = pc.fit_scipy_curves(["nope"], _raw(False), "pct", targs=["T"])
    assert df.empty and {"Compound", "Target", "LogIC50", "pIC50", "hitc"} <= set(df.columns)


@pytest.mark.parametrize("with_aeid", [False, True])
@pytest.mark.parametrize("drop_logic50", [False, True])
def test_plot_works_without_aeid_and_with_pic50_only(tmp_path, with_aeid, drop_logic50):
    raw = _raw(with_aeid)
    curve = pc.fit_scipy_curves(["c1"], raw, "pct", targs=["T"], bounds=([-20, 80, -12, -5], [20, 120, -2, 5]),
                                init=np.array([0, 100, -6, 1]))
    assert curve.hitc.iloc[0] == 1
    if drop_logic50:  # tcpl_mc5_to_curve output has pIC50 but no LogIC50
        curve = curve.drop(columns="LogIC50")
    pc.plot_pct_inhib_curves("c1", raw, curve, targs=["T"], graphcol="pct", save_dir=str(tmp_path))
