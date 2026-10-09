"""Tests for the dose-response plot — specifically the annotation block, which is the
only part of the figure a reader takes numbers off."""
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from capable_plots.assay import curves, dose_response


def _real_curve(direction="descending", ec50=50.0, hill=1.0):
    """11-point 1:4 series in duplicate, like a real plate row pair."""
    doses = 1000.0 / 4.0 ** np.arange(11)
    x = np.repeat(doses, 2)
    lo, hi = (12000.0, 2000.0) if direction == "descending" else (2000.0, 12000.0)
    return x, curves.four_pl(x, lo, hi, np.log10(ec50), hill)


def _fit(ec50=50.0, **kw):
    x, y = _real_curve(ec50=ec50)
    return x, y, curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0, **kw)


def _annotation(ax):
    return ax.texts[-1].get_text() if ax.texts else ""


def teardown_function():
    plt.close("all")


def test_annotation_carries_potency_efficacy_and_fit_quality():
    """The block must let a reader judge the fit, not just read the EC50 off it."""
    x, y, f = _fit(ref_top=1000.0)
    _, ax = plt.subplots()
    dose_response(ax, x, y, f)
    txt = _annotation(ax)
    assert "EC$_{50}$" in txt and "E$_{max}$" in txt
    assert "R$^2$" in txt and "n$_H$" in txt
    assert txt.count("\n") == 1          # potency above, fit quality below


def test_annotation_omits_fields_the_fit_does_not_carry():
    """No ref_top means no Emax; the block drops it rather than printing ``nan``."""
    x, y, f = _fit()
    _, ax = plt.subplots()
    dose_response(ax, x, y, f)
    txt = _annotation(ax)
    assert "nan" not in txt and "E$_{max}$" not in txt
    assert "R$^2$" in txt


def test_annotation_fields_can_be_narrowed():
    x, y, f = _fit(ref_top=1000.0)
    _, ax = plt.subplots()
    dose_response(ax, x, y, f, fields=("EC50",))
    assert _annotation(ax) == f"EC$_{{50}}$ {f.EC50_nM:.3g} nM"


def test_annotation_can_be_switched_off():
    x, y, f = _fit(ref_top=1000.0)
    _, ax = plt.subplots()
    dose_response(ax, x, y, f, annotate=False)
    assert not ax.texts


def test_annotation_moves_off_a_curve_that_turns_early():
    """Which corners a dose-response curve occupies depends on where its EC50 sits in
    the tested range, so a fixed anchor prints the block on top of the data whenever
    the curve turns early. Both fits here are supported — only the geometry differs."""
    anchors = {}
    for ec50 in (0.05, 50.0):            # transition at the far left / mid-range
        x, y, f = _fit(ec50=ec50)
        assert f.status == "supported"
        _, ax = plt.subplots()
        dose_response(ax, x, y, f)
        anchors[ec50] = ax.texts[-1].get_position()
    assert anchors[0.05] == (0.97, 0.97)   # bottom-left is full, so top-right
    assert anchors[50.0] == (0.03, 0.03)   # conventional corner, and it is clear
