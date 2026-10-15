"""Tests for the dose-response plot — specifically the annotation block, which is the
only part of the figure a reader takes numbers off."""
import itertools

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


# ── fit_stats / curve_stats ───────────────────────────────────────────────────

from capable_plots.assay import curve_stats, fit_stats


def _dead_fit():
    """A flat row: the standard has to call it NE, and the figure has to say so."""
    x, _ = _real_curve()
    y = np.full_like(x, 12000.0) + np.random.default_rng(0).normal(0, 150, x.size)
    return curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)


def test_block_and_lines_print_identical_numbers():
    """One formatter feeds both the corner block and the stacked lines, so a figure
    that uses both can never show two versions of the same fit."""
    x, y, f = _fit(ref_top=1000.0)
    _, ax = plt.subplots()
    dose_response(ax, x, y, f)
    block = _annotation(ax)
    _, ax2 = plt.subplots()
    (line,) = curve_stats(ax2, [("cmpd", f)])
    for _, piece in fit_stats(f):
        assert piece in block and piece in line.get_text()


def test_curve_stats_one_line_per_curve_carrying_every_field():
    fits = [(name, _fit(ec50=e, ref_top=1000.0)[2])
            for name, e in (("a", 5.0), ("b", 50.0), ("c", 500.0))]
    _, ax = plt.subplots()
    texts = curve_stats(ax, fits)
    assert [t.get_text().split("  ·  ")[0] for t in texts] == ["a", "b", "c"]
    for t in texts:
        s = t.get_text()
        assert all(k in s for k in ("EC$_{50}$", "E$_{max}$", "R$^2$", "n$_H$"))
        assert "nan" not in s


def test_curve_stats_lines_take_their_curve_colour():
    f = _fit()[2]
    _, ax = plt.subplots()
    texts = curve_stats(ax, [("a", f, "#d62728"), ("b", f, "#1f77b4")])
    assert [t.get_color() for t in texts] == ["#d62728", "#1f77b4"]


def test_curve_stats_lines_do_not_overlap_when_rendered():
    """The stack is spaced in points; mathtext subscripts make each line taller than
    its point size, which is what a fraction-of-axes spacing got wrong."""
    fits = [(f"compound {k}", _fit(ec50=10.0 ** k, ref_top=1000.0)[2]) for k in range(4)]
    fig, ax = plt.subplots(figsize=(4, 3))
    texts = curve_stats(ax, fits)
    fig.canvas.draw()
    boxes = [t.get_window_extent() for t in texts]
    for a, b in itertools.pairwise(boxes):
        assert a.y0 >= b.y1, "adjacent stats lines overlap"


def test_dead_curve_says_NE_in_both_places_and_never_nan():
    """A dead compound stays on the figure (standard v1.1) — labelled, not blank."""
    f = _dead_fit()
    assert f.status in ("NE", "fit failed")
    pieces = [t for _, t in fit_stats(f)]
    assert pieces[0] == "NE"
    assert not any("E$_{max}$" in p or "n$_H$" in p for p in pieces)
    _, ax = plt.subplots()
    (line,) = curve_stats(ax, [("dead", f)])
    assert "NE" in line.get_text() and "nan" not in line.get_text()


def test_curve_stats_lower_loc_stacks_upward():
    f = _fit()[2]
    fig, ax = plt.subplots()
    texts = curve_stats(ax, [("a", f), ("b", f)], loc="lower right")
    fig.canvas.draw()
    a, b = (t.get_window_extent() for t in texts)
    assert b.y0 >= a.y1                  # second line sits above the first


def test_curve_stats_rejects_unknown_loc():
    import pytest
    with pytest.raises(ValueError, match="loc must be one of"):
        curve_stats(plt.subplots()[1], [], loc="middle")


def test_annotation_placement_respects_limits_set_before_drawing():
    """A flat curve on a shared-y grid sits in the top band of its panel. Scoring it as
    if it filled the axis called every corner crowded and dropped the block onto a
    data point; against the real limits the bottom corners are visibly empty."""
    f = _dead_fit()
    x, _ = _real_curve()
    y = np.full_like(x, 12000.0) + np.random.default_rng(0).normal(0, 150, x.size)
    _, ax = plt.subplots()
    ax.set_ylim(0, 14000)                # what sharey / headroom would impose
    dose_response(ax, x, y, f)
    _, ty = ax.texts[-1].get_position()
    assert ty < 0.5                      # bottom corner, below the flat data
