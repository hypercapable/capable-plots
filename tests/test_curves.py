"""Tests for the shared 4PL fitter — the logic that had silently diverged across
three copies. Both directions must recover known parameters and get the Emax sign
right."""
import numpy as np
import pytest

from capable_plots.assay import curves


def _doses(top_nM=10_000.0, dilution=4.0, n=8):
    return top_nM / dilution ** np.arange(n)


@pytest.mark.parametrize("direction,bottom,top", [
    ("descending", 1.0, 0.1),   # signal falls with dose (e.g. IP-One ratio)
    ("ascending", 0.1, 1.0),    # signal rises with dose (e.g. β-arrestin RLU)
])
def test_recovers_known_params(direction, bottom, top):
    x = _doses()
    true_logEC50 = np.log10(50.0)   # 50 nM
    true_hill = 1.0
    y = curves.four_pl(x, bottom, top, true_logEC50, true_hill)

    fit = curves.fit_4pl(x, y, direction=direction, ns_mean=bottom)

    assert fit.EC50_nM == pytest.approx(50.0, rel=0.05)
    assert fit.Hill == pytest.approx(true_hill, abs=0.15)
    assert fit.top == pytest.approx(top, abs=0.05)
    assert fit.R2 > 0.99
    assert not fit.flat_flag


def test_emax_sign_follows_direction():
    x = _doses()
    ns = 1.0
    # Descending: high-dose asymptote well below ns → positive activation.
    y_desc = curves.four_pl(x, ns, 0.2, np.log10(50.0), 1.0)
    f_desc = curves.fit_4pl(x, y_desc, direction="descending", ns_mean=ns, ref_top=0.2)
    assert f_desc.Emax_pct == pytest.approx(100.0, abs=5)

    # Ascending: high-dose asymptote well above ns → positive activation.
    y_asc = curves.four_pl(x, ns, 5.0, np.log10(50.0), 1.0)
    f_asc = curves.fit_4pl(x, y_asc, direction="ascending", ns_mean=ns, ref_top=5.0)
    assert f_asc.Emax_pct == pytest.approx(100.0, abs=5)


def test_flat_curve_flagged():
    x = _doses()
    y = np.full_like(x, 0.5) + np.random.default_rng(0).normal(0, 1e-4, x.size)
    fit = curves.fit_4pl(x, y, direction="descending", ns_mean=0.5)
    assert fit.flat_flag


def test_hill_bounds_respected():
    x = _doses()
    y = curves.four_pl(x, 1.0, 0.1, np.log10(50.0), 5.0)  # true hill above cap
    fit = curves.fit_4pl(x, y, direction="descending", ns_mean=1.0, hill_bounds=(0.6, 2.5))
    assert 0.6 <= fit.Hill <= 2.5


def test_bound_overrides_are_respected():
    x = _doses()
    y = curves.four_pl(x, 1.0, 0.1, np.log10(50.0), 1.0)
    # Force EC50 into a narrow window well away from the true 50 nM and confirm the
    # fitter honors the override (the fit is worse, but the constraint holds).
    fit = curves.fit_4pl(
        x, y, direction="descending", ns_mean=1.0,
        logEC50_bounds=(np.log10(200.0), np.log10(400.0)),
    )
    assert 200.0 <= fit.EC50_nM <= 400.0


def test_defaults_unchanged_by_override_support():
    # Same synthetic curve as the recovery test still recovers 50 nM with no overrides.
    x = _doses()
    y = curves.four_pl(x, 1.0, 0.1, np.log10(50.0), 1.0)
    fit = curves.fit_4pl(x, y, direction="descending", ns_mean=1.0)
    assert fit.EC50_nM == pytest.approx(50.0, rel=0.05)


def test_bad_direction_rejected():
    x = _doses()
    y = curves.four_pl(x, 1.0, 0.1, np.log10(50.0), 1.0)
    with pytest.raises(ValueError):
        curves.fit_4pl(x, y, direction="sideways")


def test_too_few_points_rejected():
    with pytest.raises(ValueError):
        curves.fit_4pl([1.0, 2.0, 3.0], [1.0, 0.5, 0.1], direction="descending")


# ── standard v1.1 behaviours ──────────────────────────────────────────────────

def _real_curve(direction="descending", ec50=50.0, hill=1.0, noise=0.0, seed=0):
    """11-point 1:4 series in duplicate, like a real plate row pair."""
    doses = 1000.0 / 4.0 ** np.arange(11)
    x = np.repeat(doses, 2)
    lo, hi = (12000.0, 2000.0) if direction == "descending" else (2000.0, 12000.0)
    y = curves.four_pl(x, lo, hi, np.log10(ec50), hill)
    if noise:
        y = y + np.random.default_rng(seed).normal(0, noise, x.size)
    return x, y


def test_ci_uses_t_not_z():
    """t(0.975, df) must be used, so t_critical exceeds 1.96 at realistic n."""
    x, y = _real_curve(noise=150.0)
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    assert f.n_obs == 22
    assert f.residual_df == 18
    assert f.t_critical == pytest.approx(2.1009, abs=1e-3)
    assert f.t_critical > 1.96


def test_ci_is_not_clamped_to_bounds():
    """An unidentifiable EC50 must report an unbounded CI, not a confident-looking one."""
    x = np.repeat(1000.0 / 4.0 ** np.arange(11), 2)
    y = 10000.0 + np.random.default_rng(1).normal(0, 2000.0, x.size)   # pure noise
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=10000.0)
    assert f.CI_width_fold > 1e3 or not np.isfinite(f.CI_width_fold)
    assert f.status != "supported"


def test_flat_curve_returns_not_raises():
    """A dead compound is a result: it must come back with a status, never an exception."""
    x = np.repeat(1000.0 / 4.0 ** np.arange(11), 2)
    y = np.full(x.size, 11500.0) + np.random.default_rng(2).normal(0, 40.0, x.size)
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=11500.0)
    assert f.status == "NE"
    assert "weak observed response" in f.gates_failed
    assert f.flat_flag


def test_clean_curve_is_supported():
    x, y = _real_curve(noise=120.0)
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    assert f.status == "supported", f.gates_failed
    assert f.gates_failed == []
    assert f.EC50_nM == pytest.approx(50.0, rel=0.1)
    assert f.CI_width_fold < 2.0


def test_median_scaling_makes_fit_units_invariant():
    """The same curve in HTRF-ratio units and in RLU units must give the same EC50."""
    x, y = _real_curve(noise=120.0)
    a = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    b = curves.fit_4pl(x, y * 1000.0, direction="descending", ns_mean=12000.0 * 1000.0)
    assert b.EC50_nM == pytest.approx(a.EC50_nM, rel=1e-6)
    assert b.Hill == pytest.approx(a.Hill, rel=1e-6)


def test_asymptotes_named_by_dose_position_not_height():
    """On a descending assay the low-dose asymptote is the LARGER value."""
    x, y = _real_curve("descending")
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    assert f.asym_lo_dose > f.asym_hi_dose
    assert f.bottom == f.asym_lo_dose          # legacy alias, reads backwards
    x, y = _real_curve("ascending")
    g = curves.fit_4pl(x, y, direction="ascending", ns_mean=2000.0)
    assert g.asym_lo_dose < g.asym_hi_dose


def test_hill_at_bound_is_a_gate_failure():
    x, y = _real_curve(hill=6.0)               # far above the 2.0 cap
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    assert "Hill at bound" in f.gates_failed
    assert f.status != "supported"


def test_offscale_ec50_is_gated_not_reported():
    x, y = _real_curve(ec50=50000.0)           # well above the 1000 nM top dose
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0)
    assert f.status != "supported"
    assert any("tested range" in r for r in f.gates_failed)


def test_reference_ok_false_adds_gate_9():
    x, y = _real_curve(noise=120.0)
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=12000.0, reference_ok=False)
    assert "same-plate reference not fitted" in f.gates_failed
    assert f.status == "provisional"           # the curve itself is fine


def test_wrong_direction_response_fails_first_gate():
    """A curve that moves the wrong way scores a negative drop and is NE."""
    x, y = _real_curve("ascending")
    f = curves.fit_4pl(x, y, direction="descending", ns_mean=2000.0)
    assert f.observed_drop_pct < 0
    assert f.status == "NE"


def test_default_hill_bounds_are_the_standard():
    assert curves.DEFAULT_HILL_BOUNDS == (0.5, 2.0)
    assert curves.DEFAULT_GATES["min_r2"] == 0.80
    assert curves.DEFAULT_GATES["max_ci_fold"] == 10.0
