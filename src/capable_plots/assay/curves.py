"""Shared 4-parameter-logistic (4PL) dose-response math.

One canonical implementation for every Capable functional assay, implementing
**Capable curve-fitting standard v1.1** (``calcium_assay/CURVE_FITTING_STANDARD.md``
in the research workspace). The only real difference between assay modalities —
whether the signal rises or falls with dose — is an explicit ``direction`` argument,
not a silently-diverged copy of the fitter.

The 4PL is written in canonical form so that, regardless of direction:

    asym_lo_dose = response as dose -> 0      asym_hi_dose = response at saturating dose

Note these are named for **where they sit on the dose axis, not by signal height**.
On a descending assay (IP-One, Gi cAMP) ``asym_lo_dose`` is the *larger* value. The
legacy aliases ``bottom``/``top`` are retained on the result object but read backwards
on a descending curve and should not be used in new output.

Key behaviours mandated by the standard:

* the response is **median-scaled** before fitting, so the asymptote bounds are
  dimensionless and the same bounds serve an HTRF ratio, an RLU, or a percentage;
* the EC50 confidence interval uses **Student's t with the residual df**, and is
  **never clamped** to the fit bounds — an unidentifiable EC50 reports as unbounded;
* a flat or noisy curve **never raises**. It returns a result carrying
  ``gates_failed`` and ``status``, because "no activity to the top dose" is a result
  and a dead compound must stay in the output table and on the figures.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t as _tdist

Direction = str  # "ascending" | "descending"

#: Hill bounds per standard v1.1 §2.1. Empirically grounded: across both Run 8
#: analyses no reportable curve had a Hill outside 0.87-1.45, while loose bounds
#: (0.1-8) let 17/80 dead curves run to the ceiling and report a potency.
DEFAULT_HILL_BOUNDS = (0.5, 2.0)

#: Acceptance gate thresholds, standard v1.1 §5.
DEFAULT_GATES = {
    "min_drop_pct": 15.0,    # observed response, 2 highest vs 3 lowest doses
    "min_r2": 0.80,
    "max_ci_fold": 10.0,
    "min_doses_each_side": 2,
    "emax_plateau_tol_pct": 15.0,
}


def four_pl(x, bottom, top, logEC50, hill):
    """Canonical 4PL. ``bottom`` = low-dose asymptote, ``top`` = high-dose asymptote.

    Argument names are the historical ones; see the module docstring on why
    ``asym_lo_dose``/``asym_hi_dose`` are the correct way to *report* them.
    """
    return bottom + (top - bottom) / (1 + 10 ** ((logEC50 - np.log10(x + 1e-12)) * hill))


@dataclass
class FitResult:
    # potency
    EC50_nM: float
    pEC50: float
    EC50_lo: float
    EC50_hi: float
    CI_width_fold: float
    Hill: float
    # asymptotes, named by dose position
    asym_lo_dose: float
    asym_hi_dose: float
    # fit quality
    R2: float
    n_obs: int
    residual_df: int
    t_critical: float
    observed_drop_pct: float
    # verdicts
    status: str                       # supported | provisional | NE | fit failed
    gates_failed: list[str]
    Emax_pct: float | None
    Emax_status: str | None           # supported | plateau unconfirmed | NE
    flat_flag: bool
    params: tuple[float, float, float, float]
    #: 1-sigma standard errors on (asym_lo_dose, asym_hi_dose, log10EC50, Hill).
    #: Asymptote errors are in raw response units; log10EC50 is in log10 nM.
    param_se: tuple[float, float, float, float] = (float("nan"),) * 4
    direction: str = ""
    fitter_version: str = field(default="capable-standard-v1.1")

    # ── legacy aliases; read backwards on a descending curve, do not use in output ──
    @property
    def bottom(self) -> float:
        return self.asym_lo_dose

    @property
    def top(self) -> float:
        return self.asym_hi_dose

    @property
    def reportable(self) -> bool:
        return self.status == "supported"


def _r2(y, y_hat) -> float:
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def _observed_drop_pct(x, y, direction) -> tuple[float, np.ndarray, np.ndarray]:
    """Fraction of the response window traversed between the low- and high-dose ends.

    Baseline = mean response at the 3 lowest doses; extreme = mean at the 2 highest.
    Signed by ``direction``, so a curve moving the wrong way scores negative and fails
    the first gate. This replaces a Spearman monotonicity test, which was tried and
    rejected: over a dose series with a long flat low-dose shoulder it is dominated by
    baseline noise and falsely rejected curves at R2 0.98-0.997 (standard v1.1 §5.1).
    """
    doses = np.unique(x)
    means = np.array([np.mean(y[x == d]) for d in doses])
    lo_end = float(np.mean(means[:3])) if means.size >= 3 else float(means[0])
    hi_end = float(np.mean(means[-2:])) if means.size >= 2 else float(means[-1])
    delta = (lo_end - hi_end) if direction == "descending" else (hi_end - lo_end)
    denom = max(abs(lo_end), abs(hi_end), 1e-12)
    return 100.0 * delta / denom, doses, means


def _failed_result(direction, n_obs, reason) -> FitResult:
    nan = float("nan")
    return FitResult(
        EC50_nM=nan, pEC50=nan, EC50_lo=nan, EC50_hi=nan, CI_width_fold=nan, Hill=nan,
        asym_lo_dose=nan, asym_hi_dose=nan, R2=nan, n_obs=n_obs, residual_df=0,
        t_critical=nan, observed_drop_pct=nan, status="fit failed",
        gates_failed=[reason], Emax_pct=None, Emax_status=None, flat_flag=True,
        params=(nan, nan, nan, nan), param_se=(nan,) * 4, direction=direction,
    )


def fit_4pl(
    x,
    y,
    *,
    direction: Direction,
    ns_mean: float | None = None,
    ref_top: float | None = None,
    hill_bounds: tuple[float, float] = DEFAULT_HILL_BOUNDS,
    bottom_bounds: tuple[float, float] | None = None,
    top_bounds: tuple[float, float] | None = None,
    logEC50_bounds: tuple[float, float] | None = None,
    n_starts: int = 6,
    flat_frac: float = 0.10,
    maxfev: int = 20000,
    gates: dict | None = None,
    reference_ok: bool = True,
) -> FitResult:
    """Bounded, multi-start 4PL fit with the standard-v1.1 acceptance ladder.

    Parameters
    ----------
    direction : "ascending" or "descending"
        Whether response rises (e.g. β-arrestin RLU) or falls (e.g. IP-One HTRF
        ratio) with dose. Sets initial guesses, the sign of the observed-response
        gate, and the sign of ``Emax_pct``. Never inferred from the data.
    ns_mean : float, optional
        Unstimulated baseline; anchors the Emax calculation and the flat test.
    ref_top : float, optional
        The reference compound's high-dose asymptote. When given, ``Emax_pct`` is
        reported relative to it; otherwise it is ``None``.
    hill_bounds : (float, float)
        Physical bounds on the Hill slope. A fit landing on a bound fails gate 6 —
        the bound is a detector, not a value.
    bottom_bounds, top_bounds, logEC50_bounds : (float, float), optional
        Override the fit bounds. ``bottom``/``top`` overrides are in **raw response
        units** and are rescaled internally. ``logEC50`` is in log10 nM. Defaults are
        a generous window around the observed data.
    n_starts : int
        Starting points on an evenly spaced log10(EC50) grid across the tested range.
    flat_frac : float
        Amplitude below this fraction of the baseline level flags the curve flat.
    gates : dict, optional
        Override acceptance thresholds; see ``DEFAULT_GATES``.
    reference_ok : bool
        Whether this curve's plate reference itself passed. False adds gate 9, since
        a fold against a failed anchor is meaningless.

    Returns
    -------
    FitResult
        Always. Convergence failure returns ``status="fit failed"`` rather than
        raising, so the compound stays in the output. ``ValueError`` is still raised
        for programmer errors (bad ``direction``, too few points).
    """
    if direction not in ("ascending", "descending"):
        raise ValueError(f"direction must be 'ascending' or 'descending', got {direction!r}")
    g = {**DEFAULT_GATES, **(gates or {})}

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[ok], y[ok]
    if x.size < 4:
        raise ValueError("need at least 4 finite positive-dose points to fit a 4PL")

    drop_pct, doses, _means = _observed_drop_pct(x, y, direction)

    # Median-scale so the asymptote bounds are dimensionless and portable across
    # modalities. Hard-coded bounds in raw signal units are the main reason the
    # legacy per-assay copies diverged.
    scale = float(np.median(y))
    scale = scale if abs(scale) > 1e-12 else 1.0
    ys = y / scale

    ymin, ymax = float(ys.min()), float(ys.max())
    span = ymax - ymin
    logx = np.log10(x)
    lo_log, hi_log = logx.min() - 1.0, logx.max() + 1.0
    min_hill, max_hill = hill_bounds

    pad = 0.5 * span if span > 0 else 1.0
    b_lo, b_hi = ((v / scale for v in bottom_bounds) if bottom_bounds is not None
                  else (ymin - pad, ymax + pad))
    t_lo, t_hi = ((v / scale for v in top_bounds) if top_bounds is not None
                  else (ymin - pad, ymax + pad))
    b_lo, b_hi = sorted((b_lo, b_hi))
    t_lo, t_hi = sorted((t_lo, t_hi))
    e_lo, e_hi = logEC50_bounds if logEC50_bounds is not None else (lo_log, hi_log)
    lo = [b_lo, t_lo, e_lo, min_hill]
    hi = [b_hi, t_hi, e_hi, max_hill]

    if direction == "descending":
        bottom0, top0 = ymax, ymin   # high at low dose, low at high dose
    else:
        bottom0, top0 = ymin, ymax
    bottom0 = min(max(bottom0, b_lo), b_hi)
    top0 = min(max(top0, t_lo), t_hi)

    best = None
    sweep_lo, sweep_hi = max(e_lo, logx.min()), min(e_hi, logx.max())
    for logec0 in np.linspace(sweep_lo, sweep_hi, n_starts):
        p0 = [bottom0, top0, logec0, 1.0]
        try:
            popt, pcov = curve_fit(four_pl, x, ys, p0=p0, bounds=(lo, hi), maxfev=maxfev)
        except (RuntimeError, ValueError):
            continue
        sse = float(np.sum((ys - four_pl(x, *popt)) ** 2))
        if best is None or sse < best[2]:
            best = (popt, pcov, sse)

    if best is None:
        return _failed_result(direction, int(x.size), "all starting points failed")

    popt, pcov, _sse = best
    b_s, t_s, logEC50, hill = (float(v) for v in popt)
    asym_lo, asym_hi = b_s * scale, t_s * scale
    r2 = _r2(y, four_pl(x, b_s, t_s, logEC50, hill) * scale)
    ec50 = 10 ** logEC50

    # ── CI: Student's t on the residual df, never clamped to the bounds ───────────
    n_obs = int(x.size)
    df = max(n_obs - 4, 1)
    tcrit = float(_tdist.ppf(0.975, df))
    finite_cov = bool(np.isfinite(pcov).all())
    se = float(np.sqrt(max(0.0, pcov[2, 2]))) if finite_cov else np.inf
    with np.errstate(invalid="ignore"):
        se_all = np.sqrt(np.clip(np.diag(pcov), 0, None))
    # asymptote errors came out of the median-scaled fit; return them in raw units
    param_se = (float(se_all[0] * scale), float(se_all[1] * scale),
                float(se_all[2]), float(se_all[3]))
    delta = tcrit * se
    with np.errstate(over="ignore"):
        ec_lo = 10.0 ** np.clip(logEC50 - delta, -300, 300)
        ec_hi = 10.0 ** np.clip(logEC50 + delta, -300, 300)
    with np.errstate(over="ignore", divide="ignore"):
        ci_fold = float(ec_hi / ec_lo) if ec_lo > 0 and np.isfinite(ec_hi) else np.inf
    if not np.isfinite(ci_fold):
        ci_fold = np.inf

    baseline = abs(ns_mean) if (ns_mean not in (None, 0)) else abs(float(np.median(y)))
    baseline = max(baseline, 1e-12)
    flat = abs(asym_hi - asym_lo) < flat_frac * baseline

    # ── Emax, native-relative and NS-anchored (standard §6) ──────────────────────
    emax = emax_status = None
    if ref_top is not None and ns_mean is not None:
        denom = (ns_mean - ref_top) if direction == "descending" else (ref_top - ns_mean)
        if denom != 0:
            num = (ns_mean - asym_hi) if direction == "descending" else (asym_hi - ns_mean)
            emax = 100.0 * num / denom
            hi_means = [float(np.mean(y[x == d])) for d in doses[-2:]]
            obs_pct = [100.0 * ((ns_mean - m) if direction == "descending" else (m - ns_mean))
                       / denom for m in hi_means]
            within = all(abs(o - emax) <= g["emax_plateau_tol_pct"] for o in obs_pct)
            emax_status = "supported" if within else "plateau unconfirmed"

    # ── acceptance ladder (standard §5) ──────────────────────────────────────────
    failed: list[str] = []
    d_min, d_max = float(doses.min()), float(doses.max())
    if drop_pct < g["min_drop_pct"]:
        failed.append("weak observed response")
    if not (r2 >= g["min_r2"]):
        failed.append(f"R2 below {g['min_r2']}")
    if not d_min < ec50 < d_max:
        failed.append("EC50 outside tested range")
    if not d_min < ec_lo < ec_hi < d_max:
        failed.append("95% CI not contained in tested range")
    if not ci_fold <= g["max_ci_fold"]:
        failed.append(f"CI wider than {g['max_ci_fold']:g}-fold")
    tol = 1e-3
    if hill <= min_hill + tol or hill >= max_hill - tol:
        failed.append("Hill at bound")
    if (b_s <= b_lo + tol or b_s >= b_hi - tol or t_s <= t_lo + tol or t_s >= t_hi - tol
            or logEC50 <= e_lo + tol or logEC50 >= e_hi - tol):
        failed.append("asymptote or EC50 at bound")
    k = g["min_doses_each_side"]
    if int(np.sum(doses < ec50)) < k or int(np.sum(doses > ec50)) < k:
        failed.append("transition sparsely bracketed")
    if not reference_ok:
        failed.append("same-plate reference not fitted")

    if not failed:
        status = "supported"
    elif drop_pct >= g["min_drop_pct"]:
        status = "provisional"
    else:
        status = "NE"

    return FitResult(
        EC50_nM=ec50,
        pEC50=-np.log10(ec50 * 1e-9),
        EC50_lo=float(ec_lo),
        EC50_hi=float(ec_hi),
        CI_width_fold=ci_fold,
        Hill=hill,
        asym_lo_dose=asym_lo,
        asym_hi_dose=asym_hi,
        R2=r2,
        n_obs=n_obs,
        residual_df=df,
        t_critical=tcrit,
        observed_drop_pct=float(drop_pct),
        status=status,
        gates_failed=failed,
        Emax_pct=emax,
        Emax_status=emax_status,
        flat_flag=flat,
        params=(asym_lo, asym_hi, logEC50, hill),
        param_se=param_se,
        direction=direction,
    )
