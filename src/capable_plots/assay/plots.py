"""The recurring assay figures, pre-styled — the "million-times" plots.

These are the only *drawing* helpers in the package. Everything else styles a figure
you draw yourself. They assume the caller has passed already-cleaned data and (for
dose-response) an already-computed :class:`~capable_plots.assay.curves.FitResult`;
parsing, controls and NS/NC handling stay in the caller.
"""
from __future__ import annotations

import numpy as np
from matplotlib.axes import Axes

from ..axes import style_axis
from .curves import FitResult, four_pl

#: Default annotation contents. A dose-response plot has to show enough for a reader
#: to judge the fit as well as read the potency off it, so the house block carries
#: potency, efficacy and the two numbers that say whether the shape is real: R2 and
#: the Hill slope. Pass a shorter tuple for a crowded panel.
ANNOTATION_FIELDS = ("EC50", "Emax", "R2", "nH")


#: Fit-quality fields, printed on the second line of the corner block.
_QUALITY = frozenset({"R2", "nH"})


def _is_dead(fit: FitResult) -> bool:
    """No reportable potency: the standard called it NE, or the fit never converged."""
    return fit.status in ("NE", "fit failed") or not np.isfinite(fit.EC50_nM)


def fit_stats(fit: FitResult, fields=ANNOTATION_FIELDS) -> list[tuple[str, str]]:
    """The house fit-statistics text, as ``(field, text)`` pairs in ``fields`` order.

    The single source for how a fit's numbers are printed, shared by the
    one-curve corner block (:func:`dose_response`) and the multi-curve lines
    (:func:`curve_stats`) so the two can never drift apart.

    Fields with no value on this fit are dropped rather than printed as ``nan`` — an
    unreferenced fit has no Emax. A dead curve (status ``NE`` or ``fit failed``)
    prints ``NE`` in place of its EC50 and drops Emax and nH, which describe a
    transition that is not there; R2 stays, because it is the evidence for the call.
    """
    dead = _is_dead(fit)
    out = []
    for f in fields:
        if f == "EC50":
            out.append((f, "NE" if dead else f"EC$_{{50}}$ {fit.EC50_nM:.3g} nM"))
        elif f == "Emax" and not dead and fit.Emax_pct is not None \
                and np.isfinite(fit.Emax_pct):
            out.append((f, f"E$_{{max}}$ {fit.Emax_pct:.0f}%"))
        elif f == "R2" and np.isfinite(fit.R2):
            out.append((f, f"R$^2$ {fit.R2:.3f}"))
        elif f == "nH" and not dead and np.isfinite(fit.Hill):
            out.append((f, f"n$_H$ {fit.Hill:.2f}"))
    return out


def _annotation(fit: FitResult, fields) -> str:
    """Two-line annotation: potency and efficacy above, fit quality below."""
    stats = fit_stats(fit, fields)
    top = [t for f, t in stats if f not in _QUALITY]
    bottom = [t for f, t in stats if f in _QUALITY]
    return "\n".join(" · ".join(part) for part in (top, bottom) if part)


#: ``loc`` -> (x, y, ha, va, stacking sign) for :func:`curve_stats`.
_STATS_LOC = {
    "upper left":  (0.025, 0.975, "left", "top", -1),
    "upper right": (0.975, 0.975, "right", "top", -1),
    "lower left":  (0.025, 0.025, "left", "bottom", +1),
    "lower right": (0.975, 0.025, "right", "bottom", +1),
}


def curve_stats(
    ax: Axes,
    entries,
    *,
    fields=ANNOTATION_FIELDS,
    loc: str = "upper left",
    fontsize: float = 8.0,
    line_spacing: float = 1.75,
    sep: str = "  ·  ",
) -> list:
    """One line of fit statistics per curve, for a panel that carries several.

    :func:`dose_response` puts its block in whichever corner the curve leaves free,
    and there are only four corners, so it cannot annotate four curves on one axis.
    This stacks one line per curve instead — ``label · EC50 · Emax · R2 · nH`` —
    coloured to match its curve, so the numbers are read off the figure without a
    separate table. Draw the curves with ``dose_response(..., annotate=False)``.

    ``entries`` is an iterable of ``(label, fit)`` or ``(label, fit, color)``.

    Lines are spaced in points, not axes fractions, so the stack keeps its spacing on
    a short panel; the mathtext subscripts need roughly 1.75x the font size. Nothing
    here moves the data: leave headroom with ``ax.set_ylim`` so the stack does not
    sit on a curve. Returns the text artists, top line first.
    """
    try:
        x0, y0, ha, va, sign = _STATS_LOC[loc]
    except KeyError:
        raise ValueError(f"loc must be one of {sorted(_STATS_LOC)}, got {loc!r}") from None
    step = fontsize * line_spacing
    texts = []
    for k, entry in enumerate(entries):
        label, fit, *rest = entry
        color = rest[0] if rest else None
        line = sep.join([str(label)] + [t for _, t in fit_stats(fit, fields)])
        texts.append(ax.annotate(
            line, xy=(x0, y0), xycoords="axes fraction",
            xytext=(0, sign * k * step), textcoords="offset points",
            ha=ha, va=va, fontsize=fontsize, color=color))
    return texts


#: Candidate annotation anchors, in preference order, as
#: ``(x, y, ha, va)`` in axes fractions. Bottom-left is first because it is the
#: conventional place for it and is free on a well-bracketed curve of either
#: direction; the rest are fallbacks for when the curve runs through it.
_ANCHORS = ((0.03, 0.03, "left", "bottom"), (0.97, 0.97, "right", "top"),
            (0.03, 0.97, "left", "top"), (0.97, 0.03, "right", "bottom"))

#: Approximate footprint of the two-line block in axes fractions. Deliberately
#: generous: a corner that is merely close to the data is worth skipping.
_BLOCK_W, _BLOCK_H = 0.52, 0.22


def _place(ax, x, y, fit, n_lines: int):
    """Pick the corner the curve does not run through.

    A dose-response curve occupies two opposite corners — which two depends on its
    direction and on where the EC50 sits in the tested range — so a fixed anchor
    collides whenever the curve turns early. Scores each candidate by how many
    plotted points and fitted-curve samples fall inside the text's footprint and
    takes the first clear one, falling back to the least-crowded.

    Points are mapped through the axes' **actual** limits, not the curve's own
    range. A panel whose y-limits are fixed before drawing — a ``sharey`` grid, or
    headroom left for stats lines — puts a flat curve in one band of the axis, and
    scoring it as if it filled the axis would call every corner crowded. Set limits
    *before* calling :func:`dose_response` for the block to see them.
    """
    xs = np.log10(x[x > 0])
    curve_x = np.logspace(xs.min(), xs.max(), 60)
    px = np.concatenate([xs, np.log10(curve_x)])
    py = np.concatenate([y, four_pl(curve_x, *fit.params)])
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    if not np.isfinite(py).all() or x0 <= 0 or x1 <= x0 or y1 <= y0:
        return _ANCHORS[0]
    fx = (px - np.log10(x0)) / (np.log10(x1) - np.log10(x0))
    fy = (py - y0) / (y1 - y0)

    h = _BLOCK_H * n_lines / 2
    counts = []
    for ax_f, ay_f, ha, _ in _ANCHORS:
        lo_x, hi_x = (ax_f, ax_f + _BLOCK_W) if ha == "left" else (ax_f - _BLOCK_W, ax_f)
        lo_y, hi_y = (ay_f, ay_f + h) if ay_f < 0.5 else (ay_f - h, ay_f)
        n = int(((fx > lo_x) & (fx < hi_x) & (fy > lo_y) & (fy < hi_y)).sum())
        if n == 0:
            return (ax_f, ay_f, ha, "bottom" if ay_f < 0.5 else "top")
        counts.append(n)
    i = int(np.argmin(counts))
    a = _ANCHORS[i]
    return (a[0], a[1], a[2], "bottom" if a[1] < 0.5 else "top")


def dose_response(
    ax: Axes,
    x,
    y,
    fit: FitResult,
    *,
    label: str | None = None,
    color: str | None = None,
    annotate: bool = True,
    fields=ANNOTATION_FIELDS,
) -> Axes:
    """Plot points + fitted 4PL on a log-dose axis with the standard annotation.

    Annotation format, two lines::

        EC50 3.21 nM · Emax 87%
        R2 0.996 · nH 1.12

    Anything the fit does not carry is omitted, so a flat curve annotates cleanly
    instead of printing ``nan``. Narrow the block with ``fields=("EC50", "Emax")``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    ax.set_xscale("log")
    ax.scatter(x, y, s=28, color=color, edgecolor="black", linewidth=0.5, zorder=3,
               label=label)

    xs = np.logspace(np.log10(x[x > 0].min()), np.log10(x.max()), 200)
    ax.plot(xs, four_pl(xs, *fit.params), color=color, zorder=2)

    ax.set_xlabel("Dose (nM)")
    if annotate:
        txt = _annotation(fit, fields)
        if txt:
            tx, ty, ha, va = _place(ax, x, y, fit, txt.count("\n") + 1)
            ax.text(tx, ty, txt, transform=ax.transAxes, va=va, ha=ha, linespacing=1.3)

    style_axis(ax)
    return ax


def group_box(
    ax: Axes,
    data,
    x: str,
    y: str,
    *,
    order=None,
    palette=None,
) -> Axes:
    """Transparent box (median/IQR) + individual points, per the house-style rules.

    Requires seaborn (install the ``seaborn`` extra). Significance annotation is left
    to the caller so the package stays free of statistical-test policy.
    """
    try:
        import seaborn as sns
    except ImportError as e:  # pragma: no cover - optional dependency
        raise ImportError("group_box needs seaborn: pip install 'capable-plots[seaborn]'") from e

    # hue=x + legend=False is seaborn's forward-compatible way to color by group.
    sns.boxplot(data=data, x=x, y=y, hue=x, legend=False, order=order, ax=ax,
                palette=palette,
                boxprops={"alpha": 0.5, "edgecolor": "black", "linewidth": 1.5},
                showfliers=False)
    sns.stripplot(data=data, x=x, y=y, hue=x, legend=False, order=order, ax=ax,
                  palette=palette,
                  alpha=1.0, edgecolor="black", linewidth=0.5, size=6, jitter=0.15)
    style_axis(ax)
    return ax
