"""Per-axis cleanup — the fiddly spine/tick work factored out of every script.

``style_axis`` enforces the shared house-style axis rules: drop the top and right
spines, force the remaining spines to a clean black intersection at the origin
(no gaps), and keep scientific-notation offset text sized to the tick labels
instead of letting matplotlib shrink it irregularly.

``smart_legend`` handles the other recurring nuisance: a panel with far more series
than a default one-column legend can hold. A 16-compound SAR arm in one column
climbs straight through the curves it is supposed to be labelling.
"""
from __future__ import annotations

import math

from matplotlib.axes import Axes

from .color import INK


def style_axis(ax: Axes, *, offset_fontsize: float | None = None) -> Axes:
    """Apply Capable spine/tick conventions to an axis already drawn on.

    Works under either theme — linewidths/fonts come from the active rcParams;
    this only handles the structural cleanup that rcParams cannot express.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(INK)
    ax.tick_params(colors=INK)

    # Keep the sci-notation offset text from shrinking out of proportion.
    if offset_fontsize is None:
        offset_fontsize = ax.yaxis.get_ticklabels()[0].get_fontsize() \
            if ax.yaxis.get_ticklabels() else 8
    ax.xaxis.get_offset_text().set_fontsize(offset_fontsize)
    ax.yaxis.get_offset_text().set_fontsize(offset_fontsize)
    ax.grid(False)
    return ax


# Entry-count thresholds for legend layout. Derived from real panels rather than
# taste: a 12.5 x 18 inch 3x4 facet grid fits about 8 legend rows in the clear space
# below a descending curve, and about 14 before a two-column legend starts colliding
# with the curves as well.
_LEGEND_ROWS_PER_COL = 8
_LEGEND_MAX_COLS = 3


def legend_layout(n_entries: int, *, base_fontsize: float = 8.0,
                  rows_per_col: int = _LEGEND_ROWS_PER_COL,
                  max_cols: int = _LEGEND_MAX_COLS) -> dict:
    """Return the legend kwargs appropriate to ``n_entries``.

    Separated from :func:`smart_legend` so callers that build their own legend — or
    that want to know the column count before laying out a figure — can ask without
    drawing anything.
    """
    if n_entries <= 0:
        return {"ncol": 1, "fontsize": base_fontsize}
    ncol = min(max_cols, max(1, math.ceil(n_entries / rows_per_col)))
    # Shrink only once the entry count actually forces it; shrinking a short legend
    # just makes it hard to read for no gain.
    if n_entries <= rows_per_col:
        fs = base_fontsize
    elif ncol < max_cols:
        fs = base_fontsize * 0.85
    else:
        fs = base_fontsize * 0.78
    return {
        "ncol": ncol,
        "fontsize": round(fs, 2),
        # Tighten the furniture as the legend grows: at 16 entries the default
        # spacing wastes more area than the text occupies.
        "handlelength": 1.2 if n_entries > rows_per_col else 1.6,
        "columnspacing": 1.0,
        "labelspacing": 0.22 if n_entries > rows_per_col else 0.4,
        "borderaxespad": 0.15,
    }


def smart_legend(ax: Axes, *, handles=None, labels=None, loc: str = "best",
                 base_fontsize: float = 8.0, rows_per_col: int = _LEGEND_ROWS_PER_COL,
                 max_cols: int = _LEGEND_MAX_COLS, max_entries: int | None = None,
                 outside: bool = False, frameon: bool = False, **kwargs):
    """Add a legend whose layout scales with how many series the axis carries.

    Wraps into extra columns and steps the font down once the entry count exceeds
    what one column can hold, so a 16-series panel stays readable instead of running
    its legend up through the data.

    Parameters
    ----------
    max_entries : int, optional
        Cap the number of labelled entries. The overflow is replaced by a single
        unobtrusive "+N more" line, which is much better than a legend that covers
        the plot. Use when a panel genuinely has more series than anyone will read
        off a key.
    outside : bool
        Place the legend to the right of the axes instead of inside it. The axis box
        is not resized — pair it with ``constrained_layout`` or a figure-level
        ``subplots_adjust``.
    frameon : bool
        House style is frameless; exposed because a legend over dense data
        occasionally needs the frame back.

    Any further keyword arguments are passed through to
    :meth:`matplotlib.axes.Axes.legend` and win over the computed values.
    """
    if handles is None or labels is None:
        h, lab = ax.get_legend_handles_labels()
        handles = h if handles is None else handles
        labels = lab if labels is None else labels

    if max_entries is not None and len(labels) > max_entries:
        hidden = len(labels) - max_entries
        handles, labels = list(handles[:max_entries]), list(labels[:max_entries])
        # An invisible handle gives the count a row of its own without a key marker.
        handles.append(ax.plot([], [], ls="none", marker="none")[0])
        labels.append(f"+{hidden} more")

    opts = legend_layout(len(labels), base_fontsize=base_fontsize,
                         rows_per_col=rows_per_col, max_cols=max_cols)
    opts["frameon"] = frameon
    if outside:
        opts["loc"] = "upper left"
        opts["bbox_to_anchor"] = (1.02, 1.0)
        opts["ncol"] = 1 if len(labels) <= rows_per_col * 2 else opts["ncol"]
    else:
        opts["loc"] = loc
    opts.update(kwargs)
    return ax.legend(handles, labels, **opts)
