"""capable-plots: Capable Labs house style for figures + shared assay curve math.

Quick start (matplotlib)::

    import capable_plots as cap

    with cap.house:
        fig, ax = plt.subplots(figsize=cap.figsize("house-slide"))
        ax.plot(...)                                   # your normal matplotlib
        cap.style_axis(ax)
    cap.save(fig, "figure1")                           # 300dpi png + editable svg

Quick start (plotly)::

    import plotly.express as px
    import capable_plots as cap

    fig = px.scatter(df, x="dose", y="response", template=cap.plotly_house)
    fig.update_layout(**cap.plotly_figsize("house-slide"))
    cap.plotly_save(fig, "figure1")                    # 300dpi png + svg via kaleido
"""
from __future__ import annotations

from . import assay
from .axes import legend_layout, smart_legend, style_axis
from .color import (
    CAPABLE,
    CAPABLE_PAIR,
    CAPABLE_SEQ,
    COLORBLIND,
    INK,
    PLACEBO,
    VIBRANT,
    VIBRANT_DARK,
    VIBRANT_LIGHT,
    Gradient,
    Palette,
    colors,
)
from .io import save
from .plotly import (
    plotly_figsize,
    plotly_house,
    plotly_house_ctx,
    plotly_save,
)
from .style import FIGSIZES, Theme, active, figsize, house

__version__ = "0.2.0"

__all__ = [
    "CAPABLE",
    "CAPABLE_PAIR",
    "CAPABLE_SEQ",
    "COLORBLIND",
    "FIGSIZES",
    "INK",
    "PLACEBO",
    "VIBRANT",
    "VIBRANT_DARK",
    "VIBRANT_LIGHT",
    "Gradient",
    "Palette",
    "Theme",
    "__version__",
    "active",
    "assay",
    "colors",
    "figsize",
    "house",
    "legend_layout",
    "plotly_figsize",
    "plotly_house",
    "plotly_house_ctx",
    "plotly_save",
    "save",
    "smart_legend",
    "style_axis",
]
