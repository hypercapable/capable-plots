"""The Capable house style as an executable matplotlib theme.

``house`` — pitch-deck: serif, transparent background, thick lines.

Usable three ways::

    cap.house.apply()                       # set globally
    with cap.house:                         # scoped, auto-restored
        ...
    fig, ax = plt.subplots(figsize=cap.figsize("house-slide"))

``house`` is the default, but you can derive a tweaked theme without mutating it::

    light = cap.house.customize(background="white", font="Helvetica", font_size=10)
    with light:
        ...

See :meth:`Theme.customize` for the full set of knobs (background, font, size, line
weight, palette) plus a raw-rcParams escape hatch.

A single named theme by design — simplicity first. More styles (e.g. a journal theme)
can be added later as additional ``Theme`` instances without touching callers.
Fonts degrade to DejaVu when the preferred family is unavailable, so figures never
fail to render.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import matplotlib as mpl

from .color import CAPABLE, INK, PLACEBO

# The active theme drives transparent-background defaults in ``io.save``.
_active: Theme | None = None


@dataclass
class Theme:
    name: str
    rc: dict = field(default_factory=dict)
    transparent: bool = False

    def apply(self) -> Theme:
        """Set this theme's rcParams globally and mark it active."""
        global _active
        mpl.rcParams.update(self.rc)
        _active = self
        return self

    def __enter__(self):
        self._saved = mpl.rcParams.copy()
        self._prev_active = _active
        self.apply()
        return self

    def __exit__(self, *exc):
        global _active
        mpl.rcParams.update(self._saved)
        _active = self._prev_active
        return False

    def customize(
        self,
        *,
        name: str | None = None,
        background: str | None = None,
        font: str | list[str] | None = None,
        font_size: float | None = None,
        line_width: float | None = None,
        palette=None,
        rc: dict | None = None,
    ) -> Theme:
        """Return a NEW theme derived from this one; the original is left untouched.

        All knobs are optional — omit one to inherit this theme's value:

        background
            ``"transparent"`` / ``"none"`` → transparent background; any color
            (``"white"``, ``"#fff"``, …) → opaque with that facecolor.
        font
            A family name or a fallback list (``"Helvetica"`` or
            ``["Inter", "Arial"]``); ``"DejaVu Sans"`` is appended as a safety net so
            text always renders.
        font_size
            Base font size in points.
        line_width
            Width applied to both data lines and axis spines.
        palette
            A :class:`~capable_plots.color.Palette` or list of colors for the axes
            color cycle.
        rc
            Raw rcParams overrides, applied last — an escape hatch for anything the
            named knobs above don't cover.

        Usage::

            light = cap.house.customize(background="white", font="Helvetica")
            with light:
                ...
        """
        new_rc = dict(self.rc)
        transparent = self.transparent

        if background is not None:
            if str(background).lower() in ("transparent", "none", "clear"):
                new_rc.update({"figure.facecolor": "none", "axes.facecolor": "none",
                               "savefig.facecolor": "none"})
                transparent = True
            else:
                new_rc.update({"figure.facecolor": background, "axes.facecolor": background,
                               "savefig.facecolor": background})
                transparent = False

        if font is not None:
            families = [font] if isinstance(font, str) else list(font)
            if "DejaVu Sans" not in families:
                families = families + ["DejaVu Sans"]
            new_rc["font.family"] = families

        if font_size is not None:
            new_rc["font.size"] = font_size

        if line_width is not None:
            new_rc["lines.linewidth"] = line_width
            new_rc["axes.linewidth"] = line_width

        if palette is not None:
            cycle = palette.as_list() if hasattr(palette, "as_list") else list(palette)
            new_rc["axes.prop_cycle"] = mpl.cycler(color=cycle)

        if rc:
            new_rc.update(rc)

        return Theme(name=name or f"{self.name}+custom", rc=new_rc, transparent=transparent)


def active() -> Theme | None:
    """The theme most recently applied / entered, or ``None``."""
    return _active


# ── Shared editable-export settings (Illustrator-friendly) ──────────────────────
_EDITABLE = {
    "svg.fonttype": "none",   # SVG text stays text, not paths
    "pdf.fonttype": 42,       # embed TrueType so text is editable
    "ps.fonttype": 42,
}

# ── house: pitch-deck style ─────────────────────────────────────────────────────
house = Theme(
    name="house",
    transparent=True,
    rc={
        **_EDITABLE,
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Palatino", "DejaVu Serif"],
        # Times New Roman has no superscript-2 or subscript-5/0 glyph, so EC50 and R2
        # have to be written as mathtext. STIX is the Times-metric math face; without
        # this, mathtext silently falls back to DejaVu Sans and the sub/superscripts
        # render in a different typeface from the rest of the label.
        "mathtext.fontset": "stix",
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "axes.linewidth": 1.75,
        "axes.edgecolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": "none",
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "lines.linewidth": 1.75,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "axes.prop_cycle": mpl.cycler(color=[CAPABLE, PLACEBO]),
    },
)

# ── Standard figure sizes (inches) ──────────────────────────────────────────────
FIGSIZES = {
    "house-slide": (12.0, 7.0),  # pitch-deck panel
    "square": (5.0, 5.0),
}


def figsize(name: str) -> tuple[float, float]:
    try:
        return FIGSIZES[name]
    except KeyError:
        raise KeyError(f"unknown figsize {name!r}; have {sorted(FIGSIZES)}") from None
