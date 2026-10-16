# capable-plots

Capable Labs' house style for matplotlib/seaborn/plotly figures, plus the
shared dose-response curve math used across our functional assays.

It is a **styling** library: you draw with normal matplotlib/seaborn or
plotly, and `capable-plots` makes it look right and saves it correctly. It is
*not* a plotting wrapper — the only drawing helpers are the handful of figures
we make constantly (dose-response, group box + strip).

## Install

```bash
pip install -e .            # from a clone, for development
pip install -e '.[seaborn]' # if you use group_box
pip install -e '.[plotly]'  # if you draw with plotly
```

## Quick start (matplotlib)

```python
import matplotlib.pyplot as plt
import capable_plots as cap

with cap.house:
    fig, ax = plt.subplots(figsize=cap.figsize("house-slide"))
    ax.plot(x, y, color=cap.CAPABLE)
    cap.style_axis(ax)                              # spine/tick cleanup
cap.save(fig, "figure1")                            # 300dpi PNG + editable SVG
```

## Quick start (plotly)

```python
import plotly.express as px
import capable_plots as cap

fig = px.scatter(df, x="dose", y="response", template=cap.plotly_house)
fig.update_layout(**cap.plotly_figsize("house-slide"))
cap.plotly_save(fig, "figure1")                     # PNG + SVG via kaleido
```

Set as the plotly default globally or scope it:

```python
import plotly.io as pio
pio.templates.default = cap.plotly_house            # global

with cap.plotly_house_ctx():                        # scoped, auto-restored
    ...
```

## Theme

| Theme              | Draws with        | Look                                           | Source style guide  |
|--------------------|-------------------|------------------------------------------------|---------------------|
| `cap.house`        | matplotlib/seaborn| pitch-deck: serif, transparent bg, thick lines | Capable house style |
| `cap.plotly_house` | plotly            | same, as a plotly `Template`                   | Capable house style |

One theme by design — simplicity first; more can be added later as additional
`Theme` / template instances. For matplotlib use as a context manager
(`with cap.house:`), globally (`cap.house.apply()`), or just pull sizes/colors
(`cap.figsize(...)`, `cap.colors(...)`, `cap.CAPABLE`). For plotly, pass
`template=cap.plotly_house` per figure, set `pio.templates.default`, or use
`with cap.plotly_house_ctx():`.

### Customizing

`house` is the default; derive a tweaked theme with `customize()` — it returns a
**new** theme and never mutates `house`:

```python
light = cap.house.customize(
    background="white",        # "transparent"/"none", or any color ("white", "#fff")
    font="Helvetica",          # a name or fallback list; DejaVu Sans appended as backup
    font_size=10,
    line_width=1.0,
    palette=cap.colors("colorblind"),
    rc={"figure.dpi": 200},    # escape hatch: any raw rcParams
)
with light:
    ...
```

Available palettes: `"capable_pair"` (brand placebo/treatment), `"colorblind"` (Okabe-Ito, 8-color),
and `"vibrant"` — a 9-hue energetic qualitative set for pitch decks/marketing contexts where
`CAPABLE`/`PLACEBO` read as too muted, with parallel `"vibrant_light"` (background tints) and
`"vibrant_dark"` (text-safe shades) variants at the same index per hue.

Every knob is optional — omit one to inherit `house`'s value.

## Assay curve math

```python
from capable_plots.assay import curves, plots

fit = curves.fit_4pl(x, y, direction="descending",   # or "ascending"
                     ns_mean=ns, ref_top=ref_top)
plots.dose_response(ax, x, y, fit)       # log axis + EC50 · Emax / R2 · nH block

# several curves on one panel: draw them bare, then one stats line per curve
for name, (x, y, fit), c in zip(names, series, colors):
    plots.dose_response(ax, x, y, fit, label=name, color=c, annotate=False)
plots.curve_stats(ax, [(n, s[2], c) for n, s, c in zip(names, series, colors)])
```

**House convention: every dose-response figure shows EC50, Emax, R2 and nH** — the
potency and efficacy a reader takes off the curve, plus the two numbers that say
whether its shape is real. `dose_response` puts them in a free corner for one curve;
`curve_stats` stacks one line per curve when a panel carries several. Both print
through the same formatter, `fit_stats`, and a dead curve reads `NE`, never `nan`.

One canonical 4PL fitter for every modality. Whether the signal rises or falls
with dose is a single explicit `direction=` argument — not three diverging copies.
