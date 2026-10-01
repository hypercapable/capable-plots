# Changelog

All notable changes to `capable-plots` are documented here. This project follows
[semantic versioning](https://semver.org/).

## [0.4.0] — unreleased

### Added
- **`smart_legend(ax, ...)` and `legend_layout(n)`** — legend layout that scales with
  how many series an axis carries. One column up to 8 entries, then wraps to 2 and 3
  columns and steps the font down, tightening handle/label spacing as it grows. A
  16-compound SAR arm in a default one-column legend runs straight up through the
  curves it is labelling; this is the fix, lifted out of a per-experiment script and
  made reusable.
  - `max_entries=N` caps the key and replaces the overflow with a single `"+N more"`
    line, for panels with more series than anyone will read off a legend.
  - `outside=True` anchors the legend to the right of the axes.
  - Frameless by default, per house style; any `legend()` kwarg passed through wins
    over the computed layout.
  - `legend_layout(n)` returns the kwargs without drawing, for callers building their
    own legend or sizing a figure before laying it out.

## [0.3.0] — unreleased

Implements **Capable curve-fitting standard v1.1**. `capable_plots.assay.fit_4pl` is
now the single canonical fitter for every Capable functional assay.

### Changed — behaviour (review before upgrading an existing analysis)
- **`fit_4pl` no longer raises on a bad curve.** Convergence failure returns a
  `FitResult` with `status="fit failed"` instead of `RuntimeError`, so a dead compound
  stays in the output table and on the figures. `ValueError` is still raised for
  programmer errors (bad `direction`, fewer than 4 points).
- **EC50 confidence intervals use Student's t on the residual df**, not z = 1.96. At a
  typical n = 22 that is t = 2.10, so intervals widen by ~7% in log units.
- **CIs are never clamped to the fit bounds.** An unidentifiable EC50 now reports an
  unbounded interval rather than a confident-looking finite one.
- **Default `hill_bounds` tightened to (0.5, 2.0)** from (0.6, 2.5). Grounded in Run 8:
  no reportable curve in either of two independent analyses had a Hill outside
  0.87–1.45, while loose bounds let 17 of 80 dead curves run to the ceiling and report
  a potency. A fit landing on a bound is now a gate failure, not a value.
- **The response is median-scaled before fitting**, so asymptote bounds are
  dimensionless and identical bounds serve an HTRF ratio, an RLU, or a percentage.
  Fitted EC50 is now invariant to the units of the input.
- Start selection switched from best R² to lowest SSE (identical ranking, cheaper).

### Added
- **A nine-gate acceptance ladder** producing `status` of `supported` /
  `provisional` / `NE` / `fit failed`, with `gates_failed` naming every failure in
  readable text. Thresholds are in `DEFAULT_GATES` and overridable per call.
- `reference_ok=False` adds the same-plate-reference gate — a fold against a failed
  anchor is meaningless.
- `asym_lo_dose` / `asym_hi_dose` on `FitResult`, named by **dose-axis position rather
  than signal height**. On a descending assay `asym_lo_dose` is the *larger* value.
  `bottom` / `top` remain as properties but read backwards on a descending curve and
  should not be used in new output.
- `CI_width_fold`, `n_obs`, `residual_df`, `t_critical`, `observed_drop_pct`,
  `Emax_status`, `direction` and `fitter_version` on `FitResult`.
- `Emax_status` distinguishes a confirmed plateau from `"plateau unconfirmed"` —
  the two highest observed dose means must sit within 15 percentage points of the
  modelled asymptote.

### Rejected after testing
- A Spearman monotonicity gate. Over a dose series with a long flat low-dose shoulder
  it is dominated by baseline noise: on Run 8 it uniquely rejected five curves at
  R² 0.980–0.997 and CI 1.39–2.31×, and caught nothing the other gates missed.
  Direction is instead enforced by the sign of `observed_drop_pct`.
- A fixed-Hill (3PL) fallback for curves failing only on `Hill at bound`. It rescued
  zero curves on Run 8; not worth the complexity.

## [0.2.0] — unreleased

### Added
- Plotly parity via `capable_plots.plotly`: registers a plotly `Template`
  named `"capable_house"`, exposed as `cap.plotly_house`. Same palette,
  Times New Roman font, transparent backgrounds, hidden top/right axis lines,
  no grid, ink-black axes as the matplotlib `house` theme.
- `cap.plotly_house_ctx()` context manager sets and restores
  `plotly.io.templates.default`, mirroring `with cap.house:` for matplotlib.
- `cap.plotly_save(fig, name)` writes PNG + SVG via the kaleido backend,
  mirroring `cap.save` for matplotlib.
- `cap.plotly_figsize(name)` returns pixel dimensions matching the matplotlib
  `figsize` inch presets (12in × 7in `house-slide`, 5in × 5in `square`).
- New optional dependency extra: `pip install 'capable-plots[plotly]'`
  (installs `plotly>=5.20` + `kaleido>=0.2`).

## [0.1.0] — unreleased

Initial release.

### Added
- Universal styling core: `house` theme (context-manager / global),
  `style_axis`, `save` (300 dpi PNG + editable SVG),
  `figsize`, and `Palette`/`Gradient` color primitives with the Capable brand
  colors and a colorblind-safe default cycle.
- `assay` domain pack: canonical `four_pl` + bounded multi-start `fit_4pl` with
  an explicit `direction=` argument, `FitResult`, and the `dose_response` /
  `group_box` figure helpers.
- `fit_4pl` accepts optional per-parameter bound overrides (`bottom_bounds`,
  `top_bounds`, `logEC50_bounds`) alongside `hill_bounds`, so a single fitter can be
  tuned per assay modality instead of forking the function. Defaults unchanged.
- `Theme.customize(...)` derives a new theme from `house` with semantic knobs
  (`background`, `font`, `font_size`, `line_width`, `palette`) plus a raw-`rc` escape
  hatch. Returns a new theme; the default `house` is never mutated.
