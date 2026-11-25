# Changelog

Newest first. Version 1.1 is the first public release; the "beta" entries below
are the development builds that preceded it and keep their original numbering
(the development build "beta 1.1" is not version 1.1).

## 1.2 (2026-11-19), the estimation stage

NEW: the estimation stage, a separate step after matching.

  * Stata command `garfte` (garfte.ado, garfte.sthlp) and Python file
    `engine/garf_estimate.py`. `garfte` runs on the matched panel `garf` left in
    memory, or with `using` on a saved garf output folder; the data in memory
    are never changed. From a terminal:
    `python garf_estimate.py --from output --method did|gap`.
  * Two methods, no default. Both start from the gap of each treated unit to
    its own weighted donors, G(i,t) = Y(i,t) - sum_j w(i,j) Y(j,t).
    `method(did)`: matching difference-in-differences, G measured from its mean
    over the matching window; needs a panel matched without replacement
    (checked: stops if a donor serves more than one treated unit); standard
    errors clustered by block. `method(gap)`: G in levels; replacement allowed;
    block-clustered standard errors plus a correction for donors shared between
    blocks.
  * ATT(e) by event time and the overall ATT (equal-weight average of the
    post-treatment ATT(e)); Student-t intervals with (treated units - 1) df and
    a warning below 30 treated units; window periods labelled (not a
    pre-trends test). Weights are renormalised over the donors observed when
    some donor-periods are missing.
  * Output: table, event-study graph (.gph/.pdf; ATT(e) as a line with a
    level()% confidence band, synthetic-control style, matching window
    shaded, overall ATT in the subtitle, legend below; colours from the current
    scheme, band and window transparent in Stata 15+), e(b)/e(V) (coefplot, test),
    and garfte_<method>_{event,units,vcov,summary}.csv in
    garf_output/estimations/; the units file holds the effect of each treated
    unit in each period.
  * Options: outcome() pre() post() level() nograph outdir() enginedir() python();
    level() defaults to 95 (not affected by set level).
  * Checks: the panel must be a complete garf matched panel (donor weights of
    every block sum to 1; otherwise "matched panel was modified").
  * Checked on the demonstration panel: point estimates equal those of the
    stacked weighted regression with block-by-period fixed effects (to 1e-14),
    standard errors equal its cluster-robust ones with the G/(G-1) factor. In
    simulations, 95% intervals cover at 93-97% (gap with heavy donor reuse:
    about 90% without the reuse correction, about 95% with it).

Changes to `garf` (the matching computations are unchanged):

  * It stores its settings on the loaded matched panel as characteristics
    `_dta[garf_*]` (T0, window, includet0, outcome, treatment, id, time,
    replacement, k, wpower, output folder), read by `garfte`.
  * FIXED: the matched panel is imported with `asdouble`. Before, decimals
    were stored as float (about 7 significant digits), which rounded large
    outcomes (e.g. GDP in euros) and the weights.
  * The output ends with a short "Next steps" note: save the matched panel,
    edit DIAGNOSTICS/00_settings.do and run the diagnostics, run `garfte`. The
    reghdfe example line after the variable list is gone (the quick DiD stays
    documented in help garf).
  * Cleaner log: engine warnings print as one "[GARF] warning: ..." line
    (no file path or echoed source line); no "(note: file ... not found)"
    when garf writes its temporary config.
  * Engine: `__version__` = 1.2 (no change to the matching).
  * The characteristics also record the whole design (garf_call, garf_design,
    garf_forest, garf_tune, garf_tvarying, garf_tinvariant, garf_if, garf_in)
    and the file the data were loaded from (garf_source; also with unsaved
    changes such as sort/xtset), plus garf_nobs, for the diagnostics.
    01_first_run.do checks that file for garf's variables and row count.

DIAGNOSTICS: run right after your own garf, nothing to edit.

  * 00_settings.do reads every setting from the garf matched panel in memory
    (or a saved one); the manual block (demo panel) is used only when there
    is none. The settings stay in force for the session.
  * 01_first_run.do takes garf's sample (if/in) and files; without store it
    reruns the design once, same seed, to recreate them.
  * 06_placebo.do backdates the design actually used, and lowers H (or
    skips) when the data do not reach back far enough for a full window.
  * 11-12: the sensitivity grid adds your design ("your design" in the
    frontier) and uses your forest options; skipped by 00_master.do on pools
    of 20,000+ units unless global GARF_GRID on.
  * 04: no more overflow with long variable names; 03-04: the outcome may
    also be among the covariates; 03-06 use the actual window (includet0).
  * 10 prints the out-of-bag fit statistics from diagnostics.json.
  * 00_master.do returns to your folder with the matched panel in memory.

New document: GARF_Reading_Guide_v1.2.pdf (source:
integrations/stata/GARF_Reading_Guide.tex), how to read every output of garf,
the diagnostics and garfte.

garfte output: notes wrapped to the log width; for method(gap), the mean
window gap is printed next to the overall ATT; no extra blank lines.

## 1.1 (2026-11-17), first public release

Same computations and results as development build beta 1.5. Changes to the
Stata command only:

  * FIXED: names and paths are written to the engine's configuration file as
    quoted strings. Variables called on, off, yes, no, true, false or null, and
    paths containing " #" or ": ", no longer make the run fail.
  * FIXED: the data are passed to the engine with `export delimited, nolabel`.
    A value-labelled treatment variable no longer makes the run fail.
    Consequence: value-labelled variables come back in the matched panel as
    their numeric codes, not as label text.
  * FIXED: the temporary copy of the data and the temporary configuration file
    are erased after the engine has run; they used to stay in the temporary
    folder.
  * Help file: the removed option minprox() is no longer listed.
  * Repository: added README.md, LICENSE (MIT), CITATION.cff, .zenodo.json.

## beta 1.5
  * CHANGED: Forest Match Quality (FMQ) is normalised by the best donor
    ATTAINABLE in the data, Smax(i) = max over controls of c(i,j), instead of by
    a "forest twin" (c = 1):  FMQ(i) = (S(i) - cbar(i)) / (Smax(i) - cbar(i)).
    1 now means every donor is as close as the closest untreated unit available.
    The old scale could not approach 1 in large samples (with tens of thousands
    of units and small leaves no pair shares a leaf in most trees), so its level
    depended on N and minleaf. The old version is still reported (fmq_twin,
    fmq_twin_mean). On the demo: k=50 gives FMQ 0.63 (0.23 on the old scale).
    New: best_cooccurrence and fmq_twin in blocks.csv; fmq_undefined counts
    blocks where no control is closer than the average one.
  * DIAGNOSTICS: the widest design of the sensitivity grid is the setting KWIDE,
    default k(1000) (= the uncapped forest kernel on pools below 1,000 untreated
    units; safe on large pools, where nokcap can give thousands of donors per
    treated unit). The balance-over-time figure wraps at four panels per row.

## beta 1.4
  * No change to the computations: results are identical to beta 1.3.
  * NEW: DIAGNOSTICS/ -- commented Stata do-files for the full design-stage
    diagnostic protocol (balance by variable and component, variance ratios,
    balance over time, pre-treatment outcome path, out-of-window placebo, ESS
    and donor reuse, overlap, FMQ by block, forest importance, sensitivity grid
    and FMQ-ESS frontier). It replaces the old "MY TESTING DOs" folder.
    See DIAGNOSTICS/README.txt.

## beta 1.3
  * Leaner encoding: the descriptor of each time-varying variable no longer
    repeats the first and last level of the window (they were exact copies of
    two lags). Each time-varying variable now gives L+7 features (L+9 before),
    so results differ slightly from beta 1.2.
  * growth is now SYMMETRIC growth, (last - first) / ((|last| + |first|) / 2):
    bounded in [-2, 2] and defined when the first value is 0 or values are
    negative (beta 1.2 used last/first - 1 and median-imputed the gaps).
  * Honest forest fit statistics from out-of-bag predictions:
    oob_balanced_accuracy and oob_auc (Stata: r(oob_bacc), r(oob_auc)).
    The in-sample AUC was removed (it is ~1 on almost any data).
  * REMOVED: min_proximity (legacy synonym of caliper; old configs that set it
    now stop with a message) and, in Stata, the unused option load.

## beta 1.2
  * weight_power (Stata: wpower()) -- donor weights are proximity^weight_power,
    normalised to sum to 1 within each treated unit. Default 1 = exactly the
    beta 1.1 weights (proportional to proximity); 0 = uniform k:1 weights;
    >1 = more weight on the closest donors (lower effective sample size).
  * FIX: without replacement, treated units are now matched hardest-to-match
    first (increasing best proximity), as intended; beta 1.1 went easiest-first.
  * FIX: the matched SMD now uses only the matched treated units (it differed
    from beta 1.1 only when a caliper left some treated units without donors).
  * NEW: Forest Match Quality (FMQ) diagnostic and blocks.csv -- see section 6.

## beta 1.1

Caliper now filters on
a MIN-MAX NORMALISED proximity, on a fixed [0,1] scale, instead of the raw
leaf-proximity kernel used in beta 1.0 (see section 3 and 5 below). A caliper
value chosen under beta 1.0 does NOT carry over -- re-pick it using the new
proximity_minmax column. The raw proximity column is unchanged.

## beta 1.0

First testing release.
