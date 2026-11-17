# Changelog

Newest first. Version 1.1 is the first public release; the "beta" entries below
are the development builds that preceded it and keep their original numbering
(the development build "beta 1.1" is not version 1.1).

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
