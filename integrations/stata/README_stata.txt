=======================================================================
GARF (binary) — Stata integration
=======================================================================

garf.ado is a thin Stata command that runs the packaged Python engine and
brings the results back into Stata. It keeps ALL engine features (including
optional hyperparameter tuning). Stata never reimplements the forest — it just
drives the same Python code — so nothing is lost.

How it works: garf exports the panel in memory to a temporary CSV, writes a
config from your options, shells out to the Python engine, then reports a clean
matching summary. BY DEFAULT it loads the matched panel into memory (and saves
no files). Then estimate the effect with the companion command garfte (new in
1.2; matching DiD or gap, see below), or with your own weighted DiD.


-----------------------------------------------------------------------
REQUIREMENTS
-----------------------------------------------------------------------
- Stata 14 or newer (uses export/import delimited, shell, file, regexm, egen).
- Python 3.10+ with the package dependencies. Easiest: run the package's
  setup once (setup.sh / setup.bat) so its .venv exists — garf finds and uses
  that .venv automatically, so you don't have to configure Python in Stata.
- No need for Stata's own Python integration (this uses 'shell', not python:).


-----------------------------------------------------------------------
INSTALL (once)
-----------------------------------------------------------------------
1) Set up the package environment once (from the GARF-binary folder):
       macOS/Linux:  ./setup.sh
       Windows:      setup.bat   (double-click)

2) Tell Stata where the package is, and add the /stata folder to the ado-path
   (the ado-path is NOT recursive — it must end in /stata):
       global GARF_HOME "/path/to/GARF-binary"
       adopath + "$GARF_HOME/integrations/stata"
       which garf        // should print the path to garf.ado
       which garfte      // and to garfte.ado

   To make it permanent, copy garf.ado, garf.sthlp, garfte.ado and
   garfte.sthlp into your PERSONAL
   ado folder (type  sysdir  to find it) and put the 'global GARF_HOME ...'
   line in your profile.do. Keep the GARF-binary folder itself in a fixed
   place: garf finds the engine through GARF_HOME.

   CHECK THE VERSION: Stata uses a garf.ado in the CURRENT folder before the
   ado-path, so an old copy left in a working folder silently wins.
       which garf        // first line must read: version 1.2
       which garfte      // first line must read: version 1.2
       discard           // after replacing the files, drop the old program

3) Check the help:
       help garf
       help garfte


-----------------------------------------------------------------------
NEW IN 1.2: THE ESTIMATION STAGE, garfte
-----------------------------------------------------------------------
garfte runs AFTER garf, on the matched panel garf left in memory (or, with
using, on a saved garf output folder). It never changes the data in memory.

    garfte, method(did)    matching difference-in-differences: each treated
                           unit's gap to its own weighted donors, measured from
                           its mean over the matching window. Parallel trends
                           conditional on the match. NEEDS a panel matched
                           WITHOUT replacement (garf ..., noreplacement k(#)):
                           garfte checks it and stops otherwise. SEs clustered
                           by block (treated unit).
    garfte, method(gap)    the gap itself, in levels: treated unit minus its own
                           weighted donors. Unconfoundedness given the matched
                           history. Replacement allowed; SEs clustered by block
                           plus a correction for donors shared between blocks.
    Report both: under fairly general conditions they bracket the effect.

Output: a table of ATT(e) by event time and the overall ATT (equal-weight
average of the post-treatment ATT(e)), an event-study graph, e(b)/e(V), and the
files garfte_<method>_{event,units,vcov,summary}.csv (+ .gph/.pdf) in
garf_output/estimations/. The units file has the effect of EACH treated unit in
each period, for heterogeneity analysis.
Options: outcome() pre(#) post(#) level(#) nograph outdir() enginedir() python()
  level(#): confidence level for the table, the graph's band and the files;
  default 95. The graph draws ATT(e) as a line with a confidence band
  (synthetic-control style), the matching window shaded, legend below.
  Colours follow the current scheme (set scheme); in Stata 15+ the band and
  the window are transparent (Stata 14: fixed light colours).
Window estimates are NOT a pre-trends test (see help garfte).

DIAGNOSTICS (1.2): run them right after your own garf, nothing to edit:
    garf, ... store
    do "$GARF_HOME/DIAGNOSTICS/00_master.do"
  garf stores its whole design on the matched panel; the do-files read it.
  garf records the .dta the data were loaded from; if the data were imported
  from CSV, or variables were created / rows dropped after loading, save
  that panel and first set  global GARF_DATA "/path/to/panel.dta". The sensitivity grid is skipped
  on pools of 20,000+ units unless  global GARF_GRID on.

garf itself computes exactly as in 1.1. Three small changes: it stores its
settings on the loaded panel as characteristics _dta[garf_*] (read by garfte),
imports the matched panel in double precision (asdouble), and ends its output
with a short "Next steps" note on diagnostics and estimation.


-----------------------------------------------------------------------
NEW IN BETA 1.5
-----------------------------------------------------------------------
  * CHANGED: Forest Match Quality (FMQ) is normalised by the best donor
    ATTAINABLE in the data, Smax(i) = max over controls of c(i,j), instead of by
    a "forest twin" (c = 1):  FMQ(i) = (S(i) - cbar(i)) / (Smax(i) - cbar(i)).
    1 now means every donor is as close as the closest untreated unit available.
    The old scale could not approach 1 in large samples (with tens of thousands
    of units and small leaves no pair shares a leaf in most trees), so its level
    depended on N and minleaf. The old version is still reported (fmq_twin,
    r(fmq_twin_mean)). On the demo: k=50 gives FMQ 0.63 (0.23 on the old scale).
    New: best_cooccurrence and fmq_twin in blocks.csv; r(fmq_undefined) counts
    blocks where no control is closer than the average one.
  * DIAGNOSTICS: the widest design of the sensitivity grid is the setting KWIDE,
    default k(1000) (= the uncapped forest kernel on pools below 1,000 untreated
    units; safe on large pools, where nokcap can give thousands of donors per
    treated unit). The balance-over-time figure wraps at four panels per row.


-----------------------------------------------------------------------
NEW IN BETA 1.4
-----------------------------------------------------------------------
No change to the command or the computations (results identical to 1.3).
NEW: the package folder DIAGNOSTICS/ contains commented do-files, one per
  diagnostic, for the full design-stage protocol. From a project folder:
      do "$GARF_HOME/DIAGNOSTICS/00_master.do"
  Outputs go to ./garf_diagnostics. It replaces the old "MY TESTING DOs".
  See DIAGNOSTICS/README.txt and 'help garf##diagnostics'.


-----------------------------------------------------------------------
NEW IN BETA 1.3
-----------------------------------------------------------------------
Leaner encoding: each time-varying variable (outcome included) now gives L+7
  features instead of L+9 -- the first and last level of the window were exact
  copies of two lags and have been dropped. Results differ slightly from 1.2.
Growth is now symmetric growth: (last - first) / ((|last| + |first|) / 2),
  bounded in [-2, 2] and defined for zeros and negative values.
Honest forest fit statistics from out-of-bag predictions: r(oob_bacc) (balanced
  accuracy; 0.5 = chance) and r(oob_auc). The in-sample r(auc) was removed.
Removed: minprox() (use caliper()) and the unused option load.


-----------------------------------------------------------------------
NEW IN BETA 1.2
-----------------------------------------------------------------------
wpower(#): donor weights ~ proximity^#, summing to 1 per treated unit.
  wpower(1) (default) = the beta 1.1 weights; wpower(0) = uniform k:1;
  wpower(2) or more = more weight on the closest donors (check ESS).
Engine fixes: noreplacement now matches hardest-to-match treated units first;
the matched SMD now uses only matched treated units.
The time-invariant option now abbreviates to tin() (was ambiguous with time()).
New diagnostic: Forest Match Quality (FMQ) -- see the section below. It always
rewards fewer/closer donors: read it together with ESS, never alone.


-----------------------------------------------------------------------
DEFAULT BEHAVIOUR (since v1.1)
-----------------------------------------------------------------------
By default garf LOADS the matched panel into memory and SAVES NO FILES.

Why the matched panel (and not just extra columns on your data): with
replacement, one control can be a donor to several treated units, each time
with a different weight. So the correct object is a block-expanded panel — one
block per treated unit, donors duplicated across the blocks they serve — which
is exactly what you estimate a matching DiD on. That cannot be represented as a
few new columns on your original rows.

Loading REPLACES the data in memory. Use  preserve  or save your working data
first.

To also keep the CSV files, add  store  (saves to ./garf_output) or
export("folder"). To produce files WITHOUT touching the data in memory, add
noload (with store/export).


-----------------------------------------------------------------------
USAGE
-----------------------------------------------------------------------
Default (load matched panel into memory, save nothing):
    preserve
    garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) ///
          tvarying(x1 x2) tinvariant(region)
    garfte, method(gap)                 // estimation stage (help garfte)

Matching DiD (needs a match without replacement):
    garf, ... noreplacement k(3)
    garfte, method(did)

Your own regression instead (quick DiD with common time effects):
    reghdfe Y W [pw=weight], absorb(garf_unit year) vce(cluster garf_block)

Also save the CSV files:
    garf, ... store
    garf, ... export("C:/work/run1")

Files only, don't replace the data in memory:
    garf, ... noload store

More weight on the closest donors / uniform weights:
    garf, ... wpower(2)
    garf, ... k(10) wpower(0)

With hyperparameter tuning (slow):
    garf, ... tune tuningmetric(combined)

Full diagnostic protocol (balance, placebo, ESS, overlap, FMQ, sensitivity):
    do "$GARF_HOME/DIAGNOSTICS/00_master.do"

See example.do for a complete script, and 'help garf' for every option.


-----------------------------------------------------------------------
WHAT LOADING ADDS (variables in memory after a default run)
-----------------------------------------------------------------------
    weight        design weight ~ proximity^wpower, sums to 1 per block (use as [pw=weight])
    treated_id    the treated unit each control is matched to (= block_id)
    treated_grp   1 = treated unit, 0 = selected control
    control_sel   1 = this row is a selected control (donor)
    garf_unit     numeric id of unique_id  -> use in absorb()
    garf_block    numeric id of block_id   -> use in vce(cluster )
    post          1 if time >= T0
    cooccurrence  share of trees in which donor and treated unit share a leaf
    fmq           forest match quality of the block (0 = random, 1 = best attainable donors)
    (plus role, unique_id, donor_id, proximity, leaf_distance, rank,
     and your id/year/treatment/outcome/covariate columns)

With the standard treatment coding (1 only in the treated unit's treated
years) W is already the treated-times-post indicator, so its coefficient with
unit + time fixed effects, proximity-weighted, is the ATT.


-----------------------------------------------------------------------
STORED RESULTS r()
-----------------------------------------------------------------------
    r(n_treated) r(n_controls) r(n_treated_matched) r(n_donors)
    r(mean_smd_raw) r(mean_smd_matched) r(max_smd_raw) r(max_smd_matched)
    r(mean_ess) r(fmq_mean) r(fmq_p10) r(fmq_mean_all) r(fmq_twin_mean) r(fmq_undefined)
    r(oob_bacc) r(oob_auc) r(oob_score) r(wpower) r(cv_score, if tuned)
    r(outdir) r(matched_panel)     (only when files are saved)

garfte (e-class): e(b) e(V) e(event) e(att) e(se_att) e(lb_att) e(ub_att)
    e(p_att) e(N_treated) e(N_donors) e(donors_reused) e(mean_window_gap) ...


-----------------------------------------------------------------------
OUTPUT FILES (only when store / export / outdir is used)
-----------------------------------------------------------------------
    matched_panel.csv  matches.csv  balance.csv  blocks.csv  diagnostics.json  features.csv
garfte writes to garf_output/estimations/ (or outdir()):
    garfte_<method>_event.csv  _units.csv  _vcov.csv  _summary.csv  _messages.txt
    garfte_<method>.gph / .pdf


-----------------------------------------------------------------------
FOREST MATCH QUALITY (FMQ; normalisation changed in beta 1.5)
-----------------------------------------------------------------------
For treated unit i and control j, c(i,j) = share of trees in which they land in
the same leaf (0-1). For each treated unit:
    S(i)    = sum over i's donors of weight(i,j) * c(i,j)
    cbar(i) = mean of c(i,j) over ALL controls      (= random donors)
    Smax(i) = max  of c(i,j) over ALL controls      (= best attainable donor)
    FMQ(i)  = (S(i) - cbar(i)) / (Smax(i) - cbar(i))
0 = donors no closer than random controls; 1 = every donor as close as the
closest untreated unit in the data. Undefined (missing, counted in
r(fmq_undefined)) if Smax(i) = cbar(i). The beta 1.2-1.4 version,
(S(i) - cbar(i)) / (1 - cbar(i)), is still reported as r(fmq_twin_mean) and fmq_twin: its scale
depends on N and the forest settings. Reported as r(fmq_mean) / r(fmq_p10) (matched treated units) and
r(fmq_mean_all) (all treated, unmatched = 0).
  * FMQ STRONGLY REWARDS FEWER, CLOSER DONORS: a smaller k, a tighter caliper or
    a larger weight_power systematically raises it, and k = 1 typically
    gives the highest value. It measures
    the closeness (bias) side only. Never pick a design on FMQ alone: read it
    relative to the number of donors and together with the effective sample
    size (ESS, the variance side). Good designs give up little FMQ for a large
    gain in ESS.
  * With a caliper, compare designs on fmq_mean_all (unmatched treated units
    count as 0), not fmq_mean, which rises when hard units are dropped.
  * Relative, not absolute: no universal threshold. Designs that differ only in
    k / caliper / weight_power / replacement share the same forest and are fully
    comparable; across covariate sets or forest settings only roughly.
  * It sees only what the forest sees (dimensions that predict treatment): it
    does not replace the balance table or the pre-treatment outcome path.



-----------------------------------------------------------------------
TROUBLESHOOTING
-----------------------------------------------------------------------
"Set the package location first"   -> global GARF_HOME "/path/to/GARF-binary"
"engine not found"                 -> GARF_HOME/enginedir() is wrong.
help garf / garf not found         -> the ado-path must end in /stata (not
                                      .../integrations); check with  which garf.
"engine did not produce ..."       -> the Python side failed. Run the same call
                                      from a terminal to see the error, or pass
                                      python("/full/path/to/python") to force an
                                      interpreter that has pandas/scikit-learn.
lost my data after running garf    -> the default loads the matched panel and
                                      replaces memory; preserve/save first, or
                                      use noload to keep your data untouched.
progress bar looks messy in log    -> it's off by default here; 'verbose' turns
                                      the engine's own progress back on.
garfte: "no GARF matched panel in memory"
                                   -> run garf (1.2+) first, or use
                                      garfte using "garf_output", method(...)
garfte: "matched panel was modified"
                                   -> donors were dropped from the matched panel;
                                      rerun garf (dropping whole blocks or whole
                                      periods is fine).
garfte: "method(did) needs a panel matched WITHOUT replacement"
                                   -> garf ..., noreplacement k(#) with a k the
                                      donor pool can supply; or method(gap).

Notes:
- The treatment column must be 1 for a treated unit in its treated years, 0
  otherwise; a unit counts as treated if it is ever 1.
- The wrapper passes paths with forward slashes (fine on Windows too) and
  converts them for the Windows shell call automatically.
