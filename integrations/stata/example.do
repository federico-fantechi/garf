*-----------------------------------------------------------------------
* GARF (binary) — Stata example
*
* Prerequisites:
*   1) Python 3.10+ installed, and the package environment set up once:
*        run  setup.sh / setup.bat  in the GARF-binary folder (builds .venv).
*   2) Edit the GARF_HOME path below to point at your GARF-binary folder.
*      Note the ado-path must end in /stata (it is not recursive).
*-----------------------------------------------------------------------
clear all
set more off

* --- point Stata at the package -----------------------------------------
global GARF_HOME "/path/to/GARF-binary"                 // <-- EDIT THIS
adopath + "$GARF_HOME/integrations/stata"               // must end in /stata
which garf                                              // should print the .ado path (version 1.2)
which garfte

* --- get some data ------------------------------------------------------
* This example uses the synthetic demo panel. Generate it once by running the
* demo (demo/run_demo.sh or run_demo.bat), which writes demo/input/panel.csv.
import delimited "$GARF_HOME/demo/input/panel.csv", clear varnames(1) case(preserve)

* Your own data instead: load your long panel (one row per unit-year) with a
* binary treatment column (1 in the treated unit's treated years, else 0).

* --- run GARF -----------------------------------------------------------
* DEFAULT behaviour: the matched panel is loaded into memory and NO files are
* saved. Loading REPLACES the data in memory, so preserve/save first if needed.
* (Here the demo panel is disposable, so we don't bother.)
garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) ///
      tvarying(x1 x2) tinvariant(region) k(50)

* Useful numbers are returned:
return list

* After the run, the matched panel is in memory with these new variables:
*   weight       proximity weight (use as [pw=weight])
*   treated_id   the treated unit each control is matched to (= block_id)
*   treated_grp  1 = treated unit, 0 = selected control
*   control_sel  1 = selected control (donor)
*   garf_unit    numeric id for absorb()
*   garf_block   numeric id for cluster()
*   post         1 if year >= 2000

* --- estimation stage: garfte (new in 1.2) --------------------------------
* garfte runs on the matched panel in memory and never changes it. Results:
* a table by event time + overall ATT, an event-study graph, e(b)/e(V), and
* files in garf_output/estimations/ (incl. the effect of EACH treated unit).
* gap = treated minus its own weighted donors, in levels (replacement allowed):
garfte, method(gap)
ereturn list

* The matching DiD needs a match WITHOUT replacement (garfte checks it).
* garf replaced the data in memory, so reload the panel before matching again:
import delimited "$GARF_HOME/demo/input/panel.csv", clear varnames(1) case(preserve)
garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) ///
      tvarying(x1 x2) tinvariant(region) noreplacement k(3)
garfte, method(did)
garfte, method(gap)                 // the bracketing estimate on the same match

* Only 5 periods either side of T0, 90% intervals, no graph:
*   garfte, method(did) pre(5) post(5) level(90) nograph
* On a saved output folder (garf ..., store / noload store):
*   garfte using "garf_output", method(gap)

* --- or estimate the matching DiD yourself -------------------------------
* With the standard treatment coding, W is the treated-times-post indicator,
* so its coefficient (with unit + time FE, proximity-weighted) is the ATT.
capture which reghdfe
if _rc ssc install reghdfe, replace
reghdfe Y W [pw=weight], absorb(garf_unit year) vce(cluster garf_block)

* If instead your treatment variable marks a treated unit in ALL years, build
* the interaction yourself:
*   reghdfe Y c.treated_grp#c.post [pw=weight], absorb(garf_unit year) vce(cluster garf_block)

* --- variants -----------------------------------------------------------
* More weight on the closest donors (default wpower(1) = proportional to proximity;
* wpower(0) = uniform k:1 weights). Check r(mean_ess) and balance when raising it:
*   garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) wpower(2)
*
* Keep the CSV files too (matched_panel.csv, matches.csv, balance.csv, ...):
*   garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) store
*
* Produce files only, without replacing the data in memory:
*   garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) noload store
*
* Turn on hyperparameter tuning (slow):
*   garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) tune

* --- full design-stage diagnostics -----------------------------------------
* Balance by component, variance ratios, balance over time, outcome path,
* placebo, ESS and donor reuse, overlap, FMQ, importance, sensitivity grid
* and frontier, as commented do-files (outputs in ./garf_diagnostics):
*   do "$GARF_HOME/DIAGNOSTICS/00_master.do"
