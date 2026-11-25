*! 00_master.do  --  GARF diagnostics, version 1.2
*  ---------------------------------------------------------------------------
*  Runs the whole diagnostic protocol of Section 5 of the GARF paper, in order.
*
*  AFTER YOUR OWN garf (nothing to edit): with garf's matched panel in memory
*  (or after  use garf_matched.dta), cd to your project folder and run
*      do "$GARF_HOME/DIAGNOSTICS/00_master.do"
*  The settings are read from the matched panel (see 00_settings.do).
*  DEMONSTRATION PANEL: the same command with no matched panel in memory.
*  Outputs go to ./garf_diagnostics; Stata returns to your folder at the end.
*  Each step can also be run on its own after 00_settings.do and 01_first_run.do.
*  ---------------------------------------------------------------------------

local here "`c(pwd)'"
if "$GARF_DIAG" == "" global GARF_DIAG "$GARF_HOME/DIAGNOSTICS"
do "$GARF_DIAG/00_settings.do"

do "$GARF_DIAG/01_first_run.do"          // raw.dta, matched.dta, garf's files
do "$GARF_DIAG/02_balance_features.do"   // balance by variable and component
do "$GARF_DIAG/03_variance_ratios.do"    // variance ratios, raw and matched
do "$GARF_DIAG/04_balance_over_time.do"  // SMD by period          -> fig_balance_time.pdf
do "$GARF_DIAG/05_outcome_path.do"       // pre-treatment path      -> fig_outcome_path.pdf
do "$GARF_DIAG/06_placebo.do"            // out-of-window placebo   -> fig_placebo.pdf
do "$GARF_DIAG/07_ess_reuse.do"          // ESS and donor reuse
do "$GARF_DIAG/08_overlap.do"            // overlap                 -> fig_overlap.pdf
do "$GARF_DIAG/09_fmq_blocks.do"         // forest match quality by block
do "$GARF_DIAG/10_forest_importance.do"  // forest fit and top features
if "$GARF_GRID" == "on" {
    do "$GARF_DIAG/11_sensitivity_grid.do"   // 19-20 designs       -> sens.dta  (slow)
    do "$GARF_DIAG/12_frontier.do"           // FMQ-ESS frontier    -> fig_frontier.pdf
}
else {
    di as txt _n "[diagnostics] sensitivity grid (11-12) skipped. To run it:"
    di as txt `"    global GARF_GRID on"'
    di as txt `"    do "$GARF_DIAG/11_sensitivity_grid.do""'
    di as txt `"    do "$GARF_DIAG/12_frontier.do""'
}

use "$GARF_WORK/matched.dta", clear      // leave the matched panel in memory
quietly cd "`here'"
di as txt _n "All diagnostics done. Outputs are in: $GARF_WORK"
di as txt "The matched panel is back in memory (garfte can run on it)."
