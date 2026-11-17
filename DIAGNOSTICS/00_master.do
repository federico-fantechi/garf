*! 00_master.do  --  GARF diagnostics, version 1.1
*  ---------------------------------------------------------------------------
*  Runs the whole diagnostic protocol of Section 5 of the GARF paper, in order.
*  cd to your project folder (outputs go to ./garf_diagnostics), then:
*      do "$GARF_HOME/DIAGNOSTICS/00_master.do"
*  On a copy of this folder elsewhere: set global GARF_DIAG to it first.
*  Each step can also be run on its own after 00_settings.do and 01_first_run.do.
*  ---------------------------------------------------------------------------

if "$GARF_DIAG" == "" global GARF_DIAG "$GARF_HOME/DIAGNOSTICS"
do "$GARF_DIAG/00_settings.do"

do "$GARF_DIAG/01_first_run.do"          // default design: raw.dta, matched.dta, garf_output/
do "$GARF_DIAG/02_balance_features.do"   // balance by variable and component
do "$GARF_DIAG/03_variance_ratios.do"    // variance ratios, raw and matched
do "$GARF_DIAG/04_balance_over_time.do"  // SMD by period          -> fig_balance_time.pdf
do "$GARF_DIAG/05_outcome_path.do"       // pre-treatment path      -> fig_outcome_path.pdf
do "$GARF_DIAG/06_placebo.do"            // out-of-window placebo   -> fig_placebo.pdf
do "$GARF_DIAG/07_ess_reuse.do"          // ESS and donor reuse
do "$GARF_DIAG/08_overlap.do"            // overlap                 -> fig_overlap.pdf
do "$GARF_DIAG/09_fmq_blocks.do"         // forest match quality by block
do "$GARF_DIAG/10_forest_importance.do"  // forest fit and top features
do "$GARF_DIAG/11_sensitivity_grid.do"   // 19 designs              -> sens.dta  (slow)
do "$GARF_DIAG/12_frontier.do"           // FMQ-ESS frontier        -> fig_frontier.pdf

di as txt _n "All diagnostics done. Outputs are in: $GARF_WORK"
