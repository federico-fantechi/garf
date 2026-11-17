GARF-binary version 1.1 -- DIAGNOSTICS
===================================

Commented do-files for the full design-stage diagnostic protocol of the GARF
working paper (Section 5; the same code is printed in its Appendix C). Running
00_master.do on the demonstration panel reproduces every table and figure of
that section. Every file uses pre-treatment information only: post-treatment
outcomes belong to the analysis stage and are never plotted or summarised here.

HOW TO RUN (demonstration panel)
  1. Once per machine: GARF_HOME set in profile.do, garf.ado/garf.sthlp version 1.1
     on the ado-path, and the demo panel generated (demo/run_demo.sh or .bat).
  2. In Stata, cd to a project folder (NOT the package folder), then
         do "$GARF_HOME/DIAGNOSTICS/00_master.do"
  3. Outputs are written to ./garf_diagnostics: raw.dta, matched.dta,
     garf_output/, runs/, sens.dta, and the five figures (.pdf).
  Step 11 (the sensitivity grid) fits 19 forests: a few minutes on the demo,
  hours on tens of thousands of units (each run builds features, fits the
  forest and computes proximities from scratch).

YOUR OWN DATA
  Copy this DIAGNOSTICS folder to your project and edit only 00_settings.do
  (data path, variable names, T0, window L, default k, placebo backdating H).
  Then, in Stata:
         global GARF_DIAG "/path/to/your/copy/DIAGNOSTICS"
         do "$GARF_DIAG/00_master.do"
  Label positions in the figures are computed from the data; on very
  different data, 12_frontier.do may need its label positions adjusted.
  The widest design of the grid is KWIDE, default k(1000): equal to the
  uncapped forest kernel on pools below 1,000 untreated units, and safe on
  large pools, where nokcap can give thousands of donors per treated unit.

RUNNING ONE STEP
  do "$GARF_DIAG/00_settings.do"     (once per session)
  do "$GARF_DIAG/01_first_run.do"    (once: creates the files the others read)
  do "$GARF_DIAG/0X_....do"          (any step)

PITFALL
  Stata uses a garf.ado in the CURRENT folder before the ado-path. Keep the
  working folder free of old copies; 00_settings.do runs  discard  and
  which garf  (first line must read version 1.1).

FILES                              PAPER (Sect. 5)   OUTPUT
  00_settings.do   settings         --                --
  00_master.do     runs all         --                --
  01_first_run.do  default design   App. C.3          raw.dta, matched.dta, garf_output/
  02_balance_features.do            feature balance   table in the log
  03_variance_ratios.do             variance ratios   table in the log
  04_balance_over_time.do           balance over time fig_balance_time.pdf
  05_outcome_path.do                outcome path      fig_outcome_path.pdf
  06_placebo.do                     placebo           fig_placebo.pdf
  07_ess_reuse.do                   ESS and reuse     table in the log
  08_overlap.do                     overlap           fig_overlap.pdf
  09_fmq_blocks.do                  FMQ by block      table in the log
  10_forest_importance.do           forest fit        table in the log
  11_sensitivity_grid.do            sensitivity       sens.dta, runs/
  12_frontier.do                    FMQ-ESS frontier  fig_frontier.pdf
