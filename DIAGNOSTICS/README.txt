GARF-binary version 1.2 -- DIAGNOSTICS
===================================

Commented do-files for the full design-stage diagnostic protocol of the GARF
working paper (Section 5; the same code is printed in its Appendix C). Running
00_master.do on the demonstration panel reproduces every table and figure of
that section. Every file uses pre-treatment information only: post-treatment
outcomes belong to the analysis stage (the command garfte, new in 1.2) and are
never plotted or summarised here.

HOW TO RUN AFTER YOUR OWN garf  (nothing to edit)
  garf (1.2 or newer) stores its whole design on the matched panel: variables,
  T0, window, k, caliper, wpower, replacement, seed, forest options, the if/in
  sample and the file it ran on (char list _dta[]). The do-files read it.
  1. Run garf on your data; store is recommended (keeps garf's files):
         use mypanel.dta, clear
         garf, treatment(...) outcome(...) id(...) time(...) t0(...) ... store
  2. With the matched panel still in memory (or after use garf_matched.dta),
     from your project folder:
         do "$GARF_HOME/DIAGNOSTICS/00_master.do"
  3. Outputs go to ./garf_diagnostics; at the end Stata is back in your
     folder with the matched panel in memory.
  Notes
  - garf records the .dta the data were loaded from (sorting or xtset after
    loading are fine). 01_first_run.do checks that the file still has garf's
    variables and the same number of rows. If not (CSV import, variables
    created or rows dropped after loading), save the panel garf ran on and
    set first:   global GARF_DATA "/path/to/panel.dta"
  - Without store/export(), 01_first_run.do reruns your design once, with the
    same seed, to recreate garf's files (same data and seed: same match).
  - 06_placebo backdates YOUR design by H = 3 periods, fewer if the data do
    not reach back far enough for a full window (it says so), and skips if
    no backdating fits.
  - The sensitivity grid (11-12: the 19 designs below plus yours) is skipped
    on pools of 20,000 units or more:  global GARF_GRID on  to force it,
    off  to skip it. It uses your forest options, not tuning.
  - The settings stay in force for the session, so single steps can be
    rerun. To start over (another run, or the demo):  macro drop GARF_*
    (then set GARF_HOME again).

HOW TO RUN (demonstration panel)
  1. Once per machine: GARF_HOME set in profile.do, garf.ado/garf.sthlp version 1.2
     on the ado-path, and the demo panel generated (demo/run_demo.sh or .bat).
  2. In Stata, with NO matched panel in memory, cd to a project folder (NOT
     the package folder), then
         do "$GARF_HOME/DIAGNOSTICS/00_master.do"
     The manual settings of 00_settings.do (the demo panel) are used.
  3. Outputs are written to ./garf_diagnostics: raw.dta, matched.dta,
     garf_output/, runs/, sens.dta, and the five figures (.pdf).
  Step 11 (the sensitivity grid) fits 19 forests: a few minutes on the demo,
  hours on tens of thousands of units (each run builds features, fits the
  forest and computes proximities from scratch).

SETTINGS BY HAND
  Copy this DIAGNOSTICS folder to your project and edit the manual-settings
  block of 00_settings.do. Then, in Stata, with no matched panel in memory:
         global GARF_DIAG "/path/to/your/copy/DIAGNOSTICS"
         do "$GARF_DIAG/00_master.do"
  Label positions in the figures are computed from the data; on very
  different data, 12_frontier.do may need its label positions adjusted.
  The widest design of the grid is KWIDE, default k(1000) (global KWIDE to
  change it): equal to the uncapped forest kernel on pools below 1,000
  untreated units, and safe on large pools, where nokcap can give thousands
  of donors per treated unit.

RUNNING ONE STEP
  do "$GARF_DIAG/00_settings.do"     (once per session)
  do "$GARF_DIAG/01_first_run.do"    (once: creates the files the others read)
  do "$GARF_DIAG/0X_....do"          (any step)

PITFALL
  Stata uses a garf.ado in the CURRENT folder before the ado-path. Keep the
  working folder free of old copies; 00_settings.do runs  discard  and
  which garf  (first line must read version 1.2).

FILES                              PAPER (Sect. 5)   OUTPUT
  00_settings.do   settings (from the matched panel, or manual)
  00_master.do     runs all         --                --
  01_first_run.do  design checked   App. C.3          raw.dta, matched.dta, garf's files
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
