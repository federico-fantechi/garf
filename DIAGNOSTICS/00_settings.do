*! 00_settings.do  --  GARF diagnostics, version 1.1
*  ---------------------------------------------------------------------------
*  Settings shared by every do-file in this folder. 00_master.do runs it first;
*  to run a single step, run this file and 01_first_run.do once, then the step.
*  Edit ONLY this file to move from the demonstration panel to your own data
*  (best on a copy of the DIAGNOSTICS folder outside the package: see README).
*  ---------------------------------------------------------------------------

version 14.0

* --- where things are ------------------------------------------------------
* GARF_HOME is the package folder (normally set once in profile.do).
if "$GARF_HOME" == "" {
    di as error "Set the package location first, e.g.  global GARF_HOME " ///
        `""/path/to/GARF-binary""'
    exit 198
}
* Folder with these do-files. Default: the copy shipped with the package.
* If you work on a copy elsewhere, set  global GARF_DIAG "/path/to/copy"
* BEFORE running 00_master.do.
if "$GARF_DIAG" == "" global GARF_DIAG "$GARF_HOME/DIAGNOSTICS"
* Working folder: raw.dta, matched.dta, garf_output/, runs/, sens.dta and the
* figures go here. Default: a subfolder of the folder Stata is in when you
* start (your project folder), never inside the package. It must not contain
* an old copy of garf.ado (Stata uses it before the ado-path).
* Set once per session; to switch project:  macro drop GARF_WORK  and re-run.
if "$GARF_WORK" == "" global GARF_WORK "`c(pwd)'/garf_diagnostics"
capture mkdir "$GARF_WORK"

* --- data --------------------------------------------------------------------
* A long panel (one row per unit and period), .csv or .dta
global GARF_DATA "$GARF_HOME/demo/input/panel.csv"

* --- variables and design ----------------------------------------------------
global Y    "Y"          // outcome
global W    "W"          // binary treatment (1 in treated periods, or ever-treated)
global ID   "id"         // unit identifier (numeric or string)
global T    "year"       // time variable
global T0   2000         // treatment onset
global L    5            // matching window length (periods T0-L ... T0-1)
global XTV  "x1 x2"      // time-varying covariates
global XTI  "region"     // time-invariant covariates (may be empty)
global K    50           // default cap on donors
global H    3            // placebo: periods by which the design is backdated
* Widest design of the sensitivity grid (11_sensitivity_grid.do). Default: a cap
* of 1,000 donors, which on small pools (fewer than 1,000 untreated units) is the
* same as no cap, i.e. the forest kernel. The uncapped "nokcap" keeps EVERY donor
* with positive proximity: on pools of tens of thousands of untreated units that
* can mean thousands of donors per treated unit and a matched panel of tens of
* millions of rows. Use  global KWIDE "nokcap"  only for small or medium pools.
global KWIDE "k(1000)"

* --- the common part of every garf call --------------------------------------
global GARF_CALL "treatment($W) outcome($Y) id($ID) time($T) t0($T0) npre($L)"
if "$XTV" != "" global GARF_CALL "$GARF_CALL tvarying($XTV)"
if "$XTI" != "" global GARF_CALL "$GARF_CALL tinvariant($XTI)"

* --- sanity check: the right garf is loaded -----------------------------------
discard                                  // drop any garf already in memory
which garf                               // first line must read: version 1.1
