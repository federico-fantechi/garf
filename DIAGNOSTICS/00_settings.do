*! 00_settings.do  --  GARF diagnostics, version 1.2
*  ---------------------------------------------------------------------------
*  Settings shared by every do-file in this folder. 00_master.do runs it first;
*  to run a single step, run this file and 01_first_run.do once, then the step.
*
*  TWO WAYS TO RUN THE DIAGNOSTICS
*  (A) AFTER YOUR OWN garf RUN -- nothing to edit.
*      Run garf (1.2 or newer) on your data, then, with the matched panel still
*      in memory (or after  use garf_matched.dta):
*          do "$GARF_HOME/DIAGNOSTICS/00_master.do"
*      Every setting below is read from the matched panel: garf stores its
*      design on it (char list _dta[]). The diagnostics then describe YOUR
*      design. garf records the .dta the data were loaded from; 01_first_run.do
*      checks that it still has garf's variables and rows. If not (data
*      imported from CSV, variables created or rows dropped after loading),
*      save the panel garf ran on and set the data path first:
*          global GARF_DATA "/path/to/the/panel.dta"
*  (B) THE DEMONSTRATION PANEL, or settings typed by hand: with no matched
*      panel in memory the "manual settings" block below is used. To use it
*      on your own data, edit that block (best on a copy of this folder, see
*      README.txt).
*  The settings of (A) stay in force for the rest of the session, so single
*  steps can be rerun. To start over (other run, or the demo):
*      macro drop GARF_*      (then set GARF_HOME again)
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
if "$GARF_DIAG" == "" global GARF_DIAG "$GARF_HOME/DIAGNOSTICS"
* Working folder: raw.dta, matched.dta, garf_output/, runs/, sens.dta and the
* figures go here. Default: a subfolder of the folder Stata is in (your project
* folder), never inside the package. It must not contain an old garf.ado.
if "$GARF_WORK" == "" global GARF_WORK "`c(pwd)'/garf_diagnostics"
capture mkdir "$GARF_WORK"

* =============================================================================
* (A) settings read from a garf matched panel in memory
* =============================================================================
local prevmode "$GARF_MODE"
local gv : char _dta[garf_version]
local gcall : char _dta[garf_call]
if "`gv'" != "" & "`gcall'" == "" {
    di as error "The matched panel in memory was made by an older garf: rerun garf 1.2 (or newer),"
    di as error "or clear the data to use the manual settings."
    exit 198
}
if "`gcall'" != "" {
    global GARF_MODE   "garf"
    global Y    : char _dta[garf_outcome]
    global W    : char _dta[garf_treatment]
    global ID   : char _dta[garf_id]
    global T    : char _dta[garf_time]
    global T0   : char _dta[garf_t0]
    global L    : char _dta[garf_npre]
    global XTV  : char _dta[garf_tvarying]
    global XTI  : char _dta[garf_tinvariant]
    global INCT0 : char _dta[garf_includet0]
    global GARF_CALL   : char _dta[garf_call]
    global GARF_DESIGN : char _dta[garf_design]
    global GARF_FOREST : char _dta[garf_forest]
    global GARF_TUNE   : char _dta[garf_tune]
    global GARF_IF     : char _dta[garf_if]
    global GARF_IN     : char _dta[garf_in]
    local  gout        : char _dta[garf_outdir]
    local  gsrc        : char _dta[garf_source]
    global GARF_NOBS   : char _dta[garf_nobs]
    * the data garf ran on: the file garf recorded; else GARF_DATA, if you set it
    if "`prevmode'" == "manual" global GARF_DATA ""      // not the demo path
    if "`gsrc'" != "" global GARF_DATA "`gsrc'"
    if "$GARF_DATA" == "" {
        di as error "garf cannot tell which file it ran on (the data were not loaded from a .dta,"
        di as error "e.g. imported from CSV). Set the path of a saved copy of the panel garf ran on"
        di as error "(with every variable garf used), then rerun:"
        di as error `"    global GARF_DATA "/path/to/the/panel.dta""'
        di as error "Next time, save the data just before garf: garf then records the file."
        exit 198
    }
    * the files garf saved (store / export() / outdir()), else 01 recreates them
    global GARF_OUTDIR ""
    if "`gout'" != "" {
        capture confirm file "`gout'/balance.csv"
        if !_rc global GARF_OUTDIR "`gout'"
    }
    * keep the matched panel: the next steps load other data
    quietly save "$GARF_WORK/matched.dta", replace
    global GARF_MATCHED_SAVED 1
}
else if "$GARF_MODE" == "garf" {
    * settings read earlier in this session from a garf run: keep them
    di as txt "[diagnostics] using the settings of the garf run read earlier in this session"
    di as txt "              (to start over:  macro drop GARF_*  and set GARF_HOME again)"
}
else {
* =============================================================================
* (B) manual settings: the demonstration panel, or your data typed by hand
* =============================================================================
    global GARF_MODE "manual"
    global GARF_NOBS ""
    * A long panel (one row per unit and period), .csv or .dta
    global GARF_DATA "$GARF_HOME/demo/input/panel.csv"
    global Y    "Y"          // outcome
    global W    "W"          // binary treatment (1 in treated periods, or ever-treated)
    global ID   "id"         // unit identifier (numeric or string)
    global T    "year"       // time variable
    global T0   2000         // treatment onset
    global L    5            // matching window length (periods T0-L ... T0-1)
    global XTV  "x1 x2"      // time-varying covariates
    global XTI  "region"     // time-invariant covariates (may be empty)
    global INCT0 0           // 1 if the window ends at T0 (garf option includet0)
    global K    50           // the design checked: a cap of K donors, no caliper
    global GARF_DESIGN "k($K)"
    global GARF_FOREST ""    // forest options (garf defaults)
    global GARF_TUNE   ""
    global GARF_IF     ""
    global GARF_IN     ""
    global GARF_OUTDIR ""    // 01_first_run.do runs the design and stores its files
    global GARF_CALL "treatment($W) outcome($Y) id($ID) time($T) t0($T0) npre($L)"
    if "$XTV" != "" global GARF_CALL "$GARF_CALL tvarying($XTV)"
    if "$XTI" != "" global GARF_CALL "$GARF_CALL tinvariant($XTI)"
    if $INCT0       global GARF_CALL "$GARF_CALL includet0"
    global GARF_MATCHED_SAVED 0
}

* --- derived settings (both ways) ----------------------------------------------
global WLAST  = $T0 - 1 + $INCT0         // last period of the matching window
global WFIRST = $WLAST - $L + 1          // first period of the matching window
* Placebo (06): periods by which the design is backdated. Default 3; 06 lowers it
* when the data do not go back far enough for a full window before the onset.
if "$H" == "" global H 3
* Sensitivity grid (11-12): 19 garf runs. "auto" runs it in the demo/manual mode
* and, after your own garf, only on pools below 20,000 units (01 decides).
* Set  global GARF_GRID "on"  or "off"  before 00_master.do to force it.
if "$GARF_GRID" == "" global GARF_GRID "auto"
* Widest design of the grid: a cap of 1,000 donors (= the uncapped forest kernel
* on pools below 1,000 untreated units). "nokcap" keeps every donor with positive
* proximity: use it only on small or medium pools.
if "$KWIDE" == "" global KWIDE "k(1000)"

* --- report ------------------------------------------------------------------------
di as txt _n "{hline 72}"
if "$GARF_MODE" == "garf" di as txt "GARF diagnostics -- settings read from the garf matched panel"
else                      di as txt "GARF diagnostics -- manual settings (demonstration panel unless edited)"
di as txt "  data      : $GARF_DATA"
if `"$GARF_IF$GARF_IN"' != "" di as txt `"  sample    : $GARF_IF $GARF_IN"'
di as txt "  outcome $Y   treatment $W   id $ID   time $T"
di as txt "  T0 = $T0   matching window $WFIRST-$WLAST"
di as txt "  covariates: $XTV $XTI"
di as txt "  design    : $GARF_DESIGN $GARF_FOREST"
if "$GARF_TUNE" != "" di as txt "  (the run was tuned: placebo and grid use the forest options above, not tuning)"
di as txt "  work folder: $GARF_WORK"
di as txt "{hline 72}"

* --- sanity check: the right garf is loaded -----------------------------------
discard                                  // drop any garf already in memory
which garf                               // first line must read: version 1.2
