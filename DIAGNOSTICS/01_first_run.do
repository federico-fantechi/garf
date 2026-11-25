*! 01_first_run.do  --  the design checked (paper: Appendix C.3, Table 6)
*  ---------------------------------------------------------------------------
*  Prepares the three things every other do-file reads:
*    raw.dta      the panel garf ran on (same sample: garf's if/in applied)
*    matched.dta  the matched panel of the design
*    the files garf writes (balance.csv, matches.csv, blocks.csv,
*                 diagnostics.json), in $GARF_OUTDIR
*  After your own garf (settings read from the matched panel): your matched
*  panel is used as it is, and garf's files are taken from the folder it saved
*  them to (store / export() / outdir()). If it saved none, your design is run
*  once more, with the same seed, to recreate them (same data, same seed: the
*  same match; tip: run garf with  store  next time).
*  Manual settings (demo): runs the design $GARF_DESIGN with store.
*  ---------------------------------------------------------------------------

if "$GARF_WORK" == "" {
    di as error "Run 00_settings.do first."
    exit 198
}
cd "$GARF_WORK"

* --- the panel garf ran on ------------------------------------------------------
if strpos(lower("$GARF_DATA"), ".csv") {
    import delimited "$GARF_DATA", clear varnames(1) case(preserve) asdouble
}
else {
    use "$GARF_DATA", clear
}
if `"$GARF_IF$GARF_IN"' != "" keep $GARF_IF $GARF_IN      // garf's sample
* the same data garf ran on? (variables created or rows dropped after loading
* the file are not in it: then save that version and set GARF_DATA to it)
capture confirm variable $Y $W $ID $T $XTV $XTI
if _rc {
    di as error "$GARF_DATA lacks some of the variables garf used ($Y $W $ID $T $XTV $XTI)."
    di as error "They were probably created after the file was loaded. Recreate them, save that"
    di as error `"panel, set  global GARF_DATA "/path/to/it.dta"  and rerun 00_master.do."'
    exit 111
}
if "$GARF_NOBS" != "" & _N != 0$GARF_NOBS {
    di as error "$GARF_DATA has " _N " rows in garf's sample, garf ran on $GARF_NOBS."
    di as error "The data were changed after loading (rows dropped or added). Save the panel"
    di as error `"garf ran on, set  global GARF_DATA "/path/to/it.dta"  and rerun 00_master.do."'
    exit 459
}
save raw.dta, replace

* --- sensitivity grid: decide "auto" ---------------------------------------------
if "$GARF_GRID" == "auto" {
    tempvar first
    bysort $ID: gen byte `first' = _n == 1
    quietly count if `first'
    local nunits = r(N)
    drop `first'
    if "$GARF_MODE" == "manual" | `nunits' < 20000 global GARF_GRID "on"
    else {
        global GARF_GRID "off"
        di as txt "[diagnostics] `nunits' units: the sensitivity grid (11-12, 20 garf runs) is"
        di as txt "              skipped by 00_master.do. To run it:  global GARF_GRID on"
    }
}

* --- matched panel and garf's files ---------------------------------------------------
if "$GARF_MODE" == "garf" & "$GARF_OUTDIR" != "" {
    di as txt "[diagnostics] using your matched panel and the files of your garf run in"
    di as txt "              $GARF_OUTDIR"
}
else if "$GARF_MODE" == "garf" {
    di as txt "[diagnostics] your garf run saved no files: rerunning your design once (same seed,"
    di as txt "              same match) to recreate them in $GARF_WORK/garf_output."
    di as txt "              Next time, run garf with  store  to skip this step."
    garf, $GARF_CALL $GARF_DESIGN $GARF_FOREST $GARF_TUNE export("$GARF_WORK/garf_output")
    save matched.dta, replace            // the panel that goes with these files
    global GARF_OUTDIR "$GARF_WORK/garf_output"
}
else {
    * manual settings: run the design, keep the files (store -> ./garf_output)
    garf, $GARF_CALL $GARF_DESIGN $GARF_FOREST store
    save matched.dta, replace            // garf replaced the data in memory
    global GARF_OUTDIR "$GARF_WORK/garf_output"
    return list                          // all values of the summary table
}
