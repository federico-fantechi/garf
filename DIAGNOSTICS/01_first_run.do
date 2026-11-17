*! 01_first_run.do  --  the default design (paper: Appendix C.3, Table 6)
*  ---------------------------------------------------------------------------
*  Saves the original panel (raw.dta), runs garf with the default design and
*  keeps both the matched panel (matched.dta) and the output files
*  (garf_output/). Every other do-file in this folder reads these three.
*  ---------------------------------------------------------------------------

if "$GARF_WORK" == "" {
    di as error "Run 00_settings.do first."
    exit 198
}
cd "$GARF_WORK"

* --- load the panel ----------------------------------------------------------
if strpos(lower("$GARF_DATA"), ".csv") {
    import delimited "$GARF_DATA", clear varnames(1) case(preserve)
}
else {
    use "$GARF_DATA", clear
}
save raw.dta, replace                    // keep the original panel

* --- default design: k = $K, no caliper, gamma = 1, with replacement ----------
garf, $GARF_CALL k($K) store             // store: also write garf_output/

save matched.dta, replace                // garf replaced the data in memory
return list                              // all values of the summary table
