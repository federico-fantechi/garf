*! version 1.2  19nov2026  garfte.ado — estimation stage for GARF (binary): matching DiD and gap
*! Federico Fantechi, PhD — University of Palermo — federico.fantechi@unipa.it
*! Runs after garf. Shells out to the packaged Python engine garf_estimate.py (same Python as garf).
*!
*!   garfte, method(did|gap) [options]                 on the matched panel garf left in memory
*!   garfte using "folder", method(did|gap) [options]  on a saved garf output folder
*!
*! did : matching difference-in-differences, base = mean gap over the matching window;
*!       needs a panel matched WITHOUT replacement; SEs clustered by block (treated unit).
*! gap : treated minus its own weighted donors, in levels; replacement allowed; SEs clustered
*!       by block plus a correction for donors shared between blocks.
*! Results: table on screen, e(b)/e(V) (event times + overall ATT), event-study graph (ATT(e)
*! line with a level()% confidence band, default 95%), and the
*! files garfte_<method>_{event,units,vcov,summary}.csv in garf_output/estimations/.
*! The data in memory are never changed.

program define garfte, eclass
    version 14.0

    syntax [using/], Method(string) ///
        [ Outcome(string) PRE(integer -1) POST(integer -1) Level(real 95) NOGraph ///
          OUTdir(string) ENGINEdir(string) PYthon(string) ]

    local method = lower(strtrim("`method'"))
    if !inlist("`method'", "did", "gap") {
        di as error "method() must be did or gap"
        di as error "  did : matching difference-in-differences (needs garf ..., noreplacement)"
        di as error "  gap : treated unit minus its own weighted donors, in levels"
        exit 198
    }
    if `level' < 10 | `level' > 99.99 {
        di as error "level() must be between 10 and 99.99 (default 95)"
        exit 198
    }
    if `pre' < -1 | `post' < -1 {
        di as error "pre() and post() must be >= 0 (omit them to keep all event times)"
        exit 198
    }

    // -------- locate the engine and Python (as garf does) -------------------
    if "`enginedir'" == "" & "$GARF_HOME" != "" local enginedir "$GARF_HOME/engine"
    if "`enginedir'" == "" {
        di as error "Set the package location first, e.g.:  global GARF_HOME " ///
            `""/path/to/GARF-binary""' "   (or pass enginedir())."
        exit 198
    }
    local engine "`enginedir'/garf_estimate.py"
    capture confirm file "`engine'"
    if _rc {
        di as error "estimation engine not found: `engine'"
        di as error "garfte needs GARF-binary version 1.2 or newer; check GARF_HOME / enginedir()."
        exit 601
    }
    if "`python'" == "" {
        if "`c(os)'" == "Windows"  local venvpy "`enginedir'/../.venv/Scripts/python.exe"
        else                        local venvpy "`enginedir'/../.venv/bin/python"
        capture confirm file "`venvpy'"
        if !_rc local python "`venvpy'"
        else if "`c(os)'" == "Windows" local python "python"
        else                            local python "python3"
    }

    // -------- source: matched panel in memory, or a saved garf folder --------
    tempfile pcsv
    local pcsv "`pcsv'.csv"
    if `"`using'"' == "" {
        local gv : char _dta[garf_version]
        if "`gv'" == "" {
            di as error "no GARF matched panel in memory."
            di as error "Run garf (version 1.2 or newer) first, or point garfte at a saved garf output"
            di as error `"folder:  garfte using "garf_output", method(`method')"'
            exit 459
        }
        local t0     : char _dta[garf_t0]
        local npre   : char _dta[garf_npre]
        local inct0  : char _dta[garf_includet0]
        local mout   : char _dta[garf_outcome]
        local idv    : char _dta[garf_id]
        local timev  : char _dta[garf_time]
        local gout   : char _dta[garf_outdir]
        if "`outcome'" == "" local outcome "`mout'"
        foreach v in `idv' `timev' block_id role donor_id weight `outcome' {
            capture confirm variable `v', exact
            if _rc {
                di as error "variable `v' not found: garfte needs the matched panel as garf left it" ///
                    " (outcome() must be a variable in memory)."
                exit 111
            }
        }
        capture confirm numeric variable `outcome', exact
        if _rc {
            di as error "outcome `outcome' must be numeric"
            exit 109
        }
        local keepv `idv' `timev' block_id role donor_id weight `outcome'
        capture confirm variable fmq, exact
        if !_rc local keepv `keepv' fmq
        local keepv : list uniq keepv
        preserve
            quietly keep `keepv'
            quietly export delimited using "`pcsv'", replace
        restore
        if "`outdir'" == "" {
            if "`gout'" != "" local outdir "`gout'/estimations"
            else              local outdir "`c(pwd)'/garf_output/estimations"
        }
        local src "matched panel in memory"
    }
    else {
        local ufwd = subinstr(`"`using'"', "\", "/", .)
        capture confirm file "`ufwd'/matched_panel.csv"
        if _rc {
            di as error `"`using'/matched_panel.csv not found: is "`using'" a garf output folder"'
            di as error "(saved with store, export() or outdir())?"
            exit 601
        }
        if "`outdir'" == "" local outdir "`ufwd'/estimations"
        local src `"folder `ufwd'"'
    }
    local outdir = subinstr("`outdir'", "\", "/", .)
    local pfx "`outdir'/garfte_`method'_"
    foreach f in event.csv units.csv vcov.csv summary.csv messages.txt error.txt {
        capture erase "`pfx'`f'"
    }

    // -------- run the estimation engine ---------------------------------------
    if `"`using'"' == "" {
        local args `"--panel "`pcsv'" --t0 `t0' --npre `npre' --includet0 `inct0' --id "`idv'" --time "`timev'" --matched-outcome "`mout'""'
    }
    else {
        local args `"--from "`ufwd'""'
    }
    if "`outcome'" != "" local args `"`args' --outcome "`outcome'""'
    local args `"`args' --method `method' --pre `pre' --post `post' --level `level' --outdir "`outdir'" --quiet"'
    if "`c(os)'" == "Windows" {
        local s_py  = subinstr(`"`python'"', "/", "\", .)
        local s_eng = subinstr("`engine'",   "/", "\", .)
    }
    else {
        local s_py  `"`python'"'
        local s_eng "`engine'"
    }
    di as txt "[garfte] method `method' on the `src' ..."
    shell "`s_py'" "`s_eng'" `args'
    capture erase "`pcsv'"

    // -------- errors from the engine --------------------------------------------
    capture confirm file "`pfx'error.txt"
    if !_rc {
        tempname eh
        file open `eh' using "`pfx'error.txt", read text
        file read `eh' line
        di as error "[garfte] estimation stopped:"
        while r(eof) == 0 {
            di as error `"  `line'"'
            file read `eh' line
        }
        file close `eh'
        exit 498
    }
    capture confirm file "`pfx'summary.csv"
    if _rc {
        di as error "[garfte] the estimation engine produced no results."
        di as error "Run the same call from a terminal to see the Python error, or check that the"
        di as error "Python interpreter (`python') has pandas, numpy and scipy installed."
        exit 601
    }

    // -------- summary (key,value) -------------------------------------------------
    tempname sh
    file open `sh' using "`pfx'summary.csv", read text
    file read `sh' line
    while r(eof) == 0 {
        local p = strpos(`"`line'"', ",")
        if `p' > 0 {
            local k = substr(`"`line'"', 1, `p' - 1)
            local s_`k' = substr(`"`line'"', `p' + 1, .)
        }
        file read `sh' line
    }
    file close `sh'

    // -------- event table, b, V, graph (data in memory untouched) ----------------
    tempname b V EV
    preserve
        quietly import delimited using "`pfx'vcov.csv", clear varnames(1) case(preserve) asdouble
        local names
        forvalues i = 1/`=_N' {
            local names `names' `=name[`i']'
        }
        quietly drop name
        mkmat _all, matrix(`V')
        matrix `V' = (`V' + `V'') / 2
        matrix rownames `V' = `names'
        matrix colnames `V' = `names'

        quietly import delimited using "`pfx'event.csv", clear varnames(1) case(preserve) asdouble
        mkmat att, matrix(`b')
        matrix `b' = `b''
        matrix colnames `b' = `names'
        mkmat e year att se t p lb ub n_treated, matrix(`EV') rownames(name)

        // ---- table -------------------------------------------------------------
        di as txt "{hline 78}"
        di as txt "GARF estimation (garfte) — method " as res "`method'" ///
            as txt "   outcome " as res "`s_outcome'" as txt "   T0 = " as res "`s_t0'"
        if "`method'" == "did" di as txt "  matching DiD: gap measured from its mean over the matching window " ///
            as res "`s_window_first'-`s_window_last'"
        else                   di as txt "  gap: treated unit minus its own weighted donors, in levels"
        di as txt "  treated units = " as res "`s_n_treated'" as txt "   distinct donors = " as res "`s_n_donors'" ///
            as txt "   donors reused = " as res "`s_donors_reused'" as txt "   t df = " as res "`s_df'"
        di as txt "{hline 78}"
        di as txt %6s "e" %7s "year" "  " %-13s "period" %10s "ATT" %10s "Std.err." %8s "P>|t|" ///
            "   [`level'% conf. interval]" %5s "N"
        di as txt "{hline 78}"
        forvalues i = 1/`=_N' {
            local nm = name[`i']
            local pr = period[`i']
            if "`nm'" == "ATT" {
                di as txt "{hline 78}"
                di as txt %13s "ATT" "  " %-13s "overall post" ///
                    as res %10.4f att[`i'] %10.4f se[`i'] %8.3f p[`i'] "  " %10.4f lb[`i'] %10.4f ub[`i'] ///
                    as txt %5.0f n_treated[`i']
            }
            else {
                di as txt %6.0f e[`i'] %7.0f year[`i'] "  " %-13s "`pr'" ///
                    as res %10.4f att[`i'] %10.4f se[`i'] %8.3f p[`i'] "  " %10.4f lb[`i'] %10.4f ub[`i'] ///
                    as txt %5.0f n_treated[`i']
            }
        }
        di as txt "{hline 78}"
        di as txt "  ATT = equal-weight average of the post-treatment ATT(e) shown."
        if "`method'" == "gap" & "`s_mean_window_gap'" != "." {
            local mwg : di %10.4f `s_mean_window_gap'
            local mwa : di %10.4f `s_att'
            di as txt "  mean gap in the matching window = " as res strtrim("`mwg'") ///
                as txt "   (overall post ATT = " as res strtrim("`mwa'") as txt ")"
            di as txt "  The gap is credible only if the window gap is small relative to the ATT;"
            di as txt "  otherwise read method(did), which removes it."
        }

        // ---- graph: ATT(e) line with a confidence band (SCM-style) ------------------
        if "`nograph'" == "" {
            quietly drop if name == "ATT"
            sort e
            quietly summarize lb
            local ymin = min(r(min), 0)
            quietly summarize ub
            local ymax = max(r(max), 0)
            local pad  = 0.06 * (`ymax' - `ymin')
            local ymin = `ymin' - `pad'
            local ymax = `ymax' + `pad'
            quietly summarize e
            local emin = r(min)
            local emax = r(max)
            local wx1 = max(`s_window_first' - `s_t0' - 0.5, `emin')
            local wx2 = min(`s_window_last'  - `s_t0' + 0.5, `emax')
            local onset = `s_window_last' - `s_t0' + 0.5
            // overall ATT for the subtitle
            local a_att : di %9.3f `s_att'
            local a_lb  : di %9.3f `s_lb'
            local a_ub  : di %9.3f `s_ub'
            local a_att = strtrim("`a_att'")
            local a_lb  = strtrim("`a_lb'")
            local a_ub  = strtrim("`a_ub'")
            // colours: Stata 15+ takes the line and band colour from the current scheme
            // (plot style p1) and draws band and window with transparency; Stata 14 has
            // no transparency, so it falls back to fixed light colours
            if c(stata_version) >= 15 {
                local vpre "version 15:"
                local wsty "fcolor(gs8%12) lwidth(none)"
                local bsty "pstyle(p1) fcolor(%30) lwidth(none)"
                local lsty "pstyle(p1) lwidth(medthick)"
            }
            else {
                local vpre
                local wsty "fcolor(gs15) lcolor(gs15) lwidth(none)"
                local bsty "fcolor(eltblue) lcolor(eltblue) lwidth(none)"
                local lsty "lcolor(navy) lwidth(medthick)"
            }
            local plots
            local np = 0
            local wleg
            if `wx2' > `wx1' {
                local plots `plots' (function y = `ymax', range(`wx1' `wx2') recast(area) ///
                    base(`ymin') `wsty')
                local ++np
                local wleg `np'
            }
            local plots `plots' (rarea lb ub e, `bsty' cmissing(n))
            local ++np
            local bleg `np'
            local plots `plots' (line att e, `lsty' cmissing(n))
            local ++np
            local lleg `np'
            local leg `"`lleg' "ATT(e)" `bleg' "`level'% confidence band""'
            if "`wleg'" != "" local leg `"`leg' `wleg' "matching window""'
            if "`method'" == "did" {
                local ttl "GARF matching DiD"
                local sub "Effect on `s_outcome', measured from the matching-window mean"
            }
            else {
                local ttl "GARF gap"
                local sub "Effect on `s_outcome': treated minus its own weighted donors, in levels"
            }
            capture noisily `vpre' twoway `plots', ///
                xline(`onset', lcolor(gs6) lpattern(dash) lwidth(thin)) ///
                yline(0, lcolor(gs9) lwidth(thin)) ///
                title("`ttl'", size(medlarge) position(11) justification(left)) ///
                subtitle("`sub'" "Overall ATT (post average) = `a_att'   [`a_lb', `a_ub']", ///
                    size(small) position(11) justification(left)) ///
                xtitle("Periods since treatment (T0 = `s_t0')", size(small)) ///
                ytitle("ATT(e)", size(small)) ///
                xlabel(, labsize(small) nogrid) ylabel(, labsize(small) angle(0) nogrid) ///
                legend(order(`leg') rows(1) position(6) size(small) symxsize(8) ///
                    region(lstyle(none) fcolor(none))) ///
                graphregion(margin(medsmall)) plotregion(margin(small)) ///
                xsize(7) ysize(4.5) name(garfte_`method', replace)
            if _rc {
                di as txt "[garfte] note: the graph could not be drawn (error " _rc "); results are unaffected."
                local nograph "nograph"
            }
            else {
                quietly graph save "`outdir'/garfte_`method'.gph", replace
                capture graph export "`outdir'/garfte_`method'.pdf", replace
                if _rc di as txt "[garfte] note: could not export the PDF (the .gph was saved)."
            }
        }
    restore

    // -------- notes and warnings --------------------------------------------------------
    capture confirm file "`pfx'messages.txt"
    if !_rc {
        tempname mh
        file open `mh' using "`pfx'messages.txt", read text
        file read `mh' line
        local iswarn 0
        while r(eof) == 0 {
            // a message starts at column 1; its wrapped continuation lines start with "  "
            if substr(`"`line'"', 1, 2) != "  " local iswarn = (substr(`"`line'"', 1, 8) == "warning:")
            if `iswarn' di as error `"`line'"'
            else        di as txt   `"`line'"'
            file read `mh' line
        }
        file close `mh'
    }
    di as txt "[garfte] files written to " as res "`outdir'/" as txt ":"
    di as txt "         garfte_`method'_event.csv  (ATT by event time and overall)"
    di as txt "         garfte_`method'_units.csv  (effect of each treated unit in each period)"
    di as txt "         garfte_`method'_vcov.csv   garfte_`method'_summary.csv"
    if "`nograph'" == "" di as txt "         garfte_`method'.gph / .pdf  (event-study graph)"
    di as txt "[garfte] the data in memory were not changed."

    // -------- stored results ---------------------------------------------------------------
    ereturn post `b' `V', dof(`s_df')
    ereturn scalar att             = `s_att'
    ereturn scalar se_att          = `s_se'
    ereturn scalar lb_att          = `s_lb'
    ereturn scalar ub_att          = `s_ub'
    ereturn scalar p_att           = `s_p'
    ereturn scalar level           = `level'
    ereturn scalar N_treated       = `s_n_treated'
    ereturn scalar N_donors        = `s_n_donors'
    ereturn scalar donors_reused   = `s_donors_reused'
    ereturn scalar t0              = `s_t0'
    ereturn scalar npre            = `s_npre'
    ereturn scalar window_first    = `s_window_first'
    ereturn scalar window_last     = `s_window_last'
    ereturn scalar mean_window_gap = `s_mean_window_gap'
    ereturn matrix event           = `EV'
    ereturn local  outdir          "`outdir'"
    ereturn local  outcome         "`s_outcome'"
    ereturn local  matched_outcome "`s_matched_outcome'"
    ereturn local  method          "`method'"
    ereturn local  properties      "b V"
    ereturn local  title           "GARF estimation (`method')"
    ereturn local  cmd             "garfte"
end
