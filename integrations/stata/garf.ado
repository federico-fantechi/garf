*! version 1.1  17nov2026  garf.ado — Stata wrapper for the GARF (binary) forest-proximity matching engine
*! Federico Fantechi, PhD — University of Palermo — federico.fantechi@unipa.it
*! Shells out to the packaged Python engine (garf_binary.py); no Python-in-Stata setup needed.
*! Requires: Python 3.10+ and the GARF-binary package (with its .venv or a Python that has the deps).
*!
*! Point the wrapper at the package once, e.g.:
*!     global GARF_HOME "C:/Users/you/GARF-binary"      // Windows
*!     global GARF_HOME "/Users/you/GARF-binary"        // macOS / Linux
*!
*! Donor selection: k nearest (k()) and/or a proximity caliper (caliper()); they compose.
*! beta 1.1: caliper() now filters on MIN-MAX NORMALISED proximity (range [0,1]; see
*! proximity_minmax in the outputs and `help garf'), not the raw kernel used in beta 1.0.
*! beta 1.2: wpower(#) sets donor weights ~ proximity^# (default 1 = as in beta 1.1;
*! 0 = uniform k:1 weights; >1 = more weight on the closest donors). Engine fixes:
*! noreplacement now matches the hardest-to-match treated units first; the matched SMD
*! now uses only the matched treated units. New Forest Match Quality (FMQ): r(fmq_mean),
*! r(fmq_p10), r(fmq_mean_all); per-pair cooccurrence and per-block fmq in the outputs.
*! beta 1.3: the trajectory descriptor drops the duplicated end-points (L+7 features per
*! time-varying variable) and uses symmetric growth; out-of-bag balanced accuracy and AUC
*! are returned as r(oob_bacc) and r(oob_auc) (in-sample r(auc) removed); the legacy
*! options minprox() and load were removed.
*! beta 1.4: no change to the command; the package ships DIAGNOSTICS/, commented do-files
*! for the full diagnostic protocol (see help garf, section Diagnostics).
*! beta 1.5: FMQ is normalised by the best ATTAINABLE donor (1 = every donor as close as the
*! closest untreated unit in the data), so its level no longer depends on N and minleaf();
*! the beta 1.4 twin-scaled version is returned as r(fmq_twin_mean). r(fmq_undefined) counts
*! blocks where no control is closer than the average one. DIAGNOSTICS: grid's widest design
*! is the setting KWIDE, default k(1000).
*! version 1.1 (first public release; = development build beta 1.5 plus): the config file
*! passed to the engine quotes all names and paths (variables called on/off/yes/no/true/false/
*! null, and paths containing ' #' or ': ', now work); the data are exported with nolabel
*! (a value-labelled treatment no longer breaks the engine); the temporary copy of the data
*! and the temporary config file are erased after the engine has run.
*! By DEFAULT the matched panel is loaded into memory (no files saved); pass store or
*! export() to also write the CSV outputs; noload to skip loading.

program define garf, rclass
    version 14.0

    syntax [if] [in], ///
        TREatment(varname) Outcome(varname) ID(varname) TIme(varname) T0(integer) ///
        [ TVarying(varlist) TINvariant(varlist) NPre(integer 5) INCLUDEt0 ///
          K(integer 50) NOKcap CALiper(real 0) NOReplacement WPower(real 1) ///
          TUNE TUNINGmetric(string) NITer(integer 40) CVfolds(integer 5) ///
          Nestimators(integer 1000) MINLeaf(integer 5) MAXFeatures(string) MAXDepth(string) ///
          SEed(integer 42) ENGINEdir(string) PYthon(string) ///
          EXPORT(string) STORE OUTdir(string) NOLoad VERBose ]

    // -------- defaults for string options -----------------------------------
    if "`tuningmetric'" == "" local tuningmetric "combined"
    if "`maxfeatures'"  == "" local maxfeatures  "sqrt"
    if "`maxdepth'"     == "" local maxdepth     "null"

    if `wpower' < 0 {
        di as error "wpower() must be >= 0 (1 = proportional to proximity, 0 = uniform weights)"
        exit 198
    }
    if !inlist("`tuningmetric'", "accuracy", "f1_treated", "combined") {
        di as error "tuningmetric() must be accuracy, f1_treated, or combined"
        exit 198
    }

    // -------- locate the engine ---------------------------------------------
    if "`enginedir'" == "" & "$GARF_HOME" != "" local enginedir "$GARF_HOME/engine"
    if "`enginedir'" == "" {
        di as error "Set the package location first, e.g.:  global GARF_HOME " ///
            `""/path/to/GARF-binary""' "   (or pass enginedir())."
        exit 198
    }
    local engine "`enginedir'/garf_binary.py"
    capture confirm file "`engine'"
    if _rc {
        di as error "engine not found: `engine'"
        di as error "check GARF_HOME / enginedir() points at the GARF-binary folder."
        exit 601
    }

    // -------- choose the Python interpreter ---------------------------------
    if "`python'" == "" {
        if "`c(os)'" == "Windows"  local venvpy "`enginedir'/../.venv/Scripts/python.exe"
        else                        local venvpy "`enginedir'/../.venv/bin/python"
        capture confirm file "`venvpy'"
        if !_rc local python "`venvpy'"
        else if "`c(os)'" == "Windows" local python "python"
        else                            local python "python3"
    }

    // -------- decide: load into memory? persist the CSV files? ---------------
    // Default: load the matched panel into memory, DO NOT keep any files.
    local doload = cond("`noload'" != "", 0, 1)
    local dopersist = 0
    if "`export'" != "" {
        local dopersist = 1
        local outdir "`export'"
    }
    else if "`outdir'" != "" {
        local dopersist = 1
    }
    else if "`store'" != "" {
        local dopersist = 1
        local outdir "`c(pwd)'/garf_output"
    }
    if `dopersist' == 0 & `doload' == 0 {
        // nothing would be kept -> save to garf_output so the run isn't wasted
        local dopersist = 1
        local outdir "`c(pwd)'/garf_output"
        di as txt "[garf] noload without store/export: outputs will be saved to garf_output/."
    }
    if `dopersist' == 0 {
        // temporary directory, deleted after the panel is loaded
        local outdir "`c(tmpdir)'/garf_`=runiformint(1,1000000000)'"
    }
    local outdir_fwd = subinstr("`outdir'", "\", "/", .)

    // -------- export the (sub)sample to a temporary CSV ---------------------
    marksample touse, novarlist
    tempfile datacsv
    local datacsv "`datacsv'.csv"
    preserve
        quietly keep if `touse'
        quietly keep `id' `time' `treatment' `outcome' `tvarying' `tinvariant'
        // nolabel: write the numeric codes of value-labelled variables, not their labels
        quietly export delimited using "`datacsv'", replace nolabel
    restore
    local datacsv_fwd = subinstr("`datacsv'", "\", "/", .)

    // -------- boolean flags -> yaml -----------------------------------------
    local repl  = cond("`noreplacement'" != "", "false", "true")
    local inct0 = cond("`includet0'"     != "", "true",  "false")
    local dotun = cond("`tune'"          != "", "true",  "false")
    local prog  = cond("`verbose'"       != "", "true",  "false")

    // -------- donor-selection cap: k-NN cap unless nokcap (caliper-only) -----
    if "`nokcap'" != "" {
        local kwrite "null"
        local kdisp  "all"
    }
    else {
        local kwrite "`k'"
        local kdisp  "`k'"
    }

    // -------- write the config file -----------------------------------------
    tempfile cfg
    local cfg "`cfg'.yaml"
    tempname fh
    file open `fh' using "`cfg'", write replace text
    file write `fh' "# auto-generated by garf.ado" _n
    // names and paths are written as double-quoted YAML strings, so that a variable
    // called on/off/yes/no/true/false/null, or a path containing " #" or ": ", is read
    // as text (paths use forward slashes, so no YAML escape is triggered)
    file write `fh' `"data_path: "`datacsv_fwd'""' _n
    file write `fh' `"id_col: "`id'""' _n
    file write `fh' `"year_col: "`time'""' _n
    file write `fh' `"treat_col: "`treatment'""' _n
    file write `fh' `"outcome_col: "`outcome'""' _n
    if "`tvarying'" == "" {
        file write `fh' "x_time_varying: []" _n
    }
    else {
        file write `fh' "x_time_varying:" _n
        foreach v of local tvarying {
            file write `fh' `"  - "`v'""' _n
        }
    }
    if "`tinvariant'" == "" {
        file write `fh' "x_time_invariant: []" _n
    }
    else {
        file write `fh' "x_time_invariant:" _n
        foreach v of local tinvariant {
            file write `fh' `"  - "`v'""' _n
        }
    }
    file write `fh' "treatment_year: `t0'" _n
    file write `fh' "n_pre_years: `npre'" _n
    file write `fh' "include_treatment_year: `inct0'" _n
    file write `fh' "k_matches: `kwrite'" _n
    file write `fh' "caliper: `caliper'" _n
    file write `fh' "with_replacement: `repl'" _n
    file write `fh' "weight_power: `wpower'" _n
    file write `fh' "tune: `dotun'" _n
    file write `fh' "tuning_metric: `tuningmetric'" _n
    file write `fh' "n_iter_search: `niter'" _n
    file write `fh' "cv_folds: `cvfolds'" _n
    file write `fh' "n_estimators: `nestimators'" _n
    file write `fh' "min_samples_leaf: `minleaf'" _n
    file write `fh' "max_features: `maxfeatures'" _n
    file write `fh' "max_depth: `maxdepth'" _n
    file write `fh' "random_state: `seed'" _n
    file write `fh' `"out_dir: "`outdir_fwd'""' _n
    file write `fh' "progress: `prog'" _n
    file close `fh'

    // -------- run the engine -------------------------------------------------
    if "`c(os)'" == "Windows" {
        local s_py  = subinstr(`"`python'"', "/", "\", .)
        local s_eng = subinstr("`engine'",   "/", "\", .)
        local s_cfg = subinstr("`cfg'",      "/", "\", .)
    }
    else {
        local s_py  `"`python'"'
        local s_eng "`engine'"
        local s_cfg "`cfg'"
    }
    di as txt "[garf] engine : `engine'"
    di as txt "[garf] python : `python'"
    di as txt "[garf] tuning : " cond("`dotun'"=="true","ON (slow)","off") ///
        "   k=`kdisp'   caliper=`caliper' (normalised)   wpower=`wpower'   replacement=`repl'"
    di as txt "[garf] running... (this can take a while when tuning is ON)"
    shell "`s_py'" "`s_eng'" --config "`s_cfg'"

    // -------- erase the temporary copy of the data and the config -----------
    // (Stata does not delete these itself: the .csv/.yaml suffix was added to the tempfile name)
    capture erase "`datacsv'"
    capture erase "`cfg'"

    // -------- check the engine produced results -----------------------------
    local mp "`outdir_fwd'/matched_panel.csv"
    capture confirm file "`mp'"
    if _rc {
        di as error "[garf] engine did not produce matched_panel.csv."
        di as error "Run the same call from a terminal to see the Python error, or check that"
        di as error "the Python interpreter above has pandas/scikit-learn installed."
        exit 601
    }

    // -------- parse the flat diagnostics ------------------------------------
    foreach s in n_treated n_controls n_treated_matched n_unique_donors_used ///
                 mean_abs_smd_raw mean_abs_smd_matched max_abs_smd_raw max_abs_smd_matched ///
                 mean_ESS_per_treated oob_score oob_balanced_accuracy oob_auc cv_score ///
                 fmq_mean fmq_p10 fmq_mean_all fmq_twin_mean fmq_undefined {
        local `s' = .
    }
    tempname jh
    capture file open `jh' using "`outdir_fwd'/diagnostics.json", read text
    if !_rc {
        file read `jh' line
        while r(eof) == 0 {
            foreach s in n_treated n_controls n_treated_matched n_unique_donors_used ///
                         mean_abs_smd_raw mean_abs_smd_matched max_abs_smd_raw max_abs_smd_matched ///
                         mean_ESS_per_treated oob_score oob_balanced_accuracy oob_auc cv_score ///
                         fmq_mean fmq_p10 fmq_mean_all fmq_twin_mean fmq_undefined {
                if regexm(`"`line'"', `""`s'": ([0-9eE.+-]+)"') local `s' = real(regexs(1))
            }
            file read `jh' line
        }
        file close `jh'
    }

    // -------- report ---------------------------------------------------------
    di as txt _n "{hline 62}"
    di as txt "GARF (binary) — matching summary"
    di as txt "{hline 62}"
    di as txt "Treated units              = " as res %9.0f `n_treated'
    di as txt "Control units              = " as res %9.0f `n_controls'
    di as txt "Treated matched            = " as res %9.0f `n_treated_matched'
    di as txt "Distinct donors used       = " as res %9.0f `n_unique_donors_used'
    di as txt "Mean |SMD| (raw -> matched)= " as res %5.3f `mean_abs_smd_raw' ///
        as txt "  ->  " as res %5.3f `mean_abs_smd_matched'
    di as txt "Max  |SMD| (raw -> matched)= " as res %5.3f `max_abs_smd_raw' ///
        as txt "  ->  " as res %5.3f `max_abs_smd_matched'
    di as txt "Mean ESS per treated       = " as res %5.1f `mean_ESS_per_treated' as txt " / `kdisp'"
    di as txt "Match quality FMQ mean/p10 = " as res %5.3f `fmq_mean' as txt " / " as res %5.3f `fmq_p10' ///
        as txt "   (all treated: " as res %5.3f `fmq_mean_all' as txt ")"
    di as txt "  FMQ: 0 = donors no closer than random controls, 1 = as close as the best attainable donor."
    di as txt "  Fewer/closer donors systematically raise it: compare designs together with ESS."
    if `fmq_undefined' > 0 & `fmq_undefined' < . ///
        di as txt "  FMQ undefined (no control closer than average) for " as res `fmq_undefined' as txt " block(s)."
    di as txt "Forest OOB bal.acc / AUC   = " as res %5.3f `oob_balanced_accuracy' ///
        as txt " / " as res %5.3f `oob_auc' as txt "   (OOB accuracy " as res %5.3f `oob_score' as txt ")"
    di as txt "  Near 0.5 = treatment hard to predict (good overlap); near 1 = groups separable."
    di as txt "  With rare treatment read the AUC first: balanced accuracy can stay near 0.5."
    if "`dotun'" == "true" ///
        di as txt "CV score (`tuningmetric')       = " as res %5.3f `cv_score'
    di as txt "{hline 62}"

    // -------- returns --------------------------------------------------------
    return scalar n_treated           = `n_treated'
    return scalar n_controls          = `n_controls'
    return scalar n_treated_matched   = `n_treated_matched'
    return scalar n_donors            = `n_unique_donors_used'
    return scalar mean_smd_raw        = `mean_abs_smd_raw'
    return scalar mean_smd_matched    = `mean_abs_smd_matched'
    return scalar max_smd_raw         = `max_abs_smd_raw'
    return scalar max_smd_matched     = `max_abs_smd_matched'
    return scalar mean_ess            = `mean_ESS_per_treated'
    return scalar fmq_mean            = `fmq_mean'
    return scalar fmq_p10             = `fmq_p10'
    return scalar fmq_mean_all        = `fmq_mean_all'
    return scalar fmq_twin_mean       = `fmq_twin_mean'
    return scalar fmq_undefined       = `fmq_undefined'
    return scalar oob_score           = `oob_score'
    return scalar oob_bacc            = `oob_balanced_accuracy'
    return scalar oob_auc             = `oob_auc'
    return scalar wpower              = `wpower'
    if "`dotun'" == "true" return scalar cv_score = `cv_score'

    // -------- load the matched panel into memory (default) ------------------
    if `doload' {
        di as txt "[garf] loading matched panel into memory (current data replaced)."
        import delimited using "`mp'", clear varnames(1) case(preserve)

        // numeric identifiers for fixed effects / clustering (string ids -> groups)
        capture confirm variable unique_id
        if !_rc {
            quietly egen long garf_unit = group(unique_id)
            label var garf_unit "GARF: numeric id of unique_id (donor-in-block); use in absorb()"
        }
        capture confirm variable block_id
        if !_rc {
            quietly egen long garf_block = group(block_id)
            label var garf_block "GARF: numeric id of block_id (the treated unit); use in cluster()"
        }
        // group / control dummies from role
        capture confirm string variable role
        if !_rc {
            quietly gen byte treated_grp = (role == "treated")
            quietly gen byte control_sel = (role == "donor")
            label var treated_grp "GARF: 1=treated unit, 0=selected control (donor)"
            label var control_sel "GARF: 1=untreated unit selected as a control (donor)"
        }
        // post indicator
        capture confirm variable `time'
        if !_rc {
            quietly gen byte post = (`time' >= `t0')
            label var post "GARF: 1 if period >= T0 (post-treatment)"
        }

        di as txt "[garf] matched panel is now in memory. New variables:"
        di as txt "       weight            design weight ~ proximity^`wpower', sums to 1 per block (use as [pw=weight])"
        di as txt "       proximity_minmax  match quality, 0-1 (1 = closest pair in this run)"
        di as txt "       cooccurrence  share of trees in which donor and treated unit share a leaf"
        di as txt "       fmq           forest match quality of the block (0 = random controls, 1 = best attainable donors)"
        di as txt "       treated_id    the treated unit each control is matched to (= block_id)"
        di as txt "       treated_grp   1=treated unit, 0=selected control"
        di as txt "       control_sel   1=this row is a selected control (donor)"
        di as txt "       garf_unit garf_block   numeric ids for absorb()/cluster()"
        di as txt "       post          1 if `time' >= `t0'"
        di as txt "[garf] estimate the matching DiD, e.g.:"
        di as txt `"       reghdfe `outcome' `treatment' [pw=weight], absorb(garf_unit `time') vce(cluster garf_block)"'
    }

    // -------- keep or clean up the CSV files --------------------------------
    if `dopersist' {
        di as txt "[garf] output files saved in: `outdir_fwd'"
        return local outdir        "`outdir_fwd'"
        return local matched_panel "`mp'"
    }
    else {
        foreach f in matched_panel.csv matches.csv balance.csv blocks.csv features.csv diagnostics.json {
            capture erase "`outdir_fwd'/`f'"
        }
        capture rmdir "`outdir_fwd'"
    }
end
