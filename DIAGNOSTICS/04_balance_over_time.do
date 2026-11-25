*! 04_balance_over_time.do  --  balance over time (paper: 5.4, Figure 2)
*  ---------------------------------------------------------------------------
*  Standardised difference treated - comparison, period by period, for the
*  outcome and the time-varying covariates, raw and GARF-matched, on a FIXED
*  scale (pooled pre-treatment SD in the raw data), so raw and matched are
*  comparable. Pre-treatment periods only: post-treatment differences belong
*  to the analysis stage. The vertical line marks the start of the window;
*  earlier periods were not used for matching (an out-of-sample check).
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
local X  "$Y $XTV"
local X : list uniq X                    // the outcome may also be among the covariates
local nv : word count `X'
local stubs ""
forvalues i = 1/`nv' {
    local stubs "`stubs' _v`i'_"
}
* variables are renamed _v1_ ... _vN_ inside this file, so that long variable
* names (up to 32 characters) do not overflow the generated names

* --- raw: treated vs all untreated -------------------------------------------
use raw.dta, clear
bysort $ID: egen byte _everT = max($W)  // ever-treated indicator
keep if $T <= $WLAST                     // design stage: pre-treatment only
keep $T _everT `X'
forvalues i = 1/`nv' {
    local v : word `i' of `X'
    rename `v' _v`i'_
    quietly summarize _v`i'_ if _everT==1 // fixed scale per variable
    local vt = r(Var)
    quietly summarize _v`i'_ if _everT==0
    local s`i' = sqrt((`vt' + r(Var))/2)
}
collapse (mean) `stubs', by($T _everT)
reshape wide `stubs', i($T) j(_everT)
forvalues i = 1/`nv' {
    gen smd_raw_`i' = (_v`i'_1 - _v`i'_0)/`s`i''
}
keep $T smd_raw_*
tempfile rawsmd
save `rawsmd'

* --- matched: treated vs donors weighted by w_ij -----------------------------
use matched.dta, clear
keep if $T <= $WLAST
keep $T treated_grp weight `X'
forvalues i = 1/`nv' {
    local v : word `i' of `X'
    rename `v' _v`i'_
}
collapse (mean) `stubs' [aw=weight], by($T treated_grp)
reshape wide `stubs', i($T) j(treated_grp)
forvalues i = 1/`nv' {
    gen smd_mat_`i' = (_v`i'_1 - _v`i'_0)/`s`i''
}
keep $T smd_mat_*
merge 1:1 $T using `rawsmd', nogen
sort $T
forvalues i = 1/`nv' {
    local v : word `i' of `X'
    label variable smd_raw_`i' "`v' (raw)"
    label variable smd_mat_`i' "`v' (matched)"
}
list, noobs                              // the numbers behind the figure (see labels)

* --- figure: one panel per variable, common y scale ---------------------------
local wstart = $WFIRST - 0.5             // start of the matching window
local xtxt   = ($WFIRST + $WLAST)/2      // label position (middle of the window)
local ytxt   = 0
forvalues i = 1/`nv' {                   // place the label above all lines
    quietly summarize smd_raw_`i'
    local ytxt = max(`ytxt', r(max))
    quietly summarize smd_mat_`i'
    local ytxt = max(`ytxt', r(max))
}
local ytxt = `ytxt' + 0.3
local glist ""
forvalues i = 1/`nv' {
    local v : word `i' of `X'
    local ttl "`v'"
    if "`v'" == "$Y" local ttl "Outcome `v'"
    twoway (line smd_raw_`i' $T, lcolor(black) lpattern(solid) lwidth(medthick)) ///
           (line smd_mat_`i' $T, lcolor(black) lpattern(dash)  lwidth(medthick)), ///
        yline(-0.1 0.1, lcolor(gs10) lpattern(shortdash)) yline(0, lcolor(gs13))  ///
        xline(`wstart', lcolor(gs6) lpattern(solid))                              ///
        text(`ytxt' `xtxt' "matching window", size(small) color(gs6))             ///
        xlabel(#6, nogrid) ylabel(, nogrid angle(0))                              ///
        title("`ttl'", size(medium) pos(11)) ytitle("SMD") xtitle("")            ///
        plotregion(style(none)) graphregion(color(white)) legend(off)             ///
        name(g_`i', replace) nodraw
    local glist "`glist' g_`i'"
}
local nr = ceil(`nv'/4)                  // at most four panels per row
local ys = 3.2*`nr'
graph combine `glist', rows(`nr') ycommon xsize(9) ysize(`ys') graphregion(color(white)) ///
    note("Solid: raw (treated vs all untreated). Dashed: GARF-matched (treated vs weighted donors). Dotted: {&plusminus}0.1.", size(small))
graph export "fig_balance_time.pdf", replace
