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

* --- raw: treated vs all untreated -------------------------------------------
use raw.dta, clear
bysort $ID: egen byte _everT = max($W)  // ever-treated indicator
foreach v of local X {                   // fixed scale per variable
    quietly summarize `v' if _everT==1 & $T < $T0
    local vt = r(Var)
    quietly summarize `v' if _everT==0 & $T < $T0
    local s_`v' = sqrt((`vt' + r(Var))/2)
}
collapse (mean) `X', by($T _everT)
reshape wide `X', i($T) j(_everT)
foreach v of local X {
    gen smd_raw_`v' = (`v'1 - `v'0)/`s_`v''
}
keep $T smd_raw_*
tempfile rawsmd
save `rawsmd'

* --- matched: treated vs donors weighted by w_ij -----------------------------
use matched.dta, clear
collapse (mean) `X' [aw=weight], by($T treated_grp)
reshape wide `X', i($T) j(treated_grp)
foreach v of local X {
    gen smd_mat_`v' = (`v'1 - `v'0)/`s_`v''
}
keep $T smd_mat_*
merge 1:1 $T using `rawsmd', nogen
keep if $T < $T0                         // design stage: pre-treatment only
list, noobs                              // the numbers behind the figure

* --- figure: one panel per variable, common y scale ---------------------------
local wstart = $T0 - $L - 0.5            // start of the matching window
local xtxt   = $T0 - $L/2 - 0.5          // label position (middle of the window)
local ytxt   = 0
foreach v of local X {                   // place the label above all lines
    quietly summarize smd_raw_`v'
    local ytxt = max(`ytxt', r(max))
    quietly summarize smd_mat_`v'
    local ytxt = max(`ytxt', r(max))
}
local ytxt = `ytxt' + 0.3
local glist ""
foreach v of local X {
    local ttl "`v'"
    if "`v'" == "$Y" local ttl "Outcome `v'"
    twoway (line smd_raw_`v' $T, lcolor(black) lpattern(solid) lwidth(medthick)) ///
           (line smd_mat_`v' $T, lcolor(black) lpattern(dash)  lwidth(medthick)), ///
        yline(-0.1 0.1, lcolor(gs10) lpattern(shortdash)) yline(0, lcolor(gs13))  ///
        xline(`wstart', lcolor(gs6) lpattern(solid))                              ///
        text(`ytxt' `xtxt' "matching window", size(small) color(gs6))             ///
        xlabel(#6, nogrid) ylabel(, nogrid angle(0))                              ///
        title("`ttl'", size(medium) pos(11)) ytitle("SMD") xtitle("")            ///
        plotregion(style(none)) graphregion(color(white)) legend(off)             ///
        name(g_`v', replace) nodraw
    local glist "`glist' g_`v'"
}
local nv : word count `X'
local nr = ceil(`nv'/4)                  // at most four panels per row
local ys = 3.2*`nr'
graph combine `glist', rows(`nr') ycommon xsize(9) ysize(`ys') graphregion(color(white)) ///
    note("Solid: raw (treated vs all untreated). Dashed: GARF-matched (treated vs weighted donors). Dotted: {&plusminus}0.1.", size(small))
graph export "fig_balance_time.pdf", replace
