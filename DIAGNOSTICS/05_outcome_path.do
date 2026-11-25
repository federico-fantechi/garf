*! 05_outcome_path.do  --  pre-treatment outcome path (paper: 5.5, Figure 3)
*  ---------------------------------------------------------------------------
*  Mean outcome of the treated units, of their GARF counterfactuals (donors
*  weighted by w_ij) and of all untreated units, by period. Pre-treatment
*  only: plotting post-treatment periods would show the effect at the design
*  stage. With good overlap the treated and counterfactual paths nearly
*  coincide inside the window; a persistent gap shows where the donor pool
*  runs out of comparable histories.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"

use matched.dta, clear
collapse (mean) $Y [aw=weight], by($T treated_grp)
gen sample = "matched"
tempfile m
save `m'

use raw.dta, clear
bysort $ID: egen byte treated_grp = max($W)
collapse (mean) $Y, by($T treated_grp)
gen sample = "raw"
append using `m'
keep if $T <= $WLAST                     // post-treatment belongs to the analysis

list $T $Y treated_grp sample, sepby(sample) noobs

local wstart = $WFIRST - 0.5
local xtxt   = ($WFIRST + $WLAST)/2
quietly summarize $Y
local ytxt = r(max) + 0.05*(r(max) - r(min))

twoway (line $Y $T if treated_grp==1 & sample=="raw",     lcolor(black) lwidth(medthick))   ///
       (line $Y $T if treated_grp==0 & sample=="matched", lcolor(black) lpattern(dash)      ///
            lwidth(medthick))                                                               ///
       (line $Y $T if treated_grp==0 & sample=="raw",     lcolor(gs9) lwidth(medthick)),    ///
    xline(`wstart', lcolor(gs6)) text(`ytxt' `xtxt' "matching window", size(small) color(gs6)) ///
    xlabel(#10, nogrid) ylabel(, nogrid angle(0)) ytitle("Mean outcome") xtitle("")         ///
    legend(order(1 "Treated" 2 "GARF counterfactual" 3 "All untreated") cols(1) ring(0)    ///
           pos(7) size(small) region(lstyle(none)))                                        ///
    plotregion(style(none)) graphregion(color(white)) xsize(6) ysize(3.6)
graph export "fig_outcome_path.pdf", replace
