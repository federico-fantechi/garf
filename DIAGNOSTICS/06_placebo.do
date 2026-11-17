*! 06_placebo.do  --  out-of-window placebo (paper: 5.6, Figure 4)
*  ---------------------------------------------------------------------------
*  The design is backdated by H periods and run on pre-treatment data only:
*    pseudo-onset  T0p = T0 - H
*    matching window  T0p-L ... T0p-1
*    pseudo-post      T0p   ... T0-1   (true effect = 0)
*  A gap between treated units and their counterfactuals that opens in the
*  pseudo-post window means the matched history does not predict how the
*  treated evolve without treatment (doubt on unconfoundedness).
*  The raw gap (treated - all untreated) is shown for reference.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
local T0p = $T0 - $H

use raw.dta, clear
* In the pre-treatment sample the treatment variable is 0 everywhere, so the
* ever-treated indicator defines who is treated.
bysort $ID: egen byte _everT = max($W)

* --- raw gap first: garf replaces the data in memory --------------------------
preserve
    keep if $T < $T0
    collapse (mean) $Y, by($T _everT)
    reshape wide $Y, i($T) j(_everT)
    gen gap_raw = ${Y}1 - ${Y}0
    keep $T gap_raw
    tempfile r
    save `r'
restore

* --- backdated design (default k, no files kept) ------------------------------
local call : subinstr global GARF_CALL "treatment($W)" "treatment(_everT)"
local call : subinstr local call "t0($T0)" "t0(`T0p')"
garf if $T < $T0, `call' k($K)

collapse (mean) $Y [aw=weight], by($T treated_grp)
reshape wide $Y, i($T) j(treated_grp)
gen gap_garf = ${Y}1 - ${Y}0             // should stay near its in-window level
merge 1:1 $T using `r', nogen
list $T gap_garf gap_raw, noobs

* --- figure ----------------------------------------------------------------
local wstart = `T0p' - $L - 0.5          // start of the matching window
local pstart = `T0p' - 0.5               // start of the pseudo-post window
local xw = `T0p' - $L/2 - 0.5
local xp = (`T0p' + $T0 - 1)/2
quietly summarize gap_raw
local ymax = max(r(max), 0)
local ymin = min(r(min), 0)
quietly summarize gap_garf
local ymax = max(`ymax', r(max))
local ymin = min(`ymin', r(min))
local ytxt = `ymax' + 0.12*(`ymax' - `ymin')      // labels above the lines, inside the axis
twoway (line gap_garf $T, lcolor(black) lpattern(dash)  lwidth(medthick))        ///
       (line gap_raw  $T, lcolor(gs9)   lpattern(solid) lwidth(medthick)),       ///
    yline(0, lcolor(gs10)) xline(`wstart' `pstart', lcolor(gs6) lpattern(solid)) ///
    text(`ytxt' `xw' "matching window", size(small) color(gs6))                  ///
    text(`ytxt' `xp' "pseudo-post", size(small) color(gs6))                      ///
    xscale(range(`wstart' .)) yscale(range(`ymin' `ytxt'))                       ///
    xlabel(#10, nogrid) ylabel(, nogrid angle(0))                                ///
    ytitle("Treated minus comparison") xtitle("")                                ///
    legend(order(1 "GARF counterfactual" 2 "All untreated") cols(1) ring(0) pos(7) ///
           size(small) region(lstyle(none)))                                     ///
    plotregion(style(none)) graphregion(color(white)) xsize(6) ysize(3.6)
graph export "fig_placebo.pdf", replace
