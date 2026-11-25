*! 06_placebo.do  --  out-of-window placebo (paper: 5.6, Figure 4)
*  ---------------------------------------------------------------------------
*  The design checked (your garf design, or $GARF_DESIGN) is backdated by H
*  periods (default 3; lowered automatically if the data do not reach back far
*  enough for a full window) and run on pre-treatment data only:
*    pseudo-onset  T0p = T0 - H
*    matching window  T0p-L ... T0p-1
*    pseudo-post      T0p   ... T0-1   (true effect = 0)
*  A gap between treated units and their counterfactuals that opens in the
*  pseudo-post window means the matched history does not predict how the
*  treated evolve without treatment (doubt on unconfoundedness).
*  The raw gap (treated - all untreated) is shown for reference.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
use raw.dta, clear

* --- how far back: the backdated window must fit in the data -----------------
quietly summarize $T if !missing($Y)
local Hmax = $T0 - $L + $INCT0 - r(min)  // largest backdating with a full window
local H = min($H, `Hmax')
if `H' < 1 + $INCT0 {
    di as txt "[06_placebo] skipped: the data start in " r(min) ", too late for a full matching"
    di as txt "             window before a backdated onset (needs " $L " periods before " ///
        $T0 - 1 - $INCT0 ")."
    exit
}
if `H' < $H di as txt "[06_placebo] backdating H lowered from $H to `H': the data start in " r(min) "."
local T0p = $T0 - `H'

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
garf if $T < $T0, `call' $GARF_DESIGN $GARF_FOREST    // your design, backdated

collapse (mean) $Y [aw=weight], by($T treated_grp)
reshape wide $Y, i($T) j(treated_grp)
gen gap_garf = ${Y}1 - ${Y}0             // should stay near its in-window level
merge 1:1 $T using `r', nogen
list $T gap_garf gap_raw, noobs

* --- figure ----------------------------------------------------------------
local wstart = `T0p' - $L + $INCT0 - 0.5 // start of the matching window
local pstart = `T0p' + $INCT0 - 0.5      // start of the pseudo-post window
local xw = (`wstart' + `pstart')/2
local xp = (`pstart' + $T0 - 0.5)/2
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
