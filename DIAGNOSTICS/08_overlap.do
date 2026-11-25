*! 08_overlap.do  --  overlap (paper: 5.8, Figure 5)
*  ---------------------------------------------------------------------------
*  Two measures per treated unit, from blocks.csv:
*    best_proximity_minmax : normalised proximity of its best donor
*                            (0-1, on the run's common scale)
*    overlap_mass          : share of trees in which it shares a leaf with at
*                            least one untreated unit
*  Low values on either: few or no close untreated units. Unmatched treated
*  units (e.g. under a caliper) have matched == 0.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
import delimited "$GARF_OUTDIR/blocks.csv", clear varnames(1)

tabstat best_proximity_minmax overlap_mass, stats(mean p10 p50 p90 min max) format(%5.3f)

* the five treated units with the weakest best donor
sort best_proximity_minmax
list treated_id best_proximity_minmax overlap_mass n_donors fmq in 1/5, noobs

* treated units left without donors (none without a caliper)
count if matched == 0
if r(N) > 0 list treated_id overlap_mass if matched == 0, noobs

* --- figure: two histograms on [0,1] -----------------------------------------
histogram best_proximity_minmax, width(0.05) start(0) frequency fcolor(gs11) lcolor(gs6) ///
    xscale(range(0 1)) xlabel(0(0.2)1, nogrid) ylabel(, nogrid angle(0))             ///
    title("Best-donor normalised proximity", size(medium) pos(11))                    ///
    xtitle("") ytitle("Treated units") graphregion(color(white)) name(h1, replace) nodraw
histogram overlap_mass, width(0.05) start(0) frequency fcolor(gs11) lcolor(gs6)       ///
    xscale(range(0 1)) xlabel(0(0.2)1, nogrid) ylabel(, nogrid angle(0))             ///
    title("Overlap mass {it:m}{sub:i}", size(medium) pos(11))                         ///
    xtitle("") ytitle("Treated units") graphregion(color(white)) name(h2, replace) nodraw
graph combine h1 h2, rows(1) ycommon xsize(9) ysize(3.4) graphregion(color(white))
graph export "fig_overlap.pdf", replace
* with a caliper c, add  xline(c)  to the first histogram
