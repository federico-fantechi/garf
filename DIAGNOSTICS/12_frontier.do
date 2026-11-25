*! 12_frontier.do  --  the FMQ-ESS frontier (paper: 5.11, Figure 6)
*  ---------------------------------------------------------------------------
*  One point per design of 11_sensitivity_grid.do: FMQ (all-treated mean, so
*  that dropping hard-to-match units under a caliper is not rewarded) against
*  mean ESS (log scale). The connected points vary the cap k; designs above
*  that line give more match quality for the same effective sample.
*  Requires sens.dta. Label positions may need adjusting on other data
*  (variable pos: 12 = above, 6 = below, 3 = right, 9 = left).
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
use sens.dta, clear

gen fmq_plot = fmq_all
gen str12 grp = "k"
replace grp = "gamma" if strpos(design,"wpower")
replace grp = "cal"   if strpos(design,"caliper")
replace grp = "norep" if strpos(design,"noreplacement")
replace grp = "seed"  if strpos(design,"seed")
replace grp = "yours" if strpos(design,"yours:") == 1

gen str16 lab = ""
replace lab = "k=1"    if design=="k(1)"
replace lab = "k=5"    if design=="k(5)"
replace lab = "k=10"   if design=="k(10)"
replace lab = "k=50"   if design=="k(50)"
replace lab = "kernel" if design=="nokcap"
replace lab = subinstr(subinstr(design,"k(","k=",1),")","",1) if design=="$KWIDE" & design!="nokcap"
replace lab = "{&gamma}=0" if design=="k(50) wpower(0)"
replace lab = "{&gamma}=2" if design=="k(50) wpower(2)"
replace lab = "c=0.1"  if design=="k(50) caliper(0.1)"
replace lab = "c=0.2"  if design=="k(50) caliper(0.2)"
replace lab = "c=0.3"  if design=="k(50) caliper(0.3)"
replace lab = "no repl., k=1" if design=="k(1) noreplacement"
replace lab = "no repl., k=3" if design=="k(3) noreplacement"
replace lab = "your design"   if grp=="yours"

gen byte pos = 12
replace pos = 6 if inlist(design,"k(5)","k(50)","k(3) noreplacement")
replace pos = 3 if inlist(design,"k(50) wpower(2)","k(50) caliper(0.2)")
replace pos = 9 if inlist(design,"k(50) caliper(0.3)","k(50) wpower(0)","k(10)")
replace pos = 12 if inlist(design,"nokcap","$KWIDE")
quietly summarize ess
local xmax = r(max)
local xl "1 2 5 10 20 50"
foreach t in 100 200 500 1000 2000 5000 {
    if `t' <= 1.2*`xmax' local xl "`xl' `t'"
}

* your design (after your own garf): a filled square
local yplot ""
local yleg  ""
quietly count if grp == "yours"
if r(N) {
    local yplot `"(scatter fmq_plot ess if grp=="yours", mcolor(black) msymbol(S) msize(medlarge) mlabel(lab) mlabvposition(pos) mlabsize(vsmall) mlabcolor(black))"'
    local yleg `"6 "Your design""'
}
sort grp ess
twoway (connected fmq_plot ess if grp=="k", lcolor(black) mcolor(black) msymbol(O) ///
            mlabel(lab) mlabvposition(pos) mlabsize(vsmall) mlabcolor(black))         ///
       (scatter fmq_plot ess if grp=="gamma", mcolor(black) msymbol(D)                 ///
            mlabel(lab) mlabvposition(pos) mlabsize(vsmall) mlabcolor(black))          ///
       (scatter fmq_plot ess if grp=="cal", mcolor(black) msymbol(T)                   ///
            mlabel(lab) mlabvposition(pos) mlabsize(vsmall) mlabcolor(black))          ///
       (scatter fmq_plot ess if grp=="norep", mcolor(black) msymbol(X) msize(large)    ///
            mlabel(lab) mlabvposition(pos) mlabsize(vsmall) mlabcolor(black))          ///
       (scatter fmq_plot ess if grp=="seed", mcolor(gs10) msymbol(Oh) msize(small))    ///
       `yplot',                                                                        ///
    xscale(log range(0.8 `=1.3*`xmax'')) xlabel(`xl', nogrid) ylabel(, nogrid angle(0)) ///
    xtitle("Mean ESS per treated unit (log scale)") ytitle("FMQ (all-treated mean)")   ///
    legend(order(1 "Cap k ({&gamma}=1)" 2 "Weight exponent" 3 "Caliper"                ///
                 4 "Without replacement" 5 "Seeds 1-5" `yleg') cols(1) ring(0) pos(1) ///
           size(small) region(lstyle(none)))                                           ///
    plotregion(style(none)) graphregion(color(white)) xsize(6) ysize(4)
graph export "fig_frontier.pdf", replace
