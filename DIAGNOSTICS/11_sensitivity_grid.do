*! 11_sensitivity_grid.do  --  sensitivity of the design (paper: 5.11, Table 12)
*  ---------------------------------------------------------------------------
*  Runs garf on 19 designs, each changing ONE choice relative to the default
*  (k=50, no caliper, gamma=1, with replacement, seed 42), and stores the
*  diagnostics in sens.dta (used by 12_frontier.do):
*    cap k: 1 5 10 25 50 100 and the widest design $KWIDE (default k(1000):
*    the forest kernel on pools below 1,000 untreated units; see 00_settings.do)
*    weight exponent: gamma = 0, 2
*    caliper on normalised proximity: 0.1 0.2 0.3 (with k=50)
*    without replacement: k = 1, 3
*    seeds 1-5 (algorithmic variability)
*  Choose the grid BEFORE looking at any outcome. Slow: 19 forest fits.
*  Each run writes to runs/dI with export(): without it, noload would write
*  to garf_output/ and overwrite the files of 01_first_run.do.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
use raw.dta, clear

local d1  "k(1)"
local d2  "k(5)"
local d3  "k(10)"
local d4  "k(25)"
local d5  "k(50)"
local d6  "k(100)"
local d7  "$KWIDE"
local d8  "k(50) wpower(0)"
local d9  "k(50) wpower(2)"
local d10 "k(50) caliper(0.1)"
local d11 "k(50) caliper(0.2)"
local d12 "k(50) caliper(0.3)"
local d13 "k(1) noreplacement"
local d14 "k(3) noreplacement"
forvalues s = 1/5 {
    local d`=14+`s'' "k(50) seed(`s')"
}
local ndes 19

capture mkdir runs
tempname S
postfile `S' str40 design smd_mean smd_max ess fmq fmq_all fmq_twin ndonors nmatched ///
    using sens.dta, replace
forvalues i = 1/`ndes' {
    di as txt _n "[grid] design `i'/`ndes': `d`i''"
    garf, $GARF_CALL `d`i'' noload export("runs/d`i'")
    post `S' ("`d`i''") (r(mean_smd_matched)) (r(max_smd_matched)) (r(mean_ess)) ///
             (r(fmq_mean)) (r(fmq_mean_all)) (r(fmq_twin_mean)) (r(n_donors)) (r(n_treated_matched))
}
postclose `S'

use sens.dta, clear
format smd_mean smd_max fmq fmq_all fmq_twin %5.3f
format ess %6.2f
list, noobs sep(0)
