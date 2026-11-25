*! 07_ess_reuse.do  --  effective sample size and donor reuse (paper: 5.7, Table 9)
*  ---------------------------------------------------------------------------
*  Per treated unit: ESS_i = 1 / sum_j w_ij^2 and the number of donors.
*  Per donor: reuse count r_j (treated units served) and implied weight
*  omega_j = sum_i w_ij (sums to the number of treated units).
*  Pool concentration: effective number of donors (sum omega)^2 / sum omega^2
*  and the share of total weight carried by the top 1, 2 and 10 donors.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
import delimited "$GARF_OUTDIR/matches.csv", clear varnames(1)

* --- per treated unit ----------------------------------------------------------
preserve
    gen w2 = weight^2
    collapse (sum) w2 (count) ndonors=rank, by(treated_id)
    gen ess = 1/w2
    tabstat ess ndonors, stats(mean p10 p50 p90 min max) format(%6.2f)
restore

* --- per donor ----------------------------------------------------------------
collapse (count) reuse=rank (sum) omega=weight, by(donor_id)
tabstat reuse omega, stats(mean p10 p50 p90 min max) format(%6.2f)

quietly summarize omega
local tot = r(sum)
gen o2 = omega^2
quietly summarize o2
di as txt "Effective number of donors (pool) = " as res %5.1f `tot'^2/r(sum)

gsort -omega
gen cum = sum(omega)/`tot'
local n10 = min(10, _N)
di as txt "Share of total implied weight: top 1 = " as res %5.3f cum[1] ///
   as txt "   top 2 = " as res %5.3f cum[min(2,_N)]                       ///
   as txt "   top 10 = " as res %5.3f cum[`n10']
list donor_id reuse omega in 1/`n10', noobs   // the most influential donors
