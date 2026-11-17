*! 03_variance_ratios.do  --  variance ratios (paper: 5.3, Table 8)
*  ---------------------------------------------------------------------------
*  Variance of each variable among treated units over the matching window
*  (unit-periods pooled), divided by its variance among
*    raw     : all untreated units
*    matched : the donors, weighted by the design weights w_ij
*  Usual tolerance: ratios between 0.5 and 2.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
local X  "$Y $XTV"
local nv : word count `X'
local lo = $T0 - $L                      // first period of the window
local hi = $T0 - 1                       // last period of the window

tempname R
matrix `R' = J(`nv', 2, .)
local r 0
foreach v of local X {
    local ++r
    * raw: treated vs all untreated
    use raw.dta, clear
    bysort $ID: egen byte _everT = max($W)
    quietly summarize `v' if _everT==1 & inrange($T, `lo', `hi')
    local vt = r(Var)
    quietly summarize `v' if _everT==0 & inrange($T, `lo', `hi')
    matrix `R'[`r',1] = `vt'/r(Var)
    * matched: treated vs donors weighted by w_ij
    use matched.dta, clear
    quietly summarize `v' if treated_grp==1 & inrange($T, `lo', `hi')
    local vt = r(Var)
    quietly summarize `v' [aw=weight] if treated_grp==0 & inrange($T, `lo', `hi')
    matrix `R'[`r',2] = `vt'/r(Var)
}
matrix rownames `R' = `X'
matrix colnames `R' = raw matched
di as txt _n "Variance ratios (treated / comparison), window `lo'-`hi'"
matrix list `R', format(%5.2f) noheader
