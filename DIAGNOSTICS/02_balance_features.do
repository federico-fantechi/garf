*! 02_balance_features.do  --  feature-level balance (paper: 5.3, Table 7)
*  ---------------------------------------------------------------------------
*  Summarises balance.csv by variable and by component (levels vs dynamics):
*  mean and max |SMD| before and after matching, number of features with
*  matched |SMD| > 0.10 and > 0.25, and number of features that got WORSE.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"
import delimited "garf_output/balance.csv", clear varnames(1)

* feature names are  variable__component  (e.g. Y__slope, x1__L3)
gen variable  = cond(strpos(feature,"__")>0, substr(feature,1,strpos(feature,"__")-1), feature)
gen component = substr(feature, strpos(feature,"__")+2, .) if strpos(feature,"__")

* dynamics = trajectory descriptors; levels = lags and the window mean
gen byte dynamic = inlist(component,"delta","growth","slope","vol","curv","dip")

gen byte above10 = abs_matched > 0.10    // conventional threshold
gen byte above25 = abs_matched > 0.25    // Stuart's upper bound for regression adjustment
gen byte worse   = abs_matched > abs_raw // balance worsened by matching

collapse (count) n=abs_raw (mean) abs_raw abs_matched                           ///
         (max) max_raw=abs_raw max_matched=abs_matched (sum) above10 above25 worse, ///
         by(variable dynamic)
label define dyn 0 "levels" 1 "dynamics", replace
label values dynamic dyn
list, sepby(variable) noobs

* full feature-by-feature table, worst matched imbalance first:
* import delimited "garf_output/balance.csv", clear varnames(1)
* gsort -abs_matched
* list feature abs_raw abs_matched, noobs
