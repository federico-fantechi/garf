*! 10_forest_importance.do  --  forest fit and top features (paper: 5.10, Table 11)
*  ---------------------------------------------------------------------------
*  Forest fit (out-of-bag): balanced accuracy and AUC (0.5 = chance; close to 1
*  = treated and untreated nearly separable: an overlap warning), plain
*  accuracy (read against the majority-class share).
*  The out-of-bag statistics and the fifteen most important features are read
*  from diagnostics.json.
*  Read them by group: which variables, and whether their LEVELS or their
*  DYNAMICS, drive selection.
*  ---------------------------------------------------------------------------

cd "$GARF_WORK"

* --- fit statistics (from the first run; re-run it to refresh r()) -----------
use raw.dta, clear
bysort $ID: egen byte _everT = max($W)
bysort $ID: keep if _n == 1
quietly summarize _everT
local share = 1 - r(mean)
di as txt "Majority-class share (untreated) = " as res %5.3f `share'
* --- forest fit and top features: read diagnostics.json ------------------------
tempname fh
file open `fh' using "$GARF_OUTDIR/diagnostics.json", read text
file read `fh' line
local on 0
while r(eof) == 0 {
    local l = strtrim(`"`macval(line)'"')
    local isoob = strpos(`"`l'"', `""oob_auc""') | strpos(`"`l'"', `""oob_balanced_accuracy""') ///
        | strpos(`"`l'"', `""oob_score""')
    if !`on' & `isoob' {
        local l = subinstr(`"`l'"', `"""', "", .)
        local l = subinstr(`"`l'"', ",", "", .)
        di as txt "Out-of-bag " as res "`l'"
    }
    else if strpos(`"`l'"', `""top_features""') {
        local on 1
        di as txt _n "Top features (mean decrease in impurity):"
    }
    else if `on' & substr(`"`l'"', 1, 1) == "}" {
        local on 0
    }
    else if `on' {
        local l = subinstr(`"`l'"', `"""', "", .)
        local l = subinstr(`"`l'"', ",", "", .)
        di as res "   `l'"
    }
    file read `fh' line
}
file close `fh'
