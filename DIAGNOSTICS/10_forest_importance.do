*! 10_forest_importance.do  --  forest fit and top features (paper: 5.10, Table 11)
*  ---------------------------------------------------------------------------
*  Forest fit (out-of-bag): balanced accuracy and AUC (0.5 = chance; close to 1
*  = treated and untreated nearly separable: an overlap warning), plain
*  accuracy (read against the majority-class share).
*  The fifteen most important features are read from diagnostics.json.
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
di as txt "OOB balanced accuracy, AUC and accuracy are in the [garf] summary of"
di as txt "01_first_run.do, and in r(oob_bacc), r(oob_auc), r(oob_score)."

* --- top features: read the "top_features" block of diagnostics.json --------
tempname fh
file open `fh' using "garf_output/diagnostics.json", read text
file read `fh' line
local on 0
di as txt _n "Top features (mean decrease in impurity):"
while r(eof) == 0 {
    local l = strtrim(`"`macval(line)'"')
    if strpos(`"`l'"', `""top_features""') {
        local on 1
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
