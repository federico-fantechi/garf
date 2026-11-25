{smcl}
{* garfte.sthlp — help for garfte, the GARF (binary) estimation stage — version 1.2}{...}
{title:Title}

{p 4 8}{bf:garfte} {hline 2} Treatment effects on a GARF matched panel: matching DiD and gap{p_end}
{p 4 8}{bf:version 1.2} (part of the GARF-binary package; see {help garf}){p_end}
{p 4 8}Federico Fantechi, PhD {hline 1} University of Palermo {hline 1} federico.fantechi@unipa.it{p_end}


{title:Syntax}

{p 4 8}On the matched panel that {cmd:garf} left in memory:{p_end}
{p 8 15}{bf:garfte}{cmd:,} {opt m:ethod(did|gap)} [{it:options}]{p_end}

{p 4 8}On a saved {cmd:garf} output folder (runs with {opt store}, {opt export()}, {opt outdir()}
or {opt noload}):{p_end}
{p 8 15}{bf:garfte using} {it:folder}{cmd:,} {opt m:ethod(did|gap)} [{it:options}]{p_end}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt m:ethod(did|gap)}}estimator; there is no default (the two rest on different assumptions){p_end}

{syntab:Estimation}
{synopt:{opt o:utcome(name)}}outcome to estimate on; default: the outcome {cmd:garf} matched on{p_end}
{synopt:{opt pre(#)}}keep event times e >= -#; default: all periods in the data{p_end}
{synopt:{opt post(#)}}keep event times e <= #; default: all periods in the data{p_end}
{synopt:{opt l:evel(#)}}confidence level, in per cent, for the table, the graph's confidence band and
the files; default {cmd:level(95)} (not affected by {cmd:set level}){p_end}

{syntab:Output}
{synopt:{opt nog:raph}}do not draw the event-study graph{p_end}
{synopt:{opt out:dir(str)}}folder for the result files; default {cmd:garf_output/estimations}
(see {help garfte##files:Files}){p_end}

{syntab:Plumbing (as in garf)}
{synopt:{opt engine:dir(str)}}path to the engine folder (else {cmd:$GARF_HOME}/engine){p_end}
{synopt:{opt py:thon(str)}}Python interpreter (else the package .venv, else system python){p_end}
{synoptline}
{p 4 8}{cmd:garfte} has no {cmd:if}/{cmd:in}: dropping rows of a matched panel would change the
counterfactual (see {help garfte##checks:Checks}).{p_end}


{title:Description}

{p 4 8}{cmd:garf} builds the counterfactual: for every treated unit {it:i} a block of donors {it:j}
with weights {it:w_ij} that sum to 1. {cmd:garfte} is the estimation stage that comes after it. It
runs the package's Python file {cmd:garf_estimate.py} (the same Python that runs {cmd:garf}),
prints a table of effects by event time and overall, draws an event-study graph, stores
{cmd:e(b)} and {cmd:e(V)}, and writes the effect of every treated unit in every period to a file.
{bf:The data in memory are never changed.}{p_end}

{p 4 8}Everything starts from the {bf:gap} between a treated unit and its own GARF counterfactual:{p_end}
{p 12 12}{cmd:G(i,t) = Y(i,t) - sum over j of w~(i,j,t) * Y(j,t)}{p_end}
{p 4 8}where {cmd:w~} are the {cmd:garf} weights, renormalised over the donors observed in period
{it:t} (they equal the {cmd:garf} weights when no donor is missing; see
{help garfte##missing:Missing donor-periods}). The two methods use the gap in two ways.{p_end}


{marker did}{...}
{title:method(did): matching difference-in-differences}

{p 12 12}{cmd:tau(i,e) = G(i, T0+e) - Gbar(i)},  {cmd:Gbar(i)} = mean of {cmd:G(i,t)} over the matching window{p_end}

{p 4 8}The window is the one {cmd:garf} matched on: [T0-npre, T0-1], or [T0-npre+1, T0] with
{opt includet0} (then T0 is a pre-treatment period and the effects start at T0+1). Measuring the
gap from its window mean removes any level difference left after matching. Identifying
assumption: {bf:parallel trends conditional on the match} (treated units and their donors would have
moved in parallel without treatment).{p_end}

{p 4 8}This is the estimand of Callaway and Sant'Anna's ATT(g,t) with a single cohort g = T0, and of
trajectory balancing with demeaning (Hazlett and Xu), with the GARF weights in place of
inverse-probability or balancing weights. Using the window mean rather than the last
pre-treatment period as the base averages out period-specific noise (as in synthetic DiD).{p_end}

{p 4 8}{bf:Requires a panel matched without replacement.} With replacement, a donor that serves several
treated units makes their blocks dependent and the block-clustered standard errors invalid.
{cmd:garfte} checks the panel itself: if any donor appears in more than one block,
{cmd:method(did)} stops and says so. Rerun {cmd:garf} with {opt noreplacement} and a {opt k()}
the donor pool can supply (about k x the number of treated units controls are needed; otherwise
the treated units matched last get fewer or worse donors). Standard errors are then clustered by
block, which is valid after matching without replacement (Abadie and Spiess 2022).{p_end}


{marker gap}{...}
{title:method(gap): treated minus its own donors, in levels}

{p 12 12}{cmd:tau(i,e) = G(i, T0+e)}{p_end}

{p 4 8}The gap itself, with no differencing: each treated unit against its own weighted donors, period
by period, as in the post-treatment gap of a synthetic control. Identifying assumption:
{bf:unconfoundedness given the matched pre-treatment history} (once the history is matched, the donors'
outcome is what the treated unit's would have been). Replacement is allowed: each treated unit
keeps its best donors even when another treated unit uses them too.{p_end}

{p 4 8}Gap and did answer the same question (the effect on the treated) under different
assumptions. Under fairly general conditions the two estimates bracket the true effect (Ding and
Li 2019): {bf:report both}. The gap is credible only if the window gaps are small: check the window
rows of the table and {cmd:e(mean_window_gap)}; a large level difference in the window means the
level comparison is not credible, and the did is the one to read.{p_end}


{title:Event times and the overall ATT}

{p 4 8}Event time e = period - T0 (e = 0 is T0). For each e:{p_end}
{p 12 12}{cmd:ATT(e) = mean of tau(i,e)} over the treated units with a defined gap at e{p_end}
{p 12 12}{cmd:ATT    = average of ATT(e)} over the post-treatment event times shown, each with the same
weight (Callaway and Sant'Anna's dynamic aggregation){p_end}
{p 4 8}With {opt pre()}/{opt post()} the table, the graph, {cmd:e(b)} and the overall ATT use only the
event times kept. The table labels each period:{p_end}
{p 8 10}{bf:pre}           before the matching window (if the data go back that far){p_end}
{p 8 10}{bf:window}        inside the matching window, for the outcome {cmd:garf} matched on{p_end}
{p 8 10}{bf:pre (window)}  inside the window, for an outcome {cmd:garf} did not match on{p_end}
{p 8 10}{bf:post}          after the window: the treatment effects{p_end}

{p 4 8}{bf:Window estimates are not a pre-trends test.} With {cmd:did} they sum to zero by construction;
with {cmd:gap} the outcome was matched on in the window, so they are small by construction. Evidence
on pre-trends comes from periods {it:before} the window ({bf:pre}), from an outcome that was not
matched on, and from the backdated placebo in {cmd:DIAGNOSTICS/06_placebo.do}.{p_end}


{marker se}{...}
{title:Standard errors}

{p 4 8}Inference is conditional on the matched design: the forest and the weights are taken as given (as
in propensity-score matching, first-stage uncertainty is not propagated). Each treated unit's
effect is a contrast between that unit and its block, so blocks are the clusters:{p_end}
{p 12 12}{cmd:V = G/(G-1) * sum over blocks of d(i) d(i)'},  {cmd:d(i,e) = (tau(i,e) - ATT(e)) / N(e)}{p_end}
{p 4 8}with G the number of treated units and N(e) the number with a gap at e: the block-clustered
variance. It equals the cluster-robust variance of the stacked regression shown in
{help garfte##replicate:Replicating garfte}, with the G/(G-1) correction.{p_end}

{p 4 8}{bf:Donors shared between blocks} ({cmd:gap} on a panel matched with replacement) make blocks
dependent: a shock to a donor used by ten treated units moves ten gaps. {cmd:garfte} then adds a
correction for each shared donor {it:j},{p_end}
{p 12 12}{cmd:sum over j of [A(j,e) A(j,f) - sum over i of a(i,j,e) a(i,j,f)] * S(j,e,f)},  {cmd:a(i,j,e) = w~(i,j) / N(e)},  {cmd:A(j,e) = sum over i of a(i,j,e)}{p_end}
{p 4 8}the covariance that donor {it:j} creates across the blocks that share it. {cmd:S(j)} is the
variance-covariance (across periods) of donor outcomes around the counterfactual, estimated inside the
blocks that use {it:j} (a local, matching-based variance, as in Abadie and Imbens). Without shared
donors the correction is zero. In simulations with heavy reuse the uncorrected block-clustered
intervals under-cover (about 90% for a nominal 95%); with the correction coverage is close to 95%.{p_end}

{p 4 8}Critical values and p-values use Student's t with G - 1 degrees of freedom. {bf:With fewer than
about 30 treated units} cluster-type standard errors can be unreliable: {cmd:garfte} prints a warning,
and the confidence intervals should be read with care.{p_end}


{marker missing}{...}
{title:Missing donor-periods}

{p 4 8}If some donors of a treated unit have no outcome in a period (unbalanced panel, or an
{opt outcome()} with gaps), the weights are renormalised over the donors observed in that period, so
the counterfactual is still a weighted mean of observed donors. The number of donors observed in each
cell is in {cmd:n_donors_obs} ({it:units} file) and its mean in {cmd:mean_donors_obs} ({it:event}
file), and the number of renormalised cells is printed. A cell with no donor observed, or with the
treated unit's outcome missing, is skipped. With {cmd:did}, a treated unit with no gap in the
window has no base and is left out.{p_end}


{marker checks}{...}
{title:Checks}

{p 4 8}Before estimating, {cmd:garfte} verifies that the panel is a complete GARF matched panel:{p_end}
{p 8 10}- in memory: it carries the settings that {cmd:garf} (version 1.2 or newer) stores as
characteristics {cmd:_dta[garf_*]} (they survive {cmd:save}/{cmd:use}); otherwise: {it:no GARF
matched panel in memory};{p_end}
{p 8 10}- every block has one treated unit, and the donor weights of every block still sum to 1
(tolerance 1e-4); otherwise: {it:matched panel was modified} {hline 1} typically because donors were
dropped. Dropping whole periods is allowed (they are simply not estimated), and so is dropping whole
blocks (a treated unit with all its donors, e.g. to estimate a subgroup ATT); dropping donors is
not, because it changes the counterfactual;{p_end}
{p 8 10}- {cmd:method(did)}: no donor appears in more than one block.{p_end}


{marker files}{...}
{title:Files}

{p 4 8}Written to {cmd:garf_output/estimations/} in the current folder; if {cmd:garf} saved its
output to a folder, to {it:that folder}{cmd:/estimations/}; with {cmd:using}, to
{it:folder}{cmd:/estimations/}; or to {opt outdir()}. Earlier results of the same method in that
folder are replaced.{p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:_event.csv}   name, e, year, period, att, se, t, p, lb, ub, n_treated,
mean_donors_obs: one row per event time and a last row {cmd:ATT} (overall){p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:_units.csv}   the effect of {bf:each treated unit in each period}: treated_id,
block_id, year, e, period, tau, pre_gap_mean (mean window gap of the unit), n_donors_obs, n_donors, ess,
fmq{p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:_vcov.csv}    the full covariance matrix of the event-time ATTs and the ATT{p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:_summary.csv}  key numbers (key,value){p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:_messages.txt} the notes and warnings printed on screen{p_end}
{p 8 10}{cmd:garfte_}{it:method}{cmd:.gph} and {cmd:.pdf}  the event-study graph (not with {opt nograph}){p_end}

{p 4 8}The graph, in the style of synthetic-control output, draws ATT(e) as a line with its
{opt level()}% confidence band, shades the matching window, marks the start of the treatment with a
dashed vertical line and the zero effect with a horizontal one, and reports the overall ATT with its
interval in the subtitle; the legend sits below the plot. Colours follow the current graph
scheme ({cmd:set scheme}): the line and the band take the scheme's first plot colour, and in Stata 15
or newer the band and the window are drawn with transparency, so the line, the band and the window
stay visible where they overlap. Stata 14 has no transparency: there the graph uses fixed light
colours. It is saved as a .gph, so it can be edited further in Stata's graph editor or with
{cmd:graph use}.{p_end}


{marker het}{...}
{title:Heterogeneous effects: the units file}

{p 4 8}{cmd:garfte_}{it:method}{cmd:_units.csv} holds one effect per treated unit and period, so the effects can
be related to characteristics of the treated units. Two cautions:{p_end}
{p 8 10}1. {bf:Unit-level effects are noisy.} The donor average removes noise from the counterfactual, not
from the treated unit's own outcome, so the spread of {cmd:tau} across units overstates the spread of
the true effects. Regressing {cmd:tau} on characteristics is fine (noise in the dependent variable
does not bias the coefficients); reading the variance of {cmd:tau} as effect heterogeneity is not.{p_end}
{p 8 10}2. {bf:With shared donors the unit effects are correlated}, and plain OLS standard errors of such a
regression are too small. Without shared donors ({cmd:did}, or a panel matched without replacement),
the units are independent across treated units: cluster by {cmd:treated_id} when you use several
periods.{p_end}
{p 4 8}Subgroup ATTs can also be obtained, without any re-estimation of the matching, by keeping whole
blocks: tag each block by a pre-treatment characteristic of its treated unit (see {help garf}, Studying
heterogeneous effects), keep the blocks of one group, save that panel, and run {cmd:garfte} on it.{p_end}


{marker replicate}{...}
{title:Replicating garfte with a regression}

{p 4 8}The point estimates are those of a stacked regression with block-by-period fixed effects and the
GARF weights (treated rows weight 1, donor weights summing to 1, so every block counts equally):{p_end}
{phang2}{cmd:. reghdfe Y ibn.year#c.treated_grp [pw=weight], absorb(garf_block#year) vce(cluster garf_block)}{p_end}
{p 4 8}reproduces {cmd:method(gap)}; the same regression with the outcome measured from each unit's mean over
the matching window reproduces {cmd:method(did)}. With no shared donors the standard errors agree up to
reghdfe's small-sample factor. Note that {cmd:absorb(garf_unit year)}, the quick DiD shown in {help garf} for a
quick DiD, uses time effects common to all blocks; it is close to, but not the same as, the block-wise
comparison of {cmd:garfte}.{p_end}


{title:Examples}

{pstd}Gap on the default (with-replacement) match:{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) tinvariant(region)}{p_end}
{phang2}{cmd:. garfte, method(gap)}{p_end}

{pstd}Matching DiD (needs a match without replacement), only 5 years before and after T0:{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) noreplacement k(3)}{p_end}
{phang2}{cmd:. garfte, method(did) pre(5) post(5)}{p_end}
{phang2}{cmd:. garfte, method(did) level(90)}     // 90% confidence band and intervals{p_end}
{phang2}{cmd:. garfte, method(gap) pre(5) post(5)}     // the bracketing estimate{p_end}

{pstd}Another outcome, merged into the matched panel (match on Y, estimate on employment):{p_end}
{phang2}{cmd:. merge m:1 id year using employment.dta, keep(master match) nogenerate}{p_end}
{phang2}{cmd:. garfte, method(did) outcome(employment)}{p_end}

{pstd}On a saved output folder (the data in memory are not needed):{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) noload store}{p_end}
{phang2}{cmd:. garfte using "garf_output", method(gap)}{p_end}

{pstd}Coefficient plot and tests from the stored results:{p_end}
{phang2}{cmd:. coefplot, vertical keep(e*)}{p_end}
{phang2}{cmd:. test em10 em9 em8 em7 em6}     // periods before the window jointly zero{p_end}
{phang2}{cmd:. matrix list e(event)}{p_end}


{title:Stored results}

{p 4 8}{cmd:garfte} is an e-class command.{p_end}
{p 4 8}Matrices: {cmd:e(b)}, {cmd:e(V)} (ATT(e) named {cmd:em}{it:k} for e = -{it:k} and {cmd:e}{it:k} for
e = {it:k} >= 0, then {cmd:ATT}); {cmd:e(event)} (rows: event times and ATT; columns: e, year, att, se,
t, p, lb, ub, n_treated).{p_end}
{p 4 8}Scalars: {cmd:e(att)}, {cmd:e(se_att)}, {cmd:e(lb_att)}, {cmd:e(ub_att)}, {cmd:e(p_att)},
{cmd:e(level)}, {cmd:e(df_r)}, {cmd:e(N_treated)}, {cmd:e(N_donors)}, {cmd:e(donors_reused)},
{cmd:e(t0)}, {cmd:e(npre)}, {cmd:e(window_first)}, {cmd:e(window_last)}, {cmd:e(mean_window_gap)}.{p_end}
{p 4 8}Macros: {cmd:e(cmd)}, {cmd:e(method)}, {cmd:e(outcome)}, {cmd:e(matched_outcome)}, {cmd:e(outdir)},
{cmd:e(title)}, {cmd:e(properties)}.{p_end}


{title:From a terminal}

{phang2}{cmd:python garf_estimate.py --from garf_output --method did [--outcome Y --pre 5 --post 5 --level 95]}{p_end}
{p 4 8}runs the same estimation on a saved output folder (inside the package's {cmd:engine} folder, with
the package's Python). Results go to {cmd:garf_output/estimations}.{p_end}


{title:References}

{p 4 8}Abadie, A. and J. Spiess (2022). Robust post-matching inference. {it:Journal of the American
Statistical Association} 117: 983-995.{p_end}
{p 4 8}Callaway, B. and P. H. C. Sant'Anna (2021). Difference-in-differences with multiple time periods.
{it:Journal of Econometrics} 225: 200-230.{p_end}
{p 4 8}Ding, P. and F. Li (2019). A bracketing relationship between difference-in-differences and
lagged-dependent-variable adjustment. {it:Political Analysis} 27: 605-615.{p_end}
{p 4 8}Hazlett, C. and Y. Xu (2018). Trajectory balancing: a general reweighting approach to causal
inference with time-series cross-sectional data. Working paper, SSRN 3214231.{p_end}


{title:Also see}

{p 4 8}{help garf}{p_end}
