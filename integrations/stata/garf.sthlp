{smcl}
{* garf.sthlp — help for the GARF (binary) Stata wrapper — version 1.1}{...}
{title:Title}

{p 4 8}{bf:garf} {hline 2} Forest-proximity matching for panel counterfactuals (binary build){p_end}
{p 4 8}{bf:version 1.1}{p_end}
{p 4 8}Federico Fantechi, PhD {hline 1} University of Palermo {hline 1} federico.fantechi@unipa.it{p_end}


{title:What GARF does}

{p 4 8}{cmd:garf} runs the {it:Group Aware Random Forest} (GARF) engine from inside Stata: you give it
a panel and it returns a matched donor panel, including weights and treated-donor links, ready to be
used in any estimation.{p_end}

{p 4 8}GARF is a matching tool built to identify comparable (untreated) units over i) a very large set
of covariates and ii) a panel data structure, matching not only the {it:levels} of each covariate but
also their {it:dynamics}, thanks to the {it:group awareness} of the model: each variable's pre-treatment
history is encoded as a group of features (its level in every pre-treatment period plus mean, change,
growth, slope, volatility, curvature and pre-treatment dip). GARF is designed for the
settings where common strategies based on linear models (e.g. propensity score matching) often fail:
few treated units, panel data, many covariates, non-linear or threshold selection, selection on
pre-treatment trajectories. With few covariates and simple, roughly linear selection, GARF will perform
comparably to, or worse than, a well-specified propensity score.{p_end}
{p 4 8}GARF does {bf:not} estimate the treatment effect itself {hline 1} that stays in your hands.{p_end}
{p 4 8}This Stata command is a front end: it calls the packaged Python engine through Stata's {cmd:shell}
(Stata 14 or newer; Stata's own {cmd:python:} integration is {it:not} required). The forest and all
computation run in Python, using the package's own environment.{p_end}


{title:Syntax}

{p 8 15}
{bf:garf} {ifin}{cmd:,}
{opt tre:atment(varname)}
{opt o:utcome(varname)}
{opt id(varname)}
{opt ti:me(varname)}
{opt t0(#)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt tre:atment(varname)}}binary treatment column; a unit is treated if it is ever 1 (standard coding: 1 for treated units in treated years, else 0){p_end}
{synopt:{opt o:utcome(varname)}}outcome variable; enters the matching only through its pre-treatment window{p_end}
{synopt:{opt id(varname)}}unit identifier (numeric or string){p_end}
{synopt:{opt ti:me(varname)}}time variable{p_end}
{synopt:{opt t0(#)}}treatment onset period (T0){p_end}

{syntab:Design}
{synopt:{opt tv:arying(varlist)}}matching covariates that change over time (garf builds lags + trajectory features){p_end}
{synopt:{opt tin:variant(varlist)}}matching covariates that do not change (strings allowed){p_end}
{synopt:{opt np:re(#)}}pre-window length used for matching: [T0-npre, T0-1]; default 5 (minimum 2; curvature needs 4){p_end}
{synopt:{opt include:t0}}shift the pre-window to [T0-npre+1, T0], including the onset year{p_end}

{syntab:Matching (donor selection and weights)}
{synopt:{opt k(#)}}keep at most this many nearest donors per treated unit (k-NN); default 50{p_end}
{synopt:{opt cal:iper(#)}}keep only donors with {it:normalised} proximity >= # (range 0-1; 0 =
off, higher is stricter); COMBINES with k(){p_end}
{synopt:{opt nok:cap}}lift the k cap: keep ALL donors within the caliper{p_end}
{synopt:{opt nor:eplacement}}match without replacement, hardest-to-match treated units first (default is with replacement){p_end}
{synopt:{opt wp:ower(#)}}weights ~ proximity^#, summing to 1 per treated unit; default 1 (proportional); 0 = uniform k:1; >1 = more weight on the closest donors{p_end}

{syntab:Model selection (hyperparameter tuning)}
{synopt:{opt tune}}turn on cross-validated hyperparameter search (SLOW){p_end}
{synopt:{opt tuning:metric(str)}}accuracy | f1_treated | combined; default combined{p_end}
{synopt:{opt nit:er(#)}}search candidates; default 40{p_end}
{synopt:{opt cv:folds(#)}}cross-validation folds; default 5{p_end}

{syntab:Forest (used as-is when tune is off)}
{synopt:{opt n:estimators(#)}}trees; default 1000{p_end}
{synopt:{opt minl:eaf(#)}}min samples per leaf; default 5{p_end}
{synopt:{opt maxf:eatures(str)}}sqrt | log2 | a fraction; default sqrt{p_end}
{synopt:{opt maxd:epth(str)}}null or an integer; default null{p_end}
{synopt:{opt se:ed(#)}}random seed; default 42{p_end}

{syntab:Output (default: load the matched panel into memory, save no files)}
{synopt:{opt store}}also save the CSV outputs to ./garf_output{p_end}
{synopt:{opt export(str)}}also save the CSV outputs to this folder{p_end}
{synopt:{opt noload}}do NOT load the matched panel into memory (produce files only){p_end}

{syntab:Plumbing}
{synopt:{opt enginedir(str)}}path to the engine folder (else {cmd:$GARF_HOME}/engine){p_end}
{synopt:{opt python(str)}}Python interpreter to use (else the package .venv, else system python){p_end}
{synopt:{opt outdir(str)}}folder for saved files (implies saving); default ./garf_output{p_end}
{synopt:{opt verbose}}show the engine's own progress output{p_end}
{synoptline}


{title:Description}

{p 4 8}
{cmd:garf} is a thin wrapper around the packaged Python engine {cmd:garf_binary.py}. It exports
the current panel to a temporary CSV, writes a config from the options above, runs the engine,
and reports a matching summary. It does {bf:not} estimate the treatment effect: it produces a
matched donor panel that you feed to your own weighted DiD / event study.

{p 4 8}
{bf:By default the matched panel is loaded into memory and no files are saved.} This is the
right object for a matching DiD: with replacement, a control can be a donor to several treated
units with a different weight each time, so the result is a block-expanded panel (one block per
treated unit), not a few new columns on your original rows. Use {opt store} or {opt export()}
to also write the CSV outputs; use {opt noload} to skip loading. {bf:Loading replaces the data
in memory} — {cmd:preserve} or save your working data first.


{marker proximity}{...}
{title:How proximity is calculated}

{p 4 8}
GARF fits a random-forest classifier of the treatment {cmd:W} on the pre-treatment descriptor of
every unit, pushes every unit (treated and control) down every one of its {cmd:nestimators()}
trees, and reads off which terminal leaf each unit lands in. Two treated/control pairs that
repeatedly end up in the same (small) leaf are, under the forest's own learned metric,
"close" — this is Breiman's leaf-proximity idea, read as an adaptive nearest-neighbour kernel
(Lin and Jeon, 2006). Three related numbers come out of this, and it matters which one a given
part of the output is:{p_end}

{p 8 10}{bf:{cmd:proximity} (raw).} For treated unit {it:i} and control {it:j},{p_end}
{p 12 12}{cmd:proximity(i,j) = (1/B) * sum over trees b of  1{leaf_b(i)=leaf_b(j)} / (# controls sharing that leaf in tree b)}{p_end}
{p 8 10}
— the share of trees in which {it:i} and {it:j} co-occur, weighted so that sharing a small,
tight leaf counts for more than sharing a large, diffuse one. It is {bf:0} for a pair that never
once shares a leaf, and it climbs toward its (small, leaf-size-dependent) maximum as the pair
keeps landing together. {bf:Higher = closer/more similar} — the opposite convention from a
distance. Its {it:absolute} magnitude is not comparable across runs: it shrinks with
{opt minleaf()}, grows with {opt nestimators()}, and depends on how many covariates are in play,
so "proximity = 0.03" means something different in every specification.{p_end}

{p 8 10}{bf:{cmd:proximity_minmax} (normalised).} The raw proximity matrix,
rescaled {it:once, globally}, over every treated-control pair produced in that run:{p_end}
{p 12 12}{cmd:proximity_minmax(i,j) = (proximity(i,j) - min) / (max - min)}{p_end}
{p 8 10}
where {cmd:min}/{cmd:max} are the smallest/largest raw proximity value anywhere in the whole
treated x control matrix for that run (reported as {cmd:proximity_raw_min}/{cmd:proximity_raw_max}
in {cmd:diagnostics.json}). This is deliberately a {bf:global}, not a per-treated-unit,
rescaling: normalising each treated unit by its {it:own} best match would send every unit's
top donor to 1.0 regardless of how poor that donor actually is, which would erase the very
signal that tells you a treated unit has no good match anywhere in the data. With a global
rescaling, {cmd:proximity_minmax = 1} is reserved for the single closest pair observed in the
{it:whole} run, so a unit whose best available donor is mediocre still shows a mediocre score.
{bf:{opt caliper()} is applied on this scale} (range 0 to 1; 0 = off), so a caliper value
means the same thing regardless of forest size or the number of covariates — unlike a raw-scale
caliper.{p_end}

{p 8 10}{bf:{cmd:weight} (design weight).} Within each treated unit's own matched donor
set, the {it:raw} proximities, raised to the power set by {opt wpower()} (default 1), are
renormalised to sum to 1:{p_end}
{p 12 12}{cmd:weight(i,j) = proximity(i,j)^g / sum over the donors kept for i of proximity(i,*)^g},  g = {opt wpower()}{p_end}
{p 8 10}
With the default g = 1 the weight is proportional to proximity. With g = 0 every selected donor of a
treated unit gets the same weight (classic k:1 matching). With g > 1 weight is shifted toward the
closest donors (as g grows this tends to nearest-neighbour matching); this lowers the effective
sample size (ESS), so check ESS and balance when you raise it. This is what you pass to
{cmd:[pw=weight]} in the DiD step — a within-block convex combination, not a value on the 0-1 "how
good is this match" scale that {cmd:proximity_minmax} is.{p_end}


{title:Donor selection: k-NN and/or caliper}

{p 4 8}
Two rules pick each treated unit's controls, and they {bf:combine}:{p_end}
{p 8 10}{opt k(#)}       keep the # closest donors (nearest-neighbour). This is the default.{p_end}
{p 8 10}{opt caliper(#)} keep only donors with {it:normalised} proximity >= # (range 0-1).{p_end}
{p 4 8}
A donor must also have actually co-occurred with the treated unit in at least one tree (raw
{cmd:proximity > 0}); that floor applies regardless of the caliper, since a pair that never
shares a leaf carries no information and would receive weight 0 in any case.{p_end}
{p 4 8}
Examples of the two working together:{p_end}
{p 8 10}{cmd:k(50)} (default) — the 50 nearest donors, no caliper.{p_end}
{p 8 10}{cmd:caliper(0.3) nokcap} — every donor with normalised proximity >= 0.3, however many.{p_end}
{p 8 10}{cmd:caliper(0.3) k(50)} — the 50 nearest {it:among} those with normalised proximity >= 0.3.{p_end}
{p 4 8}
A caliper trades donors for quality: a tighter (larger) caliper means fewer, closer controls —
better balance but smaller effective samples, and possibly treated units left with {it:no}
donor at all (an honest signal of poor overlap, not a bug). Report how many donors survive
(see ESS below).


{title:Reading the summary}

{p 4 8}The on-screen summary and {cmd:r()} report:{p_end}
{p 8 10}{bf:Mean / Max |SMD| (raw -> matched)}: covariate balance. Matching should pull these
toward 0; |SMD| < 0.10 is the usual "well balanced" rule of thumb. The {it:raw} number is your
starting imbalance, the {it:matched} number what GARF achieved.{p_end}
{p 8 10}{bf:Match quality FMQ mean/p10 (all treated)}: Forest Match Quality, 0 = donors no closer than
random controls, 1 = every donor as close as the best donor attainable in the data; see
{help garf##fmq:Forest Match Quality} {hline 1} it systematically rewards fewer donors, so read it with ESS.{p_end}
{p 8 10}{bf:Mean ESS per treated}: the effective number of distinct donors behind each treated
unit. Well below {cmd:k} is normal (good donors get reused, or few pass the caliper) and worth
reporting; very low ESS warns the effective comparison group is thin.{p_end}
{p 8 10}{bf:Forest OOB bal.acc / AUC}: how separable treated and controls are, computed from the
forest's {it:out-of-bag} predictions (each unit is scored only by trees that did not see it). The
{bf:balanced accuracy} averages the hit rates of the two groups, so 0.5 means no better than chance
whatever the share of treated units; the {bf:AUC} is the probability that a random treated unit gets
a higher out-of-bag score than a random control. With rare treatment (a few per cent of units)
read the {bf:AUC first}: out-of-bag votes then seldom pass the 0.5 threshold, so the balanced accuracy can
stay close to 0.5 even when the AUC shows that the groups are highly separable. Values near 0.5 mean treatment is hard to predict
from the history (good for overlap); values close to 1 (say > 0.95) are an {it:overlap} warning
(treated and controls are almost perfectly distinguishable). The plain out-of-bag {it:accuracy}
({cmd:r(oob_score)}) is shown in brackets for reference: read it against the share of control
units, which a forest reaches by predicting "control" for everyone. Near-perfect separability is
the same lack-of-positivity problem that afflicts any matching method. If you see it, check that balance still improves, that ESS is not tiny, and that no covariate
mechanically encodes treatment/eligibility (look at the top feature importances in
{cmd:diagnostics.json}).{p_end}
{p 8 10}The summary is the starting point; the full design-stage protocol is in
{help garf##diagnostics:Diagnostics}.{p_end}


{marker fmq}{...}
{title:Forest Match Quality (FMQ)}

{p 4 8}FMQ summarises, on a 0-1 scale, how close each treated unit's donors are {it:in the forest's
own metric}, relative to what the data make possible. For treated unit {it:i} and control {it:j}, let{p_end}
{p 12 12}{cmd:c(i,j)} = share of trees in which {it:i} and {it:j} land in the same leaf (0-1; column {cmd:cooccurrence}){p_end}
{p 4 8}Then, for each treated unit:{p_end}
{p 12 12}{cmd:S(i)}    = sum over i's donors of {cmd:weight(i,j) * c(i,j)}  (how often, on average, i shares a leaf with its donors){p_end}
{p 12 12}{cmd:cbar(i)} = average of {cmd:c(i,j)} over {it:all} controls  (what random donors would achieve){p_end}
{p 12 12}{cmd:Smax(i)} = maximum of {cmd:c(i,j)} over {it:all} controls  (the closest untreated unit in the data){p_end}
{p 12 12}{cmd:FMQ(i)}  = (S(i) - cbar(i)) / (Smax(i) - cbar(i)){p_end}
{p 4 8}FMQ is the share of the {it:attainable} improvement, from "a random control" to "the closest untreated
unit the data contain", that the matching achieved: {bf:0} = donors no closer than random controls; {bf:1} =
every donor as close as the best one available; below 0 = worse than random. Normalising by the best
attainable donor keeps the scale comparable across sample sizes and forest settings: with
tens of thousands of units and small leaves, no pair shares a leaf in most trees, so the earlier
twin-scaled version, {cmd:(S(i) - cbar(i)) / (1 - cbar(i))}, stayed close to 0 even for good matches. It is still
returned as {cmd:r(fmq_twin_mean)} and per block as {cmd:fmq_twin}. If no control is closer to {it:i} than the
average one ({cmd:Smax(i) = cbar(i)}), FMQ is undefined for that block, excluded from the means and counted
in {cmd:r(fmq_undefined)}. FMQ uses only pre-treatment information (never the outcome effect), so it can be
used freely at the design stage.{p_end}

{p 4 8}{bf:Reported as}: {cmd:r(fmq_mean)} and {cmd:r(fmq_p10)} (mean and 10th percentile over {it:matched}
treated units) and {cmd:r(fmq_mean_all)} (mean over {it:all} treated units, with unmatched ones scored 0).
Per block: variable {cmd:fmq} in the matched panel and, in {cmd:blocks.csv}, {cmd:fmq}, {cmd:fmq_twin} and
its components {cmd:match_cooccurrence} (S), {cmd:random_cooccurrence} (cbar) and {cmd:best_cooccurrence}
(Smax); per pair: variable {cmd:cooccurrence}.{p_end}

{p 4 8}{bf:How to read it {hline 1} important.}{p_end}
{p 8 10}1. {bf:FMQ strongly rewards fewer, closer donors.} A smaller {opt k()}, a tighter {opt caliper()} or a
larger {opt wpower()} systematically raises it, and {cmd:k(1)} typically gives the highest value. It measures the {it:closeness} (bias) side of
the design only. {bf:Never choose a design by FMQ alone}: read it {it:relative to the number of donors} and
together with the effective sample size ({cmd:r(mean_ess)}), which measures the variance side. A good
design gives up little FMQ for a large gain in ESS. Comparing FMQ between designs with very different
{opt k()} or {opt wpower()} without looking at ESS is misleading.{p_end}
{p 8 10}2. {bf:With a caliper}, treated units without donors are left out of {cmd:r(fmq_mean)}, which can rise
simply because the hardest units were dropped. Compare caliper designs on {cmd:r(fmq_mean_all)} (unmatched
units count as 0), together with {cmd:r(n_treated_matched)} and ESS. Without a caliper the two means coincide.{p_end}
{p 8 10}3. {bf:Relative, not absolute.} There is no universal "good" threshold. Since 1 is now the best
the data allow, a low value means the design is far from the closest available donors, not that no close
donors exist; read the overlap diagnostics for the latter. Designs that differ only in
{opt k()}, {opt caliper()}, {opt wpower()} or replacement (same covariates, forest settings and seed) share the
same forest, so their FMQ is fully comparable. Across covariate sets or forest settings it is only roughly
comparable.{p_end}
{p 8 10}4. {bf:It sees only what the forest sees}: closeness in the dimensions that predict treatment. It does not
replace the balance table ({cmd:balance.csv}), period-by-period balance, or the pre-treatment outcome path.{p_end}
{p 8 10}5. {bf:Low-FMQ blocks} ({cmd:fmq} close to 0) flag treated units whose donors are no better than random:
list them with {cmd:random_cooccurrence}, which tells a {it:typical} unit (high {cmd:cbar}: any control is about as
close) from an {it:isolated} one (low {cmd:cbar} and low overlap mass), and say in the analysis whether they
are kept.{p_end}

{pstd}Example {hline 1} the closeness/size trade-off across the number of donors:{p_end}
{phang2}{cmd:. foreach k in 1 5 10 25 50 {c -(}}{p_end}
{phang2}{cmd:.     garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) k(`k') noload}{p_end}
{phang2}{cmd:.     di "k=`k'  FMQ=" %5.3f r(fmq_mean) "  ESS=" %5.1f r(mean_ess)}{p_end}
{phang2}{cmd:. {c )-}}{p_end}


{title:Studying heterogeneous effects}

{p 4 8}
Because GARF matches each treated unit to its {it:own} donors and keeps that link
({cmd:block_id} / {cmd:treated_id}), you can split the treated units into groups (by size,
sector, region, ...) and the right controls follow automatically: each treated unit's matched
donors travel with it into its subgroup. Tag every block by a {it:pre-treatment} characteristic
of its treated unit, then estimate the DiD within each group — no re-matching needed.{p_end}

{pstd}Conceptual example — effects by firm size:{p_end}
{phang2}{cmd:. gen s = employees if treated_grp==1 & year<2016}   // treated rows, pre-period{p_end}
{phang2}{cmd:. bysort block_id: egen size_pre = mean(s)}          // size of the whole block{p_end}
{phang2}{cmd:. gen byte big = size_pre > 50}{p_end}
{phang2}{cmd:. reghdfe Y W [pw=weight] if big==1, absorb(garf_unit year) vce(cluster garf_block)}{p_end}
{phang2}{cmd:. reghdfe Y W [pw=weight] if big==0, absorb(garf_unit year) vce(cluster garf_block)}{p_end}

{p 4 8}
Each regression uses only the treated units in that group and their own matched controls. Watch
the group sizes ({cmd:tab big if treated_grp==1}) — small groups give noisy estimates.{p_end}


{title:What loading adds}

{p 4 8}After a default run the matched panel is in memory. Two kinds of variables come along:
those carried straight through from the engine's {cmd:matched_panel.csv} ({cmd:proximity},
{cmd:proximity_minmax}, {cmd:weight}, {cmd:rank}, {cmd:leaf_distance} — see
{help garf##proximity:How proximity is calculated} above; on the treated unit's own rows
{cmd:proximity} and {cmd:proximity_minmax} are missing, while {cmd:rank} = 0, {cmd:weight} = 1 and
{cmd:leaf_distance} = 0), and a few Stata-side conveniences the command adds on top:{p_end}
{p 8 10}{cmd:weight}       design weight, ~ proximity^{opt wpower()}, summing to 1 per block (use as {cmd:[pw=weight]}){p_end}
{p 8 10}{cmd:proximity_minmax}  match quality on a 0-1 scale (1 = the closest pair in this run){p_end}
{p 8 10}{cmd:cooccurrence}  share of trees in which the donor and its treated unit share a leaf (missing on treated rows){p_end}
{p 8 10}{cmd:fmq}          Forest Match Quality of the block (same value on all rows of a block){p_end}
{p 8 10}{cmd:treated_id}   the treated unit each control is matched to (= {cmd:block_id}){p_end}
{p 8 10}{cmd:treated_grp}  1 = treated unit, 0 = selected control (donor){p_end}
{p 8 10}{cmd:control_sel}  1 = this row is a selected control (donor){p_end}
{p 8 10}{cmd:garf_unit}    numeric id of {cmd:unique_id} — use in {cmd:absorb()}{p_end}
{p 8 10}{cmd:garf_block}   numeric id of {cmd:block_id} — use in {cmd:vce(cluster )}{p_end}
{p 8 10}{cmd:post}         1 if {it:time} >= T0{p_end}


{title:Setup (once)}

{phang2}{cmd:. global GARF_HOME "/path/to/GARF-binary"}{p_end}
{p 4 8}Put {cmd:garf.ado} and {cmd:garf.sthlp} on the ado-path — note the {bf:/stata} folder:{p_end}
{phang2}{cmd:. adopath + "$GARF_HOME/integrations/stata"}{p_end}
{p 4 8}To make it permanent, put the {cmd:global} line in your {cmd:profile.do} and either the
{cmd:adopath} line too, or copy {cmd:garf.ado} and {cmd:garf.sthlp} into your PERSONAL folder
(see {cmd:sysdir}). Keep the package folder itself in a fixed place: {cmd:garf} finds the engine
through {cmd:$GARF_HOME}.{p_end}
{p 4 8}{bf:Check the version.} Stata looks for {cmd:garf.ado} in the {it:current} folder before the
ado-path, so an old copy left in your working folder silently takes precedence. {cmd:which garf}
must report {bf:version 1.1}; after replacing the files, type {cmd:discard} (or restart Stata) so that
a version already loaded in memory is dropped.{p_end}


{marker diagnostics}{...}
{title:Diagnostics}

{p 4 8}The package folder {cmd:DIAGNOSTICS} contains commented do-files, one per diagnostic, for the
full design-stage protocol of the GARF working paper: feature-level balance by variable and
component, variance ratios, balance over time, the pre-treatment outcome path, an out-of-window
placebo, effective sample size and donor reuse, overlap, Forest Match Quality by block, forest
fit and feature importance, a sensitivity grid over {opt k()}, {opt caliper()}, {opt wpower()},
{opt noreplacement} and {opt seed()}, and the FMQ-ESS frontier. All of them use pre-treatment
information only. On the demonstration panel, from a project folder:{p_end}
{phang2}{cmd:. do "$GARF_HOME/DIAGNOSTICS/00_master.do"}{p_end}
{p 4 8}Outputs (tables in the log, five figures as PDF) go to {cmd:./garf_diagnostics}. For your own
data, copy the folder, edit only {cmd:00_settings.do}, set {cmd:global GARF_DIAG} to the copy and run
its {cmd:00_master.do}. The widest design of the sensitivity grid is set by {cmd:KWIDE} in
{cmd:00_settings.do}, default {cmd:k(1000)}: on large donor pools the uncapped forest kernel ({opt nokcap}) can
give thousands of donors per treated unit and a matched panel of tens of millions of rows; set
{cmd:global KWIDE "nokcap"} for small pools. See {cmd:DIAGNOSTICS/README.txt}.{p_end}


{title:Stored results}

{p 4 8}Scalars in {cmd:r()}: {cmd:r(n_treated)}, {cmd:r(n_controls)}, {cmd:r(n_treated_matched)},
{cmd:r(n_donors)}, {cmd:r(mean_smd_raw)}, {cmd:r(mean_smd_matched)}, {cmd:r(max_smd_raw)},
{cmd:r(max_smd_matched)}, {cmd:r(mean_ess)}, {cmd:r(fmq_mean)}, {cmd:r(fmq_p10)}, {cmd:r(fmq_mean_all)},
{cmd:r(oob_bacc)}, {cmd:r(oob_auc)}, {cmd:r(oob_score)}, {cmd:r(wpower)},
{cmd:r(cv_score)} (if tuned). The matched |SMD| compares the matched treated units (those with at
least one donor) with their weighted donors; the raw |SMD| uses all treated and all controls; both
use the same (raw, pooled) standard deviation. When files are saved: {cmd:r(outdir)}, {cmd:r(matched_panel)}.{p_end}

{title:Hyperparameter Tuning}

{p 4 8}The random forest underneath GARF can be tuned to specific use cases, to
identify the best fitting model. Parameters can be set by hand (section Forest above). Hyperparameter tuning, via
the {cmd:tune} option, is a data-driven strategy to identify the best fitting model across a wide range of
parameters: a randomized search ({opt niter()} candidates) with stratified {opt cvfolds()}-fold cross-validation
over units. It is computationally intensive. Candidates are ranked by {cmd:tuningmetric()}, always computed on the
forest's prediction of {it:treatment} (never on the outcome or the effect): {cmd:accuracy} is the share of correct
guesses, which can be misleading when treated units are few (a model that always predicts "control" scores high);
{cmd:f1_treated} is the F1 score of the treated class, which targets the rare class but can favour over-eager
models; {cmd:combined} (default) is 0.5 x F1 of the treated class + 0.5 x balanced accuracy (the average of the
true-positive and true-negative rates). The weights are set in {cmd:combined_metric()} near the top of
{cmd:garf_binary.py}. Because the forest that best predicts treatment is not necessarily the one that gives the best
matched sample, re-check balance, ESS and overlap after tuning.
Tuning the hyperparameters of a Random Forest should not be done lightly, extensive documentation on the specific Random Forest
algorithm used by GARF and all possible parameters are available at: 
https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html{p_end}

{title:Examples}

{pstd}Default (k-NN) — run, load the matched panel into memory, save nothing:{p_end}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) tinvariant(region)}{p_end}
{phang2}{cmd:. reghdfe Y W [pw=weight], absorb(garf_unit year) vce(cluster garf_block)}{p_end}

{pstd}Caliper selection (all controls within a proximity of each treated; caliper is on the
0-1 normalised scale, {opt caliper(0.3)} keeps a fairly loose middle band — tighten toward 1 for
stricter matches):{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) caliper(0.3) nokcap}{p_end}

{pstd}Both together (the k nearest among those within the caliper):{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) caliper(0.3) k(50)}{p_end}

{pstd}Also keep the CSV files:{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) store}{p_end}

{pstd}More weight on the closest donors (check ESS and balance), or uniform k:1 weights:{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) wpower(2)}{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) k(10) wpower(0)}{p_end}

{pstd}With hyperparameter tuning (slow):{p_end}
{phang2}{cmd:. garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) tvarying(x1 x2) tinvariant(region) tune}{p_end}


{title:Notes}

{p 4 8}A unit is treated if its treatment column is ever 1. With the standard coding (1 for a treated
unit in its treated years, 0 otherwise) {cmd:W} is already the treated-times-post indicator, so a
regression of the outcome on {cmd:W} with unit and time fixed effects returns the ATT directly. An
ever-treated indicator (1 in all years of treated units) gives the same matching; then build the
interaction yourself ({cmd:c.treated_grp#c.post}). The design assumes one onset T0: with staggered
adoption, run {cmd:garf} cohort by cohort on the cohort's treated units plus the never-treated.
Units without a complete pre-window are dropped by the engine (with a note).{p_end}


{title:Author and feedback}

{p 4 8}
This is {bf:version 1.1}, the first public release.{p_end}

{p 4 8}Special thanks for sparking the need for this package, and for all the time spent trying to
understand and debug what I was building, go to my friend and colleague{p_end}
{p 4 8}dott. Giulio Breglia, from GSSI - Gran Sasso Science Institute.{p_end}

{p 4 8}
Comments, bugs and suggestions are very welcome:{p_end}
{p 4 8}
Federico Fantechi, PhD{p_end}
{p 4 8}
Department of Economics, Business and Statistics — University of Palermo{p_end}
{p 4 8}
{cmd:federico.fantechi@unipa.it}{p_end}
