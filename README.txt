=======================================================================
GARF (binary) — forest-proximity matching for panel counterfactuals
version 1.1
Federico Fantechi, PhD — University of Palermo — federico.fantechi@unipa.it
=======================================================================

Version history: see CHANGELOG.md.

A simplified, binary-only build of GARF, packaged so a colleague can test it
in a few minutes. You give it a panel; it fits a random forest on the binary
treatment, turns the forest into a proximity kernel, and returns a MATCHED
PANEL (each treated unit paired with its k nearest control "donors"), plus
balance diagnostics. It does NOT estimate the treatment effect — you run your
own DiD / event study on the matched panel it produces. Keeping the outcome out
of the tool is deliberate: nothing here can be tuned on the effect you are
trying to measure.


-----------------------------------------------------------------------
1. FOLDER LAYOUT
-----------------------------------------------------------------------

  README.txt          this file
  INSTALL.txt         how to install Python (the only prerequisite)
  requirements.txt    the Python packages (installed automatically)
  setup.sh / setup.bat   optional one-time environment setup

  engine/             THE TOOL — run it on your own data
    garf_binary.py      the engine (one file)
    config.yaml         the only thing you edit
    run.sh / run.bat    run on your data
    input/              <-- put your panel CSV here (default: input/panel.csv)
    output/             <-- results are written here

  DIAGNOSTICS/        Stata do-files for the full diagnostic protocol
    00_master.do        runs every diagnostic in order (see its README.txt)

  demo/               TESTING — synthetic data, to check the install works
    make_example.py     generates a synthetic panel
    config.yaml         demo config
    run_demo.sh / run_demo.bat
    input/  output/     demo data / demo results

The engine and the demo share ONE virtual environment (.venv), created in this
main folder on first run. Nothing is installed system-wide.


-----------------------------------------------------------------------
2. QUICK START
-----------------------------------------------------------------------

Prerequisite: Python 3.10+ (see INSTALL.txt if you don't have it).

STEP 1 — test the install (synthetic data, ~1 minute):
  Windows:        open the  demo  folder, double-click  run_demo.bat
  macOS / Linux:  in a terminal:   cd demo   then   ./run_demo.sh
                  (first time only:  chmod +x ../*.sh ../**/*.sh )

STEP 2 — run on your data:
  a) put your panel CSV in   engine/input/   (simplest: name it  panel.csv )
  b) open  engine/config.yaml  and set the column names / T0 / k
  c) Windows:       open  engine, double-click  run.bat
     macOS / Linux: cd engine   then   ./run.sh
  Results appear in  engine/output/ .

The first run builds the shared .venv and installs the packages; later runs
start immediately. Prefer a separate setup step? Run setup.bat / ./setup.sh once.


-----------------------------------------------------------------------
3. THE CONFIG FILE  (engine/config.yaml)
-----------------------------------------------------------------------

Your data must be a LONG-format panel CSV: one row per unit-year, with a binary
treatment column. Key settings:

  data_path            path to your CSV (default: input/panel.csv)
  id_col, year_col, treat_col, outcome_col   your column names
  x_time_varying       covariates that change over time
                       -> get lags AND trajectory features
  x_time_invariant     covariates that don't change (strings are fine)
  treatment_year       T0, the first treated period
  n_pre_years          length of the pre-window used for matching: [T0-n, T0-1]
  k_matches            keep at most this many nearest donors (k-NN); null or 0 = no cap
  caliper              keep only donors with MIN-MAX NORMALISED proximity >= caliper
                       (range [0,1]; 0 = off). Normalised proximity rescales the raw
                       proximity matrix once, globally, over this run's treated x
                       control pairs, so a caliper value means the same thing across
                       specifications (see section 5). A donor must also have actually
                       co-occurred with the treated unit at least once (raw proximity
                       > 0), regardless of the caliper.
                       k_matches and caliper COMBINE: e.g. caliper>0 with k_matches:null
                       keeps every control within the caliper; with a finite k it keeps
                       the k closest among those within the caliper.
  with_replacement     true = donors reusable across treated units (recommended);
                       false = each control used at most once, treated units matched
                       hardest-to-match first
  weight_power         donor weights ~ proximity^weight_power (default 1; 0 = uniform;
                       >1 = more weight on the closest donors)
  out_dir              where results go (default: output)
  progress             true = live progress + timing in the terminal

Treatment column: 1 for a treated unit in its treated years, 0 otherwise; a unit
is counted as "treated" if it is ever 1. Units without a complete pre-window are
dropped (with a warning saying how many).


-----------------------------------------------------------------------
4. HYPERPARAMETER TUNING  (fast vs slow)
-----------------------------------------------------------------------

The forest can either use fixed defaults (fast) or be tuned by cross-validated
model selection (slower, sometimes sharper matches). Controlled in config.yaml:

  tune: false     FAST. Uses the fixed forest below (n_estimators, min_samples_leaf,
                  max_features, max_depth). This is the default.

  tune: true      SLOW. Randomized cross-validated search over forest
                  hyperparameters, selected on the forest's OWN prediction target
                  (never on the treatment effect). Also set:
                    tuning_metric: accuracy | f1_treated | combined
                    n_iter_search: how many candidates to try (e.g. 40)
                    cv_folds: cross-validation folds (e.g. 5)

  tuning_metric = "combined" (default) blends F1 of the treated class with
  balanced accuracy — a better model-selection criterion than plain accuracy when
  treatment is rare. The exact formula is the small combined_metric() function
  near the top of garf_binary.py; edit it there if your pipeline uses a different
  one.

Rule of thumb: leave tuning OFF for a first look and for large panels; turn it ON
when you want the forest itself chosen by cross-validation. Tuning never sees the
ATT, by design.


-----------------------------------------------------------------------
5. WHAT IT BUILDS  (features)
-----------------------------------------------------------------------

For every unit, from the pre-window only:
  - raw lags of Y and each time-varying X (Y__L5 ... Y__L1);
  - trajectory features per variable: mean, delta (last minus first), growth
    (symmetric: delta over the mean of |first| and |last|), slope (OLS), vol
    (sd of first differences), curv (2nd-half slope minus 1st-half slope; with
    an odd window the 2nd half is one period longer), dip (last value minus
    mean of the earlier pre-period). The first and last levels are not repeated:
    they are already the lags Y__L5 and Y__L1 -- so L+7 features per variable;
  - time-invariant X as-is (categoricals -> integer codes; forests split on
    order statistics, so no scaling / one-hot decision is needed).

The (raw) proximity is
  proximity(i,j) = (1/B) * sum_b  1{leaf_b(i) = leaf_b(j)} / |leaf_b(i) among donors|.
The leaf-size normalisation is not cosmetic: it is what turns co-occurrence into
the adaptive-nearest-neighbour kernel of Lin & Jeon (2006) / Athey, Tibshirani &
Wager (2019). Higher = closer. Its absolute size is NOT comparable across runs
(it depends on min_samples_leaf, n_estimators, and the number of covariates).

The min-max NORMALISED proximity (new in beta 1.1; what the caliper filters on) is
  proximity_minmax(i,j) = (proximity(i,j) - min) / (max - min)
where min/max are the smallest/largest raw proximity anywhere in the WHOLE
treated x control matrix for this run (reported in diagnostics.json). This is a
GLOBAL rescaling, not a per-treated-unit one, on purpose: normalising each
treated unit by its own best match would send every unit's top donor to 1.0
regardless of how poor it actually is, hiding exactly the units with no good
match. With one global (min, max), proximity_minmax = 1 is reserved for the
single closest pair in the whole run, so a caliper value means the same thing
across specifications.

weight(i,j), used for estimation, renormalises the RAW proximities, raised to
weight_power (default 1), to sum to 1 within each treated unit's own selected
donor set:
  weight(i,j) = proximity(i,j)^g / sum_l proximity(i,l)^g,   g = weight_power.
It is a within-block convex combination, not a 0-1 match-quality score (that is
proximity_minmax).


-----------------------------------------------------------------------
6. OUTPUTS  (engine/output/)
-----------------------------------------------------------------------

  matched_panel.csv   long panel, treated + donor rows, with unique_id, block_id
                      (= treated unit), role, year, W, Y, X's, weight,
                      proximity, proximity_minmax, leaf_distance,
                      cooccurrence (donor rows), fmq (per block).
                      Ready for a weighted DiD, e.g. in R:
                      feols(Y ~ treat*post | unique_id + year, weights = weight)
  matches.csv         treated_id, donor_id, unique_id, rank, proximity,
                      proximity_minmax, weight, leaf_distance, cooccurrence.  weight =
                      proximity normalised within the treated block (sums to
                      1) -> use it in the weighted DiD.  proximity_minmax is
                      the 0-1 normalised match quality (see section 5).
                      cooccurrence = share of trees in which the pair shares a leaf.
  blocks.csv          one row per treated unit: matched (0/1), n_donors, ess,
                      match_cooccurrence S(i), random_cooccurrence cbar(i),
                      best_cooccurrence Smax(i), fmq, fmq_twin, best_proximity_minmax, overlap_mass (share of trees in which
                      the unit's leaf contains at least one control).
  balance.csv         standardized mean difference per feature, raw (all treated vs
                      all controls) vs matched (matched treated vs weighted donors),
                      both on the same raw pooled SD.
  diagnostics.json    config used + n treated/donors, mean & max |SMD| before/after,
                      mean ESS per treated, weight_power, FMQ summaries (fmq_mean,
                      fmq_p10, fmq_median, fmq_mean_all, fmq_twin_mean,
                      fmq_undefined), forest out-of-bag
                      accuracy, balanced accuracy and AUC (and CV score if tuned),
                      top-15 feature importances.
  features.csv        the unit-level design matrix, for inspection.

FOREST MATCH QUALITY (FMQ; normalisation changed in beta 1.5)
For treated unit i and control j, c(i,j) = share of trees in which they land in
the same leaf (0-1). For each treated unit:
    S(i)    = sum over i's donors of weight(i,j) * c(i,j)
    cbar(i) = mean of c(i,j) over ALL controls      (= random donors)
    Smax(i) = max  of c(i,j) over ALL controls      (= best attainable donor)
    FMQ(i)  = (S(i) - cbar(i)) / (Smax(i) - cbar(i))
0 = donors no closer than random controls; 1 = every donor as close as the
closest untreated unit in the data. Undefined (missing, counted in
fmq_undefined) if Smax(i) = cbar(i). The beta 1.2-1.4 version,
(S(i) - cbar(i)) / (1 - cbar(i)), is still reported as fmq_twin: its scale
depends on N and the forest settings. Reported as fmq_mean / fmq_p10 (matched treated units) and
fmq_mean_all (all treated, unmatched = 0).
  * FMQ STRONGLY REWARDS FEWER, CLOSER DONORS: a smaller k, a tighter caliper or
    a larger weight_power systematically raises it, and k = 1 typically
    gives the highest value. It measures
    the closeness (bias) side only. Never pick a design on FMQ alone: read it
    relative to the number of donors and together with the effective sample
    size (ESS, the variance side). Good designs give up little FMQ for a large
    gain in ESS.
  * With a caliper, compare designs on fmq_mean_all (unmatched treated units
    count as 0), not fmq_mean, which rises when hard units are dropped.
  * Relative, not absolute: no universal threshold. Designs that differ only in
    k / caliper / weight_power / replacement share the same forest and are fully
    comparable; across covariate sets or forest settings only roughly.
  * It sees only what the forest sees (dimensions that predict treatment): it
    does not replace the balance table or the pre-treatment outcome path.

While running, the engine prints each stage and how long it took, with a progress
bar over the proximity computation (the heaviest step).


-----------------------------------------------------------------------
7. TROUBLESHOOTING
-----------------------------------------------------------------------

  "Python was not found"        install it (INSTALL.txt), reopen the terminal.
  macOS: .sh opens a text editor  run it from Terminal (./run_demo.sh) or rename
                                  to .command to double-click.
  macOS: "permission denied"     chmod +x ../*.sh ../**/*.sh  once, then retry.
  "Columns not found ..."        the names in config.yaml don't match your CSV
                                 headers; the error lists what it found — copy exactly.
  everything dropped             check treatment_year, n_pre_years, and that units
                                 have complete pre-window data.
  slow                           lower n_estimators (e.g. 500) and keep tune: false.

-----------------------------------------------------------------------
8. USING IT FROM STATA / R
-----------------------------------------------------------------------
The engine can be driven from Stata or R without losing any features — they act
as front-ends while Python stays the compute engine.

  integrations/stata/   garf.ado + help + example.do + README_stata.txt
                        (a native 'garf' command; run setup once, then call it)
  DIAGNOSTICS/          the design-stage diagnostic protocol, as do-files that
                        call garf (see DIAGNOSTICS/README.txt)

An R wrapper (reticulate + shell fallback) can be added the same way.


Simplified binary build. For the full version (continuous-dose fitting, honest
splitting, and the SC/DiD fitting-target theory), use the complete GARF package.
