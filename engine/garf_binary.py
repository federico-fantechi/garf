"""
GARF (binary) — forest-proximity matching for panel counterfactuals.
====================================================================
Version: 1.1 (first public release; identical to development build beta 1.5).

beta 1.5 changes:
  * CHANGED: Forest Match Quality is now normalised by the best match ATTAINABLE
    in the data, not by forest twins. With
        c(i,j)  = share of trees in which treated i and control j share a leaf,
        S_i     = sum_j w_ij c(i,j)        (weighted co-occurrence with i's donors)
        cbar_i  = mean_j c(i,j) over ALL controls      (random-control benchmark)
        Smax_i  = max_j c(i,j) over ALL controls       (best attainable donor)
    the new
        FMQ_i   = (S_i - cbar_i) / (Smax_i - cbar_i)   in [0, 1]
    is the share of the attainable gain over random controls that the design
    captures: 0 = donors no closer than random controls, 1 = every donor as close
    as the closest untreated unit in the data. The old twin-scaled version,
        FMQ_twin_i = (S_i - cbar_i) / (1 - cbar_i),
    could not approach 1 in large samples (with tens of thousands of units and
    small leaves no pair shares a leaf in most trees), so its level depended on
    N and min_samples_leaf. It is still reported, as fmq_twin.
    New outputs: blocks.csv gains best_cooccurrence (Smax_i) and fmq_twin;
    diagnostics gain fmq_twin_mean and best_cooccurrence_mean. FMQ still rewards
    concentration (fewer / closer donors): read it together with ESS.
    If Smax_i = cbar_i (no control closer than the average one), FMQ_i is
    undefined (missing) and excluded from the summaries.
  * DIAGNOSTICS: the widest design of the sensitivity grid is now a setting,
    KWIDE, with default k(1000) (was the uncapped forest kernel, which on large
    donor pools can give thousands of donors per treated unit); the
    balance-over-time figure wraps at four panels per row.

beta 1.4 changes:
  * No change to the computations: results are identical to beta 1.3.
  * The package's testing do-files ("MY TESTING DOs") are replaced by the
    DIAGNOSTICS folder: commented Stata do-files for the full design-stage
    diagnostic protocol (balance, variance ratios, balance over time, outcome
    path, out-of-window placebo, ESS and donor reuse, overlap, FMQ, forest
    importance, sensitivity grid and FMQ-ESS frontier).

beta 1.3 changes:
  * CHANGED: the trajectory descriptor no longer contains `__first` and `__last`.
    They were exact copies of the window's end-point lags (`__L<n>` and `__L1`),
    added no information, made the end-points twice as likely to be offered at a
    split, and were counted twice in the balance summary. Each time-varying
    variable (outcome included) now yields L lags + 7 descriptors (L+7 features;
    L+9 before). Results are therefore not numerically identical to beta 1.2.
  * CHANGED: `__growth` is now the SYMMETRIC (mid-point) growth
        (last - first) / ((|last| + |first|) / 2),
    bounded in [-2, 2] and defined for zero and negative values (it is 0 when
    both end-points are 0). Beta 1.2 used last/first - 1, which was undefined
    when the first value was 0, meaningless with negative values, and then
    silently median-imputed. When both end-points are positive the new measure
    is a strictly increasing transformation of the old one, so it ranks units
    identically and a forest splits on it in the same way.
  * NEW: honest fit statistics of the forest, computed from its OUT-OF-BAG
    predictions: `oob_balanced_accuracy` (mean of the true-positive and
    true-negative rates; 0.5 = no better than chance, whatever the treated
    share) and `oob_auc`. `oob_score` (plain OOB accuracy) is still reported.
  * REMOVED: `auc_full_sample` (in-sample AUC; ~1 on almost any data, hence
    uninformative). Use `oob_auc`.
  * REMOVED: `min_proximity` (legacy synonym of `caliper`). Old configs that
    still set it now stop with a message asking to use `caliper` instead.

beta 1.2 changes:
  * NEW: weight_power (gamma >= 0, default 1). Donor weights are
        w_ij = proximity(i,j)**gamma / sum_l proximity(i,l)**gamma
    within each treated unit's selected donors. gamma = 1 reproduces beta 1.1
    exactly (proportional to proximity); gamma = 0 gives uniform k:1 weights;
    gamma > 1 concentrates weight on the closest donors (gamma -> infinity
    tends to nearest-neighbour matching).
  * FIX: matching WITHOUT replacement now visits treated units hardest-to-match
    first (increasing best proximity), as the code comment always stated;
    beta 1.1 visited them easiest-first.
  * FIX: the matched SMD now compares the matched treated units with their
    donors (beta 1.1 averaged over ALL treated units, including any left
    without donors by a caliper). Unchanged when every treated unit is matched.
  * NEW: Forest Match Quality (FMQ). For each treated unit i, with
        c(i,j) = share of trees in which i and control j share a leaf,
        S_i    = sum_j w_ij c(i,j)          (weighted co-occurrence with its donors)
        cbar_i = mean_j c(i,j) over ALL controls (random-control benchmark)
        FMQ_i  = (S_i - cbar_i) / (1 - cbar_i)
    0 = donors no closer than random controls, 1 = donors indistinguishable
    from i for the forest. FMQ mechanically rewards concentration (smaller k,
    tighter caliper, larger weight_power): always read it together with the
    effective sample size (ESS), and compare designs in relative terms.
    Reported as fmq_mean / fmq_p10 (matched treated) and fmq_mean_all (all
    treated, unmatched = 0); per-pair c(i,j) in matches.csv ("cooccurrence"),
    per-block values in the new blocks.csv and as "fmq" in matched_panel.csv.

beta 1.1 change (read this if you used beta 1.0): the caliper now filters on
a MIN-MAX NORMALISED proximity (0 = the least similar treated-control pair
observed in this run, 1 = the most similar), not on the raw leaf-proximity
kernel. Raw proximity magnitudes depend on min_samples_leaf / n_estimators /
the number of covariates and are not comparable across specifications, which
made a raw-scale caliper hard to pick sensibly. A caliper value calibrated
against beta 1.0 does NOT carry over -- re-pick it on this run's output (see
`proximity_minmax` in matches.csv / diagnostics.json). The raw `proximity`
column is unchanged and still reported alongside it.

Simplified, binary-only build of GARF, meant to be handed to colleagues to
test with as little friction as possible.

Pipeline:

    panel data  ->  unit-level pre-treatment features (levels + trajectory shape)
                ->  random forest CLASSIFIER on the binary treatment W
                    (optionally hyperparameter-tuned by cross-validation)
                ->  leaf-proximity kernel (raw), then min-max normalised to [0, 1]
                ->  donor selection: k-nearest and/or a caliper on the
                    NORMALISED proximity (with / without replacement)
                ->  matched panel + balance/ESS diagnostics, ready for DiD

It does NOT estimate the treatment effect. That stays in your hands, on purpose:
the tool never sees the outcome effect, so a hyperparameter cannot be tuned on it.

Usage
-----
    python garf_binary.py --init                 # write a template config.yaml
    python garf_binary.py --config config.yaml   # run

or from Python:

    from garf_binary import GarfConfig, run
    res = run(GarfConfig.from_yaml("config.yaml"))

Requires: pandas, numpy, scipy, scikit-learn, pyyaml (tqdm optional, for a nicer
progress bar).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    balanced_accuracy_score,
    f1_score,
    make_scorer,
    roc_auc_score,
)
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold

__version__ = "1.1"


# =============================================================================
# 0. PROGRESS / LOGGING HELPERS
# =============================================================================
# The engine prints what it is doing and how long each stage takes. The heavy
# proximity loop shows a progress bar: a real tqdm bar if tqdm is installed,
# otherwise a plain-text fallback so there is never a hard dependency.

def _log(msg: str) -> None:
    print(f"[GARF] {msg}", flush=True)


def _dt(t0: float) -> str:
    return f"{time.time() - t0:.1f}s"


try:
    from tqdm import tqdm as _tqdm

    def _progress(iterable, total=None, desc=""):
        return _tqdm(iterable, total=total, desc=f"[GARF]   {desc}",
                     ncols=72, leave=False)
except Exception:  # tqdm not installed -> lightweight fallback
    def _progress(iterable, total=None, desc=""):
        seq = iterable if total is not None else list(iterable)
        n = total if total is not None else len(seq)
        step = max(1, n // 40)

        def gen():
            for i, x in enumerate(seq, 1):
                if i % step == 0 or i == n:
                    pct = int(100 * i / n)
                    bar = "#" * (pct // 3)
                    sys.stdout.write(f"\r[GARF]   {desc} [{bar:<33}] {pct:3d}% ({i}/{n})")
                    sys.stdout.flush()
                yield x
            sys.stdout.write("\n")
            sys.stdout.flush()
        return gen()


# =============================================================================
# 1. CONFIG
# =============================================================================


@dataclass
class GarfConfig:
    # ---- data ---------------------------------------------------------------
    data_path: str = "input/panel.csv"
    id_col: str = "id"
    year_col: str = "year"
    treat_col: str = "W"                     # binary treatment (0/1)
    outcome_col: str = "Y"
    x_time_varying: Sequence[str] = field(default_factory=list)
    x_time_invariant: Sequence[str] = field(default_factory=list)

    # ---- design -------------------------------------------------------------
    treatment_year: int = 2000               # T0: first treated period
    n_pre_years: int = 5                     # pre-window length used for fitting
    include_treatment_year: bool = False     # if False, window is [T0-L, T0-1]

    # ---- forest defaults (used as-is when tuning is off) --------------------
    n_estimators: int = 1000
    min_samples_leaf: int = 5
    max_features: object = "sqrt"            # "sqrt" | "log2" | float in (0,1]
    max_depth: object = None                 # null for unlimited
    random_state: int = 42
    n_jobs: int = -1

    # ---- hyperparameter tuning (model selection) ----------------------------
    tune: bool = False                       # true -> randomized CV search (SLOW)
    tuning_metric: str = "combined"          # "accuracy" | "f1_treated" | "combined"
    n_iter_search: int = 40                  # candidates drawn in the search
    cv_folds: int = 5                        # cross-validation folds

    # ---- matching -----------------------------------------------------------
    # Two ways to pick each treated unit's controls; they COMPOSE:
    #   * k_matches : keep at most this many nearest donors (k-NN). null/0 = no cap.
    #   * caliper   : keep only donors with MIN-MAX NORMALISED proximity >= caliper
    #                 (0 = off). Normalised proximity is proximity rescaled ONCE,
    #                 GLOBALLY over every treated x control pair in this run, to
    #                 [0, 1]: 0 = least similar pair observed, 1 = most similar.
    #                 It lives in [0, 1] regardless of forest size/leaf size, which
    #                 the raw proximity kernel does not, so caliper values are
    #                 comparable across specifications. See `proximity_minmax` in
    #                 the outputs, and the guide's proximity section, for the exact
    #                 formula.
    # Defaults (k_matches=50, caliper=0) reproduce the original k-NN behaviour.
    k_matches: object = 50                    # int cap, or null / 0 = "no cap"
    caliper: float = 0.0                      # min NORMALISED proximity (0-1) to be an eligible donor
    with_replacement: bool = True
    # Donor weights within each treated unit: w_ij proportional to proximity**weight_power.
    #   1 = proportional to proximity (default; identical to beta 1.1)
    #   0 = uniform weights over the selected donors (classic k:1)
    #  >1 = more weight on the closest donors (-> nearest neighbour as it grows)
    weight_power: float = 1.0

    # ---- output / misc ------------------------------------------------------
    out_dir: str = "output"
    progress: bool = True                    # live progress + timing in terminal

    # ---- io -----------------------------------------------------------------
    @classmethod
    def from_yaml(cls, path: str) -> "GarfConfig":
        with open(path) as fh:
            raw = yaml.safe_load(fh) or {}
        known = {f for f in cls.__dataclass_fields__}
        if "min_proximity" in raw:
            raise ValueError(
                "min_proximity was removed in beta 1.3 (it was a legacy synonym of "
                "caliper). Delete it from the config and set `caliper` instead."
            )
        unknown = set(raw) - known
        if unknown:
            raise ValueError(
                f"Unknown config keys: {sorted(unknown)}.\n"
                f"Allowed keys: {sorted(known)}"
            )
        return cls(**raw)

    @property
    def forest_params(self) -> dict:
        return {
            "n_estimators": self.n_estimators,
            "min_samples_leaf": self.min_samples_leaf,
            "max_features": self.max_features,
            "max_depth": self.max_depth,
        }

    @property
    def proximity_threshold(self) -> float:
        """Effective caliper, on the [0, 1] MIN-MAX NORMALISED proximity scale.

        Donors need normalised proximity >= this (and, always, raw proximity > 0:
        a control that never once shared a leaf with the treated unit is never
        an eligible donor, caliper or not -- see select_donors()).
        """
        return float(self.caliper)

    @property
    def capped(self) -> bool:
        """True if the number of donors per treated is capped at k_matches."""
        return (self.k_matches is not None) and (int(self.k_matches) > 0)

    def validate(self) -> None:
        if self.k_matches is not None and int(self.k_matches) < 0:
            raise ValueError("k_matches must be >= 1, or null / 0 for no cap")
        if self.caliper < 0:
            raise ValueError("caliper must be >= 0")
        if self.caliper > 1:
            raise ValueError(
                "caliper must be in [0, 1]. As of beta 1.1 it "
                "filters on the MIN-MAX NORMALISED proximity (1 = the closest "
                "treated-control pair observed in this run), not the raw "
                "leaf-proximity kernel used in beta 1.0. A caliper value "
                "calibrated against the old raw scale will not transfer -- "
                "re-pick it by inspecting `proximity_minmax` in this run's "
                "matches.csv / diagnostics.json."
            )
        if float(self.weight_power) < 0:
            raise ValueError("weight_power must be >= 0 (1 = proportional to proximity, "
                             "0 = uniform weights, >1 = more weight on the closest donors)")
        if self.n_pre_years < 2:
            raise ValueError("n_pre_years must be >= 2 (trajectory features need a slope)")
        if self.tuning_metric not in {"accuracy", "f1_treated", "combined"}:
            raise ValueError(
                "tuning_metric must be 'accuracy', 'f1_treated', or 'combined'"
            )


# =============================================================================
# 2. TUNING METRICS
# =============================================================================
#
# `combined` is the plug-in slot. Replace the body of `combined_metric` with the
# exact formula from your pipeline; everything else is untouched.
#
COMBINED_F1_WEIGHT = 0.5      # weight on F1 of the treated (positive) class
COMBINED_BAL_WEIGHT = 0.5     # weight on balanced accuracy


def combined_metric(y_true, y_pred) -> float:
    """Weighted blend of F1 (positive class) and balanced accuracy.

    Rationale: treatment is usually rare, so plain accuracy is dominated by the
    control class and is a bad model-selection criterion for a *metric* that must
    discriminate treated-like from control-like units.
    """
    f1 = f1_score(y_true, y_pred, pos_label=1, zero_division=0)
    bal = balanced_accuracy_score(y_true, y_pred)
    return COMBINED_F1_WEIGHT * f1 + COMBINED_BAL_WEIGHT * bal


def _scorer(cfg: GarfConfig):
    return {
        "accuracy": "accuracy",
        "f1_treated": make_scorer(f1_score, pos_label=1, zero_division=0),
        "combined": make_scorer(combined_metric),
    }[cfg.tuning_metric]


def _param_grid():
    return {
        "n_estimators": [300, 500, 1000, 2000],
        "min_samples_leaf": [1, 2, 3, 5, 10, 20],
        "max_features": ["sqrt", "log2", 0.3, 0.5, 0.8],
        "max_depth": [None, 4, 6, 8, 12, 20],
    }


# =============================================================================
# 3. FEATURE CONSTRUCTION
# =============================================================================


def _traj_features(block: pd.DataFrame, var: str) -> dict:
    """Trajectory descriptors for one variable over the pre-window.

    These are what make the metric *panel-native*: the forest can split on the
    SHAPE of the pre-period, not just its level.
    """
    v = block[var].to_numpy(dtype=float)
    t = np.arange(len(v), dtype=float)
    out = {}
    # The window's end-points are already in the lags (__L<n> = first period,
    # __L1 = last period); beta 1.3 no longer duplicates them as __first/__last.
    out[f"{var}__mean"] = np.nanmean(v)
    out[f"{var}__delta"] = v[-1] - v[0]
    # Symmetric (mid-point) growth: bounded in [-2, 2], defined for zeros and
    # negative values; 0 when both end-points are 0.
    mid = (abs(v[-1]) + abs(v[0])) / 2.0
    out[f"{var}__growth"] = (v[-1] - v[0]) / mid if mid > 0 else 0.0
    if len(v) >= 2 and np.isfinite(v).all():
        out[f"{var}__slope"] = np.polyfit(t, v, 1)[0]
    else:
        out[f"{var}__slope"] = np.nan
    d = np.diff(v)
    out[f"{var}__vol"] = np.nanstd(d) if len(d) else np.nan
    # Curvature: slope over the second half of the window minus slope over the
    # first half. With an odd window the second half is one period longer
    # (e.g. L = 5: periods 1-2 vs 3-5).
    if len(v) >= 4 and np.isfinite(v).all():
        h = len(v) // 2
        s1 = np.polyfit(t[:h], v[:h], 1)[0]
        s2 = np.polyfit(t[h:], v[h:], 1)[0]
        out[f"{var}__curv"] = s2 - s1
    else:
        out[f"{var}__curv"] = np.nan
    out[f"{var}__dip"] = v[-1] - np.nanmean(v[:-1]) if len(v) >= 2 else np.nan
    return out


def build_features(df: pd.DataFrame, cfg: GarfConfig):
    """Panel -> one row per unit of pre-treatment features. Returns (X, meta)."""
    t0 = cfg.treatment_year
    last = t0 if cfg.include_treatment_year else t0 - 1
    first = last - cfg.n_pre_years + 1
    pre_years = np.arange(first, last + 1)

    pre = df[df[cfg.year_col].isin(pre_years)].copy()
    pre = pre.sort_values([cfg.id_col, cfg.year_col])

    tv_vars = [cfg.outcome_col] + list(cfg.x_time_varying)
    rows, meta_rows, dropped = [], [], []

    unit_W = df.groupby(cfg.id_col)[cfg.treat_col].max()

    for uid, block in pre.groupby(cfg.id_col, sort=False):
        block = block.set_index(cfg.year_col).reindex(pre_years)
        if block[tv_vars].isna().any().any():
            dropped.append(uid)
            continue
        feat = {cfg.id_col: uid}
        for var in tv_vars:
            for j, yr in enumerate(pre_years):          # raw lags
                feat[f"{var}__L{cfg.n_pre_years - j}"] = block.loc[yr, var]
            feat.update(_traj_features(block, var))
        for var in cfg.x_time_invariant:               # last observed pre value
            feat[var] = block[var].dropna().iloc[-1] if block[var].notna().any() else np.nan
        rows.append(feat)
        meta_rows.append({cfg.id_col: uid, "W": int(unit_W.loc[uid])})

    if dropped:
        warnings.warn(
            f"{len(dropped)} units dropped: incomplete pre-window "
            f"[{first}-{last}] on outcome/time-varying X."
        )

    X = pd.DataFrame(rows).set_index(cfg.id_col)
    meta = pd.DataFrame(meta_rows).set_index(cfg.id_col)

    for c in X.columns:
        if not pd.api.types.is_numeric_dtype(X[c]):
            X[c] = X[c].astype("category").cat.codes
    X = X.astype(float)
    if X.isna().any().any():
        n_bad = int(X.isna().any(axis=1).sum())
        warnings.warn(f"{n_bad} units have NaN features; median-imputing.")
        X = X.fillna(X.median())

    return X, meta


# =============================================================================
# 4. FOREST (binary; optional CV tuning)
# =============================================================================


def fit_forest(X: pd.DataFrame, meta: pd.DataFrame, cfg: GarfConfig):
    """Fit a random-forest classifier on the binary treatment. Returns (model, info).

    If cfg.tune is True, hyperparameters are selected by randomized cross-validated
    search on the forest's OWN target (never on the ATT). This is slower.
    """
    rng = cfg.random_state
    y = meta["W"].to_numpy()
    base = RandomForestClassifier(
        random_state=rng, n_jobs=cfg.n_jobs, oob_score=True,
        class_weight="balanced_subsample", **cfg.forest_params,
    )
    info = {"target": "binary", "n_fit": int(len(X)), "tuned": bool(cfg.tune)}

    if cfg.tune:
        cv = StratifiedKFold(n_splits=cfg.cv_folds, shuffle=True, random_state=rng)
        search = RandomizedSearchCV(
            base, _param_grid(), n_iter=cfg.n_iter_search, scoring=_scorer(cfg),
            cv=cv, random_state=rng, n_jobs=cfg.n_jobs, refit=True,
            verbose=1 if cfg.progress else 0,
        )
        search.fit(X, y)
        model = search.best_estimator_
        info["best_params"] = search.best_params_
        info["cv_score"] = float(search.best_score_)
        info["tuning_metric"] = cfg.tuning_metric
    else:
        model = base.fit(X, y)
        info["params"] = cfg.forest_params

    # Honest fit statistics from OUT-OF-BAG predictions. Plain OOB accuracy is
    # dominated by the majority class when treatment is rare, so balanced
    # accuracy and AUC are reported alongside it. (The in-sample AUC of
    # beta <= 1.2 was ~1 on almost any data and has been removed.)
    try:
        info["oob_score"] = float(model.oob_score_)
        oob = np.asarray(model.oob_decision_function_)[:, 1]
        ok = np.isfinite(oob)
        info["oob_n_scored"] = int(ok.sum())
        if ok.sum() and len(np.unique(y[ok])) == 2:
            info["oob_balanced_accuracy"] = float(
                balanced_accuracy_score(y[ok], (oob[ok] >= 0.5).astype(int)))
            info["oob_auc"] = float(roc_auc_score(y[ok], oob[ok]))
        else:
            info["oob_balanced_accuracy"] = float("nan")
            info["oob_auc"] = float("nan")
    except AttributeError:
        info["oob_balanced_accuracy"] = float("nan")
        info["oob_auc"] = float("nan")
    info["top_features"] = (
        pd.Series(model.feature_importances_, index=X.columns)
        .sort_values(ascending=False).head(15).round(4).to_dict()
    )
    return model, info


# =============================================================================
# 5. PROXIMITY
# =============================================================================


def leaf_proximity(model, X: pd.DataFrame, treated_mask: np.ndarray,
                   show: bool = True, return_cooccurrence: bool = False):
    """pi(i,j) = (1/B) * sum_b 1{leaf_b(i)=leaf_b(j)} / |leaf_b(i) among donors|.

    Leaf-size normalisation is NOT optional: it is what turns co-occurrence into
    the adaptive-nearest-neighbour KERNEL of Lin & Jeon (2006) / Athey et al.
    (2019). Rows sum to ~1 over the donor pool.
    Returns an (n_treated x n_control) dense array.

    If return_cooccurrence is True, also returns the RAW co-occurrence matrix
        c(i,j) = (1/B) * sum_b 1{leaf_b(i) = leaf_b(j)}   in [0, 1],
    the share of trees in which i and j share a leaf (no leaf-size weighting).
    It has a fixed, interpretable scale and is used for Forest Match Quality.
    It costs one more dense (n_treated x n_control) array.
    """
    leaves = model.apply(X)                     # (N, B)
    n, B = leaves.shape
    t_idx = np.flatnonzero(treated_mask)
    c_idx = np.flatnonzero(~treated_mask)
    P = np.zeros((len(t_idx), len(c_idx)), dtype=np.float64)
    C = np.zeros_like(P) if return_cooccurrence else None

    tree_iter = _progress(range(B), total=B, desc="proximity") if show else range(B)
    for b in tree_iter:
        col = leaves[:, b]
        _, inv = np.unique(col, return_inverse=True)
        n_leaf = inv.max() + 1
        Mt = sp.csr_matrix(
            (np.ones(len(t_idx)), (np.arange(len(t_idx)), inv[t_idx])),
            shape=(len(t_idx), n_leaf),
        )
        Mc = sp.csr_matrix(
            (np.ones(len(c_idx)), (np.arange(len(c_idx)), inv[c_idx])),
            shape=(len(c_idx), n_leaf),
        )
        size = np.asarray(Mc.sum(axis=0)).ravel()          # donors per leaf
        inv_size = np.divide(1.0, size, out=np.zeros_like(size), where=size > 0)
        P += (Mt @ sp.diags(inv_size) @ Mc.T).toarray()
        if return_cooccurrence:
            C += (Mt @ Mc.T).toarray()
    if return_cooccurrence:
        return P / B, C / B
    return P / B


def minmax_normalize(P: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Rescale the full treated x control proximity matrix to [0, 1].

    This is a GLOBAL rescaling: one shared (min, max) for every treated unit,
    not a per-treated-unit one. That is deliberate, not an oversight.
    Normalising each treated unit's row by its OWN max would send every unit's
    *best available* donor to 1.0 no matter how poor that best match actually
    is -- which would erase exactly the overlap signal the diagnostics rely on
    (a treated unit with no good donor anywhere in the control pool must still
    come out with LOW proximity, not a perfect 1.0 relative to its own weak
    options). A single global (min, max) preserves that: 1.0 is reserved for
    the single most similar treated-control pair observed anywhere in this
    run, and a unit whose best donor is mediocre gets a mediocre score.

    Returns (P_norm, raw_min, raw_max). If every pair is equally (dis)similar
    (raw_max == raw_min, a degenerate forest fit), P_norm is all zeros rather
    than dividing by zero -- the normalised proximity carries no information
    in that case, which is itself worth surfacing rather than masking.
    """
    raw_min = float(P.min())
    raw_max = float(P.max())
    spread = raw_max - raw_min
    if spread <= 0:
        return np.zeros_like(P), raw_min, raw_max
    return (P - raw_min) / spread, raw_min, raw_max


# =============================================================================
# 6. MATCHING
# =============================================================================


def selection_desc(cfg: GarfConfig) -> str:
    """One-line, human-readable description of the donor-selection rule."""
    thr = cfg.proximity_threshold
    if cfg.capped and thr > 0:
        return f"top-{int(cfg.k_matches)} nearest with normalised proximity>={thr:g} (k-NN + caliper)"
    if cfg.capped:
        return f"top-{int(cfg.k_matches)} nearest (k-NN)"
    if thr > 0:
        return f"all donors with normalised proximity>={thr:g} (caliper)"
    return "all donors that ever shared a leaf with their treated unit"


def select_donors(P: np.ndarray, P_norm: np.ndarray, treated_ids, control_ids,
                   cfg: GarfConfig) -> pd.DataFrame:
    """Pick each treated unit's donors.

    Proximity has two representations from here on:
      * raw proximity `P[i, j]` -- the leaf-co-occurrence kernel as fit; its
        magnitude depends on min_samples_leaf / n_estimators / the number of
        covariates, so it is NOT comparable across runs or specifications.
      * normalised proximity `P_norm[i, j]` -- P min-max-rescaled ONCE,
        GLOBALLY, over every treated x control pair in this run (see
        minmax_normalize()): 0 = least similar pair observed, 1 = most
        similar. The caliper is applied on THIS scale, so a given caliper
        value means the same thing regardless of forest shape.

    A donor must always have actually co-occurred with the treated unit in at
    least one tree (raw proximity > 0); a pair that never shares a leaf
    carries no information and would receive a weight of exactly 0 in any
    case, so excluding it up front is a bookkeeping/diagnostics improvement,
    not a change in what the estimator sees.

    On top of that structural floor, two composable rules:
      * caliper  -> keep only controls with NORMALISED proximity >= threshold
        (cfg.proximity_threshold = caliper; 0 = off);
      * k-NN cap -> keep at most cfg.k_matches of them (closest first -- raw
        and normalised proximity rank donors identically within a treated
        unit, since the rescaling is a common increasing affine map), unless
        k_matches is null/0, in which case ALL that pass the caliper are kept.
    With k_matches=50 and caliper=0 this reproduces the original k-NN
    behaviour (modulo dropping true zero-proximity pairs, which never
    affected estimation since their weight was always 0).
    """
    thr = cfg.proximity_threshold
    capped = cfg.capped
    kcap = int(cfg.k_matches) if capped else None
    recs = []

    if not capped and thr <= 0:
        warnings.warn(
            "No k cap and caliper=0: every control that ever shares a leaf with a "
            "treated unit becomes a donor to it (kernel matching). This can be "
            "large and slow."
        )

    if cfg.with_replacement:
        for i, tid in enumerate(treated_ids):
            cand = np.flatnonzero((P[i] > 0) & (P_norm[i] >= thr))  # co-occurred, and clears the caliper
            cand = cand[np.argsort(-P[i, cand])]              # closest first
            if capped:
                cand = cand[:kcap]                            # k-NN cap
            for r, j in enumerate(cand, start=1):
                j = int(j)
                recs.append((tid, control_ids[j], r, float(P[i, j]), float(P_norm[i, j])))
    else:
        n_t, n_c = P.shape
        kk = kcap if capped else n_c                          # effective per-treated cap
        if capped and kcap * n_t > n_c:
            warnings.warn(
                f"Without replacement, k={kcap} x {n_t} treated = {kcap * n_t} donors needed "
                f"but only {n_c} controls exist. Max feasible k = {n_c // max(n_t, 1)}. "
                f"Some treated units will be under-matched."
            )
        used, filled = set(), {i: 0 for i in range(n_t)}
        # hardest-to-match units first: increasing best available proximity
        # (beta 1.1 sorted on -P.max, i.e. easiest first, contradicting this comment)
        order_t = np.argsort(P.max(axis=1), kind="stable")
        for r in range(1, kk + 1):
            for i in order_t:
                if filled[i] >= kk:
                    continue
                cand = [(P[i, j], j) for j in range(n_c)
                        if j not in used and P[i, j] > 0 and P_norm[i, j] >= thr]
                if not cand:
                    continue
                prox, j = max(cand)
                used.add(j)
                filled[i] += 1
                recs.append((treated_ids[i], control_ids[j], r, prox, float(P_norm[i, j])))
        if capped:
            short = [treated_ids[i] for i, c in filled.items() if c < kcap]
            if short:
                warnings.warn(
                    f"Without replacement: {len(short)} treated units got fewer than "
                    f"k={kcap} donors (donor pool exhausted or below caliper)."
                )

    m = pd.DataFrame(
        recs, columns=["treated_id", "donor_id", "rank", "proximity", "proximity_minmax"]
    )
    if m.empty:
        raise RuntimeError("No donors matched. Check caliper / forest fit.")

    gamma = float(cfg.weight_power)
    if gamma == 1.0:      # exact beta 1.1 computation
        m["weight"] = m.groupby("treated_id")["proximity"].transform(lambda s: s / s.sum())
    else:                 # power kernel; gamma = 0 -> uniform (0**0 never occurs: proximity > 0)
        powered = m["proximity"] ** gamma
        m["weight"] = powered / powered.groupby(m["treated_id"]).transform("sum")
    m["leaf_distance"] = m.groupby("treated_id")["proximity"].transform(
        lambda s: 1.0 - s / s.max()
    )
    m["unique_id"] = m["donor_id"].astype(str) + "_" + m["treated_id"].astype(str)
    return m


def block_quality(P: np.ndarray, C: np.ndarray, m: pd.DataFrame,
                  treated_ids, control_ids) -> pd.DataFrame:
    """Per-treated-unit match quality; also adds m["cooccurrence"] in place.

    For treated unit i:
      overlap_mass_i     = sum_j proximity(i,j)  (share of trees where i's leaf holds a control)
      random_cooc_i      = cbar_i = mean over ALL controls of c(i,j)  (random-control benchmark)
      best_cooc_i        = Smax_i = max over ALL controls of c(i,j)   (best attainable donor)
      match_cooc_i       = S_i = sum_j w_ij c(i,j) over i's selected donors
      fmq_i              = (S_i - cbar_i) / (Smax_i - cbar_i)   (share of the attainable gain)
      fmq_twin_i         = (S_i - cbar_i) / (1 - cbar_i)        (beta 1.2-1.4 definition)
    Unmatched treated units get NaN for the match-based columns.
    """
    t_pos = {t: k for k, t in enumerate(treated_ids)}
    c_pos = {c: k for k, c in enumerate(control_ids)}
    ti = m["treated_id"].map(t_pos).to_numpy()
    cj = m["donor_id"].map(c_pos).to_numpy()
    m["cooccurrence"] = C[ti, cj]

    cbar = C.mean(axis=1)
    cmax = C.max(axis=1)
    g = m.assign(_wc=m["weight"] * m["cooccurrence"]).groupby("treated_id")
    per = pd.DataFrame({
        "n_donors": g.size(),
        "ess": g["weight"].apply(lambda w: 1.0 / (w ** 2).sum()),
        "match_cooccurrence": g["_wc"].sum(),
        "best_proximity_minmax": g["proximity_minmax"].max(),
    })
    blocks = pd.DataFrame({
        "treated_id": list(treated_ids),
        "overlap_mass": P.sum(axis=1),
        "random_cooccurrence": cbar,
        "best_cooccurrence": cmax,
    }).set_index("treated_id")
    blocks = blocks.join(per, how="left")
    blocks["matched"] = blocks["n_donors"].notna().astype(int)
    blocks["n_donors"] = blocks["n_donors"].fillna(0).astype(int)
    gain = blocks["match_cooccurrence"] - blocks["random_cooccurrence"]
    denom = blocks["best_cooccurrence"] - blocks["random_cooccurrence"]
    blocks["fmq"] = np.where(denom > 1e-12, gain / denom.where(denom > 1e-12, 1.0), np.nan)
    denom_twin = 1.0 - blocks["random_cooccurrence"]
    blocks["fmq_twin"] = np.where(denom_twin > 0, gain / denom_twin.where(denom_twin > 0, 1.0), np.nan)
    cols = ["matched", "n_donors", "ess", "match_cooccurrence", "random_cooccurrence",
            "best_cooccurrence", "fmq", "fmq_twin", "best_proximity_minmax", "overlap_mass"]
    return blocks[cols].reset_index()


def fmq_summary(blocks: pd.DataFrame) -> dict:
    """FMQ summaries. fmq_mean / fmq_p10 over MATCHED treated units; fmq_mean_all
    over ALL treated units with unmatched ones scored 0 (no better than random),
    so that a caliper cannot raise it simply by dropping hard-to-match units."""
    f = blocks.loc[blocks["matched"] == 1, "fmq"].dropna()
    f_all = blocks.loc[blocks["fmq"].notna() | (blocks["matched"] == 0), "fmq"]
    f_all = f_all.where(blocks.loc[f_all.index, "matched"] == 1, 0.0)
    ft = blocks.loc[blocks["matched"] == 1, "fmq_twin"].dropna()
    return {
        "fmq_mean": float(f.mean()) if len(f) else float("nan"),
        "fmq_p10": float(f.quantile(0.10)) if len(f) else float("nan"),
        "fmq_median": float(f.median()) if len(f) else float("nan"),
        "fmq_mean_all": float(f_all.mean()) if len(f_all) else float("nan"),
        "fmq_undefined": int(blocks.loc[blocks["matched"] == 1, "fmq"].isna().sum()),
        "fmq_twin_mean": float(ft.mean()) if len(ft) else float("nan"),
        "random_cooccurrence_mean": float(blocks["random_cooccurrence"].mean()),
        "best_cooccurrence_mean": float(blocks["best_cooccurrence"].mean()),
        "overlap_mass_mean": float(blocks["overlap_mass"].mean()),
    }


# =============================================================================
# 7. DIAGNOSTICS
# =============================================================================


def diagnostics(X, meta, m, cfg, fit_info, prox_range=None) -> tuple[dict, pd.DataFrame]:
    t_ids = meta.index[meta["W"] == 1]
    c_ids = meta.index[meta["W"] == 0]
    tm_ids = pd.Index(m["treated_id"].unique())   # matched treated units (>= 1 donor)
    w = m.groupby("donor_id")["weight"].sum()

    def smd(col):
        xt = X.loc[t_ids, col]                    # all treated: raw comparison + scale
        xt_m = X.loc[t_ids.intersection(tm_ids), col]   # matched treated: matched comparison
        xc_raw = X.loc[c_ids, col]
        wd = w.reindex(X.loc[c_ids].index).fillna(0.0).to_numpy()
        pooled = np.sqrt((xt.var() + xc_raw.var()) / 2) or 1.0   # same scale before/after
        raw = (xt.mean() - xc_raw.mean()) / pooled
        matched = (xt_m.mean() - np.average(xc_raw, weights=wd if wd.sum() else None)) / pooled
        return raw, matched

    bal = pd.DataFrame(
        [(c, *smd(c)) for c in X.columns], columns=["feature", "smd_raw", "smd_matched"]
    )
    bal["abs_raw"] = bal["smd_raw"].abs()
    bal["abs_matched"] = bal["smd_matched"].abs()

    ess = m.groupby("treated_id")["weight"].apply(lambda s: s.sum() ** 2 / (s**2).sum())

    out = {
        "n_treated": int((meta["W"] == 1).sum()),
        "n_controls": int((meta["W"] == 0).sum()),
        "n_treated_matched": int(m["treated_id"].nunique()),
        "n_unique_donors_used": int(m["donor_id"].nunique()),
        "mean_abs_smd_raw": float(bal["abs_raw"].mean()),
        "mean_abs_smd_matched": float(bal["abs_matched"].mean()),
        "max_abs_smd_raw": float(bal["abs_raw"].max()),
        "max_abs_smd_matched": float(bal["abs_matched"].max()),
        "mean_ESS_per_treated": float(ess.mean()),
        "weight_power": float(cfg.weight_power),
        "median_top_proximity": float(m[m["rank"] == 1]["proximity"].median()),
        "median_top_proximity_minmax": float(m[m["rank"] == 1]["proximity_minmax"].median()),
        "proximity_raw_min": None if prox_range is None else float(prox_range[0]),
        "proximity_raw_max": None if prox_range is None else float(prox_range[1]),
        "forest": fit_info,
    }
    return out, bal.sort_values("abs_matched", ascending=False)


# =============================================================================
# 8. DRIVER
# =============================================================================


@dataclass
class GarfResult:
    matches: pd.DataFrame
    matched_panel: pd.DataFrame
    balance: pd.DataFrame
    diagnostics: dict
    features: pd.DataFrame
    blocks: pd.DataFrame = None


def run(cfg: GarfConfig) -> GarfResult:
    t_all = time.time()
    cfg.validate()

    if not Path(cfg.data_path).exists():
        raise FileNotFoundError(
            f"data_path '{cfg.data_path}' not found. Put your panel CSV there "
            f"(default: the 'input' folder next to config.yaml), or edit data_path."
        )

    _log(f"loading data from {cfg.data_path}")
    df = pd.read_csv(cfg.data_path)
    required = [cfg.id_col, cfg.year_col, cfg.treat_col, cfg.outcome_col,
                *cfg.x_time_varying, *cfg.x_time_invariant]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(
            f"Columns not found in {cfg.data_path}: {missing}. "
            f"Available columns: {list(df.columns)}"
        )
    _log(f"  {len(df):,} rows, {df[cfg.id_col].nunique():,} units, "
         f"years {int(df[cfg.year_col].min())}-{int(df[cfg.year_col].max())}")

    t = time.time()
    _log("building pre-treatment features")
    X, meta = build_features(df, cfg)
    treated_mask = meta["W"].to_numpy() == 1
    if treated_mask.sum() == 0 or (~treated_mask).sum() == 0:
        raise RuntimeError("Need both treated and control units with a complete pre-window.")
    _log(f"  {X.shape[0]:,} units x {X.shape[1]} features "
         f"({int(treated_mask.sum())} treated, {int((~treated_mask).sum())} control) [{_dt(t)}]")

    t = time.time()
    _log(f"fitting random forest  (tuning: {'ON — this is the slow path' if cfg.tune else 'off'})")
    model, fit_info = fit_forest(X, meta, cfg)
    _log(f"  forest ready [{_dt(t)}]  OOB acc={fit_info.get('oob_score', float('nan')):.3f}  "
         f"bal.acc={fit_info.get('oob_balanced_accuracy', float('nan')):.3f}  "
         f"AUC={fit_info.get('oob_auc', float('nan')):.3f} (out-of-bag)"
         + (f"  CV({cfg.tuning_metric})={fit_info['cv_score']:.3f}" if cfg.tune else ""))
    if cfg.tune:
        _log(f"  best params: {fit_info['best_params']}")

    t = time.time()
    _log(f"computing leaf proximity over {model.n_estimators} trees")
    P, C = leaf_proximity(model, X, treated_mask, show=cfg.progress,
                          return_cooccurrence=True)
    _log(f"  proximity matrix {P.shape[0]}x{P.shape[1]} [{_dt(t)}]")

    P_norm, prox_min, prox_max = minmax_normalize(P)
    _log(f"  raw proximity range [{prox_min:.4g}, {prox_max:.4g}] -> "
         f"min-max normalised to [0, 1] (caliper is applied on this scale)")

    t = time.time()
    _log(f"selecting donors: {selection_desc(cfg)} "
         f"({'with' if cfg.with_replacement else 'without'} replacement; "
         f"weights ~ proximity^{cfg.weight_power:g})")
    t_list = X.index[treated_mask].tolist()
    c_list = X.index[~treated_mask].tolist()
    m = select_donors(P, P_norm, t_list, c_list, cfg)
    blocks = block_quality(P, C, m, t_list, c_list)
    diag, bal = diagnostics(X, meta, m, cfg, fit_info, prox_range=(prox_min, prox_max))
    diag.update(fmq_summary(blocks))
    _log(f"  {diag['n_treated_matched']} treated matched to "
         f"{diag['n_unique_donors_used']} distinct donors [{_dt(t)}]")

    # ---- matched panel (long, ready for DiD) --------------------------------
    _log("assembling matched panel and writing outputs")
    keep = [cfg.id_col, cfg.year_col, cfg.treat_col, cfg.outcome_col,
            *cfg.x_time_varying, *cfg.x_time_invariant]
    long = df[keep].copy()

    donors = m.merge(long, left_on="donor_id", right_on=cfg.id_col, how="left")
    donors["block_id"] = donors["treated_id"]
    donors["role"] = "donor"

    tre = long[long[cfg.id_col].isin(m["treated_id"].unique())].copy()
    tre["treated_id"] = tre[cfg.id_col]
    tre["donor_id"] = tre[cfg.id_col]
    tre["unique_id"] = tre[cfg.id_col].astype(str) + "_" + tre[cfg.id_col].astype(str)
    tre["block_id"] = tre[cfg.id_col]
    tre["rank"] = 0
    tre["proximity"] = np.nan
    tre["proximity_minmax"] = np.nan
    tre["weight"] = 1.0
    tre["leaf_distance"] = 0.0
    tre["cooccurrence"] = np.nan
    tre["role"] = "treated"

    cols = ["unique_id", cfg.id_col, "treated_id", "donor_id", "block_id", "role",
            cfg.year_col, cfg.treat_col, cfg.outcome_col,
            *cfg.x_time_varying, *cfg.x_time_invariant,
            "rank", "proximity", "proximity_minmax", "weight", "leaf_distance", "cooccurrence"]
    panel = pd.concat([tre[cols], donors[cols]], ignore_index=True)
    panel["fmq"] = panel["block_id"].map(blocks.set_index("treated_id")["fmq"])
    panel = panel.sort_values(["block_id", "role", "unique_id", cfg.year_col])

    out = Path(cfg.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    m.to_csv(out / "matches.csv", index=False)
    panel.to_csv(out / "matched_panel.csv", index=False)
    bal.to_csv(out / "balance.csv", index=False)
    blocks.to_csv(out / "blocks.csv", index=False)
    X.to_csv(out / "features.csv")
    with open(out / "diagnostics.json", "w") as fh:
        json.dump({"config": asdict(cfg), "diagnostics": diag}, fh, indent=2, default=str)

    kdisp = cfg.k_matches if cfg.capped else "all"
    print()
    _log(f"DONE in {_dt(t_all)}  |  k={kdisp}  caliper={cfg.proximity_threshold:g} (normalised)  "
         f"wpower={cfg.weight_power:g}  replacement={cfg.with_replacement}  "
         f"tuning={'on' if cfg.tune else 'off'}")
    _log(f"treated={diag['n_treated']} controls={diag['n_controls']} "
         f"donors used={diag['n_unique_donors_used']}")
    _log(f"mean |SMD|: {diag['mean_abs_smd_raw']:.3f} -> {diag['mean_abs_smd_matched']:.3f}   "
         f"max |SMD|: {diag['max_abs_smd_raw']:.3f} -> {diag['max_abs_smd_matched']:.3f}")
    _log(f"mean ESS per treated: {diag['mean_ESS_per_treated']:.1f}"
         + (f" / {kdisp}" if cfg.capped else ""))
    _log(f"forest match quality FMQ: mean {diag['fmq_mean']:.3f}  p10 {diag['fmq_p10']:.3f}  "
         f"all-treated {diag['fmq_mean_all']:.3f}   (0 = random controls, 1 = best attainable donors; "
         f"read with ESS: fewer/closer donors systematically raise it)")
    _log(f"  twin-scaled FMQ (beta 1.4 definition) mean {diag['fmq_twin_mean']:.3f}; "
         f"co-occurrence: donors {blocks['match_cooccurrence'].mean():.3f}, "
         f"best attainable {diag['best_cooccurrence_mean']:.3f}, random {diag['random_cooccurrence_mean']:.3f}"
         + (f"; FMQ undefined for {diag['fmq_undefined']} block(s)" if diag['fmq_undefined'] else ""))
    _log(f"results written to {out.resolve()}")
    print()

    return GarfResult(m, panel, bal, diag, X, blocks)


TEMPLATE = """\
# =====================================================================
# GARF (binary) config.  Edit the values, then run the tool.
# Paths are relative to this 'engine' folder.
# =====================================================================

# ---- data -----------------------------------------------------------
data_path: input/panel.csv   # put your panel CSV in the 'input' folder
id_col: id                   # unit identifier column
year_col: year               # time column
treat_col: W                 # binary treatment column (0/1)
outcome_col: Y               # outcome column
x_time_varying: [x1, x2]     # covariates that change over time (get lags + trajectory features)
x_time_invariant: [region]   # covariates that don't change (strings are fine, no encoding needed)

# ---- design ---------------------------------------------------------
treatment_year: 2000         # T0: first treated period
n_pre_years: 5               # pre-window length used for matching: [T0 - n, T0 - 1]
include_treatment_year: false

# ---- matching -------------------------------------------------------
# Two ways to choose each treated unit's controls; they COMBINE:
k_matches: 50                # keep at most this many nearest donors (k-NN); null or 0 = no cap
caliper: 0.0                 # keep only donors with MIN-MAX NORMALISED proximity >= caliper
                              # (range [0,1]; 0 = off). See proximity_minmax in the outputs.
with_replacement: true       # true = donors reusable across treated units (recommended)
weight_power: 1              # donor weights ~ proximity^weight_power (1 = proportional,
                             # 0 = uniform k:1, >1 = more weight on the closest donors)

# ---- hyperparameter tuning (model selection) ------------------------
tune: false                  # false = fast (uses the fixed forest below)
                             # true  = SLOW: cross-validated randomized search
tuning_metric: combined      # accuracy | f1_treated | combined
n_iter_search: 40            # candidates tried when tuning
cv_folds: 5                  # cross-validation folds when tuning

# ---- forest defaults (used when tune: false) ------------------------
n_estimators: 1000
min_samples_leaf: 5
max_features: sqrt           # sqrt | log2 | a fraction like 0.5
max_depth: null              # null = unlimited
random_state: 42

# ---- output ---------------------------------------------------------
out_dir: output              # results are written into the 'output' folder
progress: true               # live progress + timing in the terminal
"""


def main():
    ap = argparse.ArgumentParser(description="GARF (binary) forest-proximity matching")
    ap.add_argument("--config", help="path to config.yaml")
    ap.add_argument("--init", action="store_true", help="write a template config.yaml")
    a = ap.parse_args()
    if a.init:
        Path("config.yaml").write_text(TEMPLATE)
        print("wrote config.yaml — edit it, then: python garf_binary.py --config config.yaml")
        return
    if not a.config:
        ap.error("give --config config.yaml (or --init to create one)")
    run(GarfConfig.from_yaml(a.config))


if __name__ == "__main__":
    main()
