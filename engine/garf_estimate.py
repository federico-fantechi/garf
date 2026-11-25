#!/usr/bin/env python3
"""
garf_estimate.py -- estimation stage for GARF (binary), version 1.2
===================================================================

Federico Fantechi, PhD -- University of Palermo -- federico.fantechi@unipa.it

Takes the matched panel produced by garf_binary.py (one block per treated unit:
the treated unit's rows plus its donors' rows, donor weights w_ij summing to 1
within the block) and estimates the effect of the treatment on the treated.
It is a separate step: the matching engine never calls it. Run it

  * from Stata, through the command  garfte  (after garf), or
  * from a terminal, on a saved garf output folder:
        python garf_estimate.py --from output --method did
        python garf_estimate.py --from output --method gap

The gap
-------
For treated unit i (block i) and period t, the GARF counterfactual is the
weighted mean of i's donors, and the gap is

    G_it = Y_it - sum_j w~_ij(t) Y_jt

w~_ij(t) are the design weights renormalised over the donors observed in t
(identical to w_ij when no donor is missing; the number of donors observed in
each cell is reported).

Methods (no default: choose one)
--------------------------------
  did   matching difference-in-differences. tau_i(e) = G_i,T0+e - Gbar_i, where
        Gbar_i is the mean gap over the matching window [T0-L, T0-1] (or
        [T0-L+1, T0] with includet0). Identifying assumption: parallel trends
        conditional on the match. Requires a panel matched WITHOUT replacement
        (no donor in more than one block); standard errors are clustered by block
        (Abadie and Spiess 2022, JASA). Same estimand as Callaway and Sant'Anna's
        ATT(g,t) with a single cohort, and as trajectory balancing with demeaning
        (Hazlett and Xu 2018), with GARF weights in place of IPW/balancing weights.
  gap   the gap itself, in levels: tau_i(e) = G_i,T0+e. Identifying assumption:
        unconfoundedness given the matched pre-treatment history (the logic of
        matching on lagged outcomes and of the synthetic control's post-treatment
        gap). Replacement allowed. Standard errors: cluster by block plus a
        correction for donors shared between blocks (see below).
        Under fairly general conditions the did and gap estimates bracket the
        effect (Ding and Li 2019, Political Analysis): report both.

Estimates
---------
  ATT(e)  = mean over treated units with a defined tau_i(e)
  ATT     = equal-weight average of ATT(e) over the post-treatment event times
            shown (Callaway-Sant'Anna "dynamic" aggregation)

Inference (conditional on the matched design; forest uncertainty is ignored, as
in propensity-score matching)
  V = G/(G-1) * sum_i d_i d_i'        d_i(e) = (tau_i(e) - ATT(e)) / N_e
    + sum_j  M_j o Sigma_j            only when donors are shared (gap)
  M_j(e,f)    = A_j(e) A_j(f) - sum_i a_ij(e) a_ij(f),   a_ij(e) = w~_ij(e)/N_e,
                A_j(e) = sum_i a_ij(e)   (cross-block reuse of donor j)
  Sigma_j     = weighted mean, over the blocks donor j belongs to, of the block's
                weighted covariance of donor outcomes around the counterfactual.
  Without shared donors the second term is zero and V is the block-clustered
  variance. Critical values: Student t with (number of treated units - 1) df.

Outputs (prefix garfte_<method>_), written to --outdir
  event.csv     e, year, period, att, se, t, p, lb, ub, n_treated, mean_donors_obs
  units.csv     treated_id, year, e, period, tau, pre_gap_mean, n_donors_obs,
                n_donors, ess, fmq          (one row per treated unit x period)
  vcov.csv      full covariance matrix of (ATT(e) ..., ATT)
  summary.csv   key,value pairs (read by garfte)
  messages.txt  notes and warnings (printed by garfte)
  error.txt     only if the estimation stopped; the reason
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = "1.2"

METHODS = ("did", "gap")
WEIGHT_TOL = 1e-4          # tolerance on "donor weights sum to 1 in each block"
                           # (Stata stores imported decimals as float: ~7 digits)
FEW_TREATED = 30


class GarfteError(Exception):
    """An error the user can act on; its message is shown as is."""


# =============================================================================
# helpers
# =============================================================================


def ename(e: int) -> str:
    """Stata-valid coefficient name for event time e: em3, e0, e2."""
    return f"em{-e}" if e < 0 else f"e{e}"


def _tcrit(level: float, df: int) -> float:
    from scipy import stats
    return float(stats.t.ppf(0.5 + level / 200.0, df))


def _pval(t: float, df: int) -> float:
    from scipy import stats
    if not np.isfinite(t):
        return float("nan")
    return float(2 * stats.t.sf(abs(t), df))


# =============================================================================
# core
# =============================================================================


def estimate(panel: pd.DataFrame, *, method: str, id_col: str, time_col: str,
             t0: int, npre: int, includet0: bool, matched_outcome: str,
             outcome: str | None = None, pre: int | None = None,
             post: int | None = None, level: float = 95.0) -> dict:
    """Estimate ATT(e) and ATT on a GARF matched panel. Returns a dict of results."""
    method = str(method).lower()
    if method not in METHODS:
        raise GarfteError(f"method must be one of {', '.join(METHODS)} (got '{method}').")
    outcome = outcome or matched_outcome
    if not (0 < level < 100):
        raise GarfteError("level must be between 0 and 100.")
    msgs: list[str] = []

    need = [id_col, time_col, "block_id", "role", "donor_id", "weight", outcome]
    missing = [c for c in need if c not in panel.columns]
    if missing:
        raise GarfteError(
            f"column(s) {missing} not found: this does not look like a GARF matched panel "
            f"(or outcome '{outcome}' is not in it).")

    df = panel[[c for c in dict.fromkeys(need + (["fmq"] if "fmq" in panel.columns else []))]].copy()
    df["role"] = df["role"].astype(str).str.strip()
    if not set(df["role"].unique()) <= {"treated", "donor"}:
        raise GarfteError("column 'role' must contain only 'treated' and 'donor'.")
    for c in ("block_id", "donor_id"):
        df[c] = df[c].astype(str)
    df[time_col] = pd.to_numeric(df[time_col], errors="raise").astype(int)
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
    df["_y"] = pd.to_numeric(df[outcome], errors="coerce")

    tre = df[df["role"] == "treated"]
    don = df[df["role"] == "donor"]
    if tre.empty:
        raise GarfteError("no treated rows in the matched panel.")

    # ---- integrity: one treated unit per block, donor weights sum to 1 ---------
    if tre.groupby("block_id")["donor_id"].nunique().max() > 1:
        raise GarfteError("matched panel was modified: a block contains more than one treated unit.")
    if tre.duplicated(["block_id", time_col]).any() or don.duplicated(["block_id", "donor_id", time_col]).any():
        raise GarfteError("matched panel was modified: duplicated block x unit x period rows.")
    pair_w = don.groupby(["block_id", "donor_id"])["weight"].agg(["min", "max"])
    if pair_w.isna().any().any() or ((pair_w["max"] - pair_w["min"]).abs() > WEIGHT_TOL).any():
        raise GarfteError("matched panel was modified: a donor's weight differs across its rows "
                          "(or is missing).")
    wsum = pair_w["max"].groupby(level=0).sum()
    blocks_all = pd.Index(sorted(tre["block_id"].unique()))
    no_donor = blocks_all.difference(wsum.index)
    bad = wsum[(wsum - 1).abs() > WEIGHT_TOL]
    if len(bad):
        raise GarfteError(
            f"matched panel was modified: donor weights do not sum to 1 in {len(bad)} block(s) "
            f"(e.g. block {bad.index[0]}: {bad.iloc[0]:.6f}). Were donors dropped? garfte needs the "
            f"complete matched panel: rerun garf (or reload its output folder).")
    if len(no_donor):
        msgs.append(f"note: {len(no_donor)} treated unit(s) have no donor (e.g. no control within "
                    f"the caliper) and are left out.")

    # ---- donor reuse ------------------------------------------------------------
    reuse = don.groupby("donor_id")["block_id"].nunique()
    n_reused = int((reuse > 1).sum())
    if method == "did" and n_reused:
        raise GarfteError(
            f"method(did) needs a panel matched WITHOUT replacement, but {n_reused} donor(s) serve "
            f"more than one treated unit (max {int(reuse.max())}). Rerun garf with noreplacement and a "
            f"k small enough for the donor pool (the pool must hold about k x the number of treated "
            f"units), or use method(gap), which allows replacement.")

    # ---- window and event times -------------------------------------------------
    last = t0 if includet0 else t0 - 1
    first = last - npre + 1
    win_years = set(range(first, last + 1))

    # ---- gap G_it -----------------------------------------------------------------
    w_pair = pair_w["max"].rename("w")
    don = don.join(w_pair, on=["block_id", "donor_id"])
    don_obs = don[don["_y"].notna()].copy()
    don_obs["_wy"] = don_obs["w"] * don_obs["_y"]
    cf = don_obs.groupby(["block_id", time_col]).agg(sw=("w", "sum"), swy=("_wy", "sum"),
                                                    n_obs=("w", "size"))
    cf["ycf"] = cf["swy"] / cf["sw"]
    n_donors = don.groupby("block_id")["donor_id"].nunique()
    ess = w_pair.groupby(level=0).apply(lambda s: 1.0 / (s ** 2).sum())

    T = tre[["block_id", "donor_id", time_col, "_y"] + (["fmq"] if "fmq" in tre.columns else [])]
    T = T.rename(columns={"donor_id": "treated_id"}).join(cf[["ycf", "n_obs", "sw"]],
                                                         on=["block_id", time_col])
    T = T[T["block_id"].isin(wsum.index)]
    T["G"] = T["_y"] - T["ycf"]
    n_renorm = int(((T["sw"] < 1 - WEIGHT_TOL) & T["G"].notna()).sum())
    if n_renorm:
        msgs.append(f"note: in {n_renorm} treated unit x period cell(s) some donors are missing; "
                    f"weights were renormalised over the donors observed (see n_donors_obs).")
    n_undef = int(T["G"].isna().sum())
    if n_undef:
        msgs.append(f"note: {n_undef} treated unit x period cell(s) have no gap (outcome missing "
                    f"for the treated unit or for all its donors) and are skipped.")

    base = T[T[time_col].isin(win_years)].groupby("block_id")["G"].mean()
    T["pre_gap_mean"] = T["block_id"].map(base)
    win_cover = T[T[time_col].isin(win_years)].groupby("block_id")["G"].count()
    if outcome == matched_outcome and (win_cover < npre).any():
        msgs.append("note: some treated units have gaps for only part of the matching window.")
    no_base = blocks_all.intersection(wsum.index).difference(base.dropna().index)
    if len(no_base):
        msgs.append(f"note: {len(no_base)} treated unit(s) have no observed gap in the matching "
                    f"window and are left out.")
        T = T[~T["block_id"].isin(no_base)]

    T["tau"] = T["G"] - T["pre_gap_mean"] if method == "did" else T["G"]
    T["e"] = T[time_col] - t0

    def period_of(y: int) -> str:
        if y > last:
            return "post"
        if y in win_years:
            return "window" if outcome == matched_outcome else "pre (window)"
        return "pre"
    T["period"] = T[time_col].map(period_of)

    lo = -pre if pre is not None and pre >= 0 else -10**9
    hi = post if post is not None and post >= 0 else 10**9
    T = T[(T["e"] >= lo) & (T["e"] <= hi)]
    Tv = T[T["tau"].notna()]
    if Tv.empty:
        raise GarfteError("no estimable treated unit x period cell in the requested window.")
    events = sorted(int(e) for e in Tv["e"].unique())
    skipped = sorted(set(int(e) for e in T["e"].unique()) - set(events))
    if skipped:
        msgs.append(f"note: event time(s) {skipped} have no estimable cell and are skipped.")
    post_events = [e for e in events if e + t0 > last]
    if not post_events:
        raise GarfteError("no post-treatment period in the data (or in the post() window).")
    E = len(events)
    eidx = {e: k for k, e in enumerate(events)}

    blocks = sorted(Tv["block_id"].unique())
    bidx = {b: k for k, b in enumerate(blocks)}
    G = len(blocks)
    if G < 2:
        raise GarfteError("at least two treated units are needed for standard errors.")

    tau = np.full((G, E), np.nan)
    tau[Tv["block_id"].map(bidx).to_numpy(), Tv["e"].map(eidx).to_numpy()] = Tv["tau"].to_numpy()
    defined = ~np.isnan(tau)
    N_e = defined.sum(axis=0)
    att_e = np.nanmean(tau, axis=0)

    # ---- variance: block-clustered term ----------------------------------------------
    d = np.where(defined, (tau - att_e) / N_e, 0.0)
    V = (G / (G - 1)) * d.T @ d

    # ---- variance: donor-reuse correction (gap with replacement) ---------------------
    reuse_term = False
    if method == "gap" and n_reused:
        reuse_term = True
        V = V + _reuse_correction(don_obs, cf, Tv, events, eidx, N_e, time_col, t0)

    # ---- aggregate ATT over post event times (equal weights) -------------------------
    a = np.zeros(E)
    for e in post_events:
        a[eidx[e]] = 1.0 / len(post_events)
    names = [ename(e) for e in events] + ["ATT"]
    b = np.append(att_e, a @ att_e)
    A = np.vstack([np.eye(E), a])
    Vfull = A @ V @ A.T
    Vfull = (Vfull + Vfull.T) / 2
    ev = np.linalg.eigvalsh(Vfull)
    if ev.min() < -1e-10 * max(ev.max(), 1e-300):
        msgs.append("note: the variance matrix had small negative eigenvalues (reuse correction); "
                    "they were set to zero.")
        w_, Q = np.linalg.eigh(Vfull)
        Vfull = (Q * np.clip(w_, 0, None)) @ Q.T

    df_t = G - 1
    tc = _tcrit(level, df_t)
    se = np.sqrt(np.clip(np.diag(Vfull), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = np.where(se > 0, b / se, np.nan)
    pv = np.array([_pval(t, df_t) for t in tstat])

    # ---- tables --------------------------------------------------------------------------
    mean_obs = Tv.groupby("e")["n_obs"].mean()
    event = pd.DataFrame({
        "name": names,
        "e": events + [np.nan],
        "year": [e + t0 for e in events] + [np.nan],
        "period": [period_of(e + t0) for e in events] + ["overall post"],
        "att": b, "se": se, "t": tstat, "p": pv,
        "lb": b - tc * se, "ub": b + tc * se,
        "n_treated": list(N_e) + [int(defined[:, [eidx[e] for e in post_events]].any(axis=1).sum())],
        "mean_donors_obs": [mean_obs.get(e, np.nan) for e in events] + [np.nan],
    })
    units = T.assign(n_donors=T["block_id"].map(n_donors), ess=T["block_id"].map(ess))
    units = units.rename(columns={"n_obs": "n_donors_obs", time_col: "year"})
    ucols = ["treated_id", "block_id", "year", "e", "period", "tau", "pre_gap_mean",
             "n_donors_obs", "n_donors", "ess"] + (["fmq"] if "fmq" in units.columns else [])
    units = units[ucols].sort_values(["treated_id", "year"]).reset_index(drop=True)
    vcov = pd.DataFrame(Vfull, index=names, columns=names)

    if G < FEW_TREATED:
        msgs.append(f"warning: only {G} treated units. Standard errors that rely on many treated "
                    f"units (clusters) can be unreliable with fewer than about {FEW_TREATED}; "
                    f"read the confidence intervals with care.")
    if method == "did":
        msgs.append("window periods: tau is measured from the window mean, so the window "
                    "estimates sum to zero by construction. They are NOT a pre-trends test; "
                    "use periods before the window (if any) or the placebo in DIAGNOSTICS/06.")
    elif outcome == matched_outcome:
        msgs.append("window periods: the outcome was matched on in the window, so gaps there "
                    "are small by construction and are NOT a pre-trends test. A large mean "
                    "window gap (pre_gap_mean) means the level comparison is not credible.")
    if reuse_term:
        msgs.append(f"donor reuse: {n_reused} donor(s) serve more than one treated unit; the "
                    f"standard errors include the correction for shared donors. Regressions of "
                    f"the unit-level effects (units.csv) on covariates need the same care: "
                    f"plain OLS standard errors ignore the shared donors.")
    if outcome != matched_outcome:
        msgs.append(f"outcome '{outcome}' was not matched on (garf matched on '{matched_outcome}'): "
                    f"window estimates are informative about pre-treatment differences.")

    k = names.index("ATT")
    summary = {
        "version": __version__, "method": method, "outcome": outcome,
        "matched_outcome": matched_outcome, "t0": t0, "npre": npre,
        "includet0": int(includet0), "window_first": first, "window_last": last,
        "n_treated": G, "n_donors": int(don["donor_id"].nunique()),
        "donors_reused": n_reused, "reuse_correction": int(reuse_term),
        "n_events": E, "n_post": len(post_events),
        "e_min": events[0], "e_max": events[-1],
        "att": b[k], "se": se[k], "t": tstat[k], "p": pv[k], "lb": b[k] - tc * se[k],
        "ub": b[k] + tc * se[k], "df": df_t, "level": level,
        "mean_window_gap": float(base.reindex(blocks).mean()),
        "mean_abs_window_gap": float(base.reindex(blocks).abs().mean()),
        "renormalised_cells": n_renorm,
    }
    return {"event": event, "units": units, "vcov": vcov, "summary": summary, "messages": msgs}


def _reuse_correction(don_obs, cf, Tv, events, eidx, N_e, time_col, t0) -> np.ndarray:
    """sum_j M_j o Sigma_j for donors shared between blocks (gap method)."""
    E = len(events)
    # cells used in the estimate (block, period) and their event index
    cells = Tv[["block_id", time_col]].assign(k=Tv["e"].map(eidx).to_numpy())
    D = don_obs.merge(cells, on=["block_id", time_col], how="inner")
    D = D.join(cf[["sw", "ycf"]], on=["block_id", time_col])
    D["wt"] = D["w"] / D["sw"]                       # renormalised weight w~_ij(t)
    D["a"] = D["wt"] / N_e[D["k"].to_numpy()]
    D["dev"] = D["_y"] - D["ycf"]

    # pair-level vectors a_ij(.) and deviations dev_ij(.) over events
    pairs = D[["block_id", "donor_id"]].drop_duplicates().reset_index(drop=True)
    pidx = pd.Series(np.arange(len(pairs)), index=pd.MultiIndex.from_frame(pairs))
    p = pidx.loc[list(zip(D["block_id"], D["donor_id"]))].to_numpy()
    P = len(pairs)
    a = np.zeros((P, E)); a[p, D["k"].to_numpy()] = D["a"].to_numpy()
    dev = np.zeros((P, E)); obs = np.zeros((P, E))
    dev[p, D["k"].to_numpy()] = D["dev"].to_numpy(); obs[p, D["k"].to_numpy()] = 1.0
    wt = np.zeros((P, E)); wt[p, D["k"].to_numpy()] = D["wt"].to_numpy()

    # block covariance S_b(e,f) = sum_j w~ dev_e dev_f / (1 - sum_j w~^2)   (both observed)
    blocks = pairs["block_id"].astype("category")
    bcode = blocks.cat.codes.to_numpy(); B = len(blocks.cat.categories)
    wbar = np.sqrt(wt[:, :, None] * wt[:, None, :]) * obs[:, :, None] * obs[:, None, :]
    num = np.zeros((B, E, E)); w2 = np.zeros((B, E, E)); sw = np.zeros((B, E, E))
    np.add.at(num, bcode, wbar * dev[:, :, None] * dev[:, None, :])
    np.add.at(w2, bcode, wbar ** 2)
    np.add.at(sw, bcode, wbar)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(sw > 0, (num / np.where(sw > 0, sw, 1)) / np.clip(1 - w2 / np.where(sw > 0, sw, 1) ** 2, 1e-12, None), 0.0)

    # donor covariance Sigma_j = weighted mean over j's blocks of S_b (weights w_ij)
    donors = pairs["donor_id"].astype("category")
    dcode = donors.cat.codes.to_numpy(); J = len(donors.cat.categories)
    wpair = D.groupby(["block_id", "donor_id"])["w"].first().reindex(
        pd.MultiIndex.from_frame(pairs)).to_numpy()
    Sig = np.zeros((J, E, E)); wj = np.zeros(J)
    np.add.at(Sig, dcode, wpair[:, None, None] * S[bcode])
    np.add.at(wj, dcode, wpair)
    Sig /= np.where(wj > 0, wj, 1)[:, None, None]

    # M_j = A_j A_j' - sum_i a_ij a_ij'
    A = np.zeros((J, E)); np.add.at(A, dcode, a)
    corr = np.einsum("je,jf,jef->ef", A, A, Sig)
    for s in range(0, P, 20000):          # chunked: P x E x E can be large
        sl = slice(s, s + 20000)
        corr -= np.einsum("pe,pf,pef->ef", a[sl], a[sl], Sig[dcode[sl]])
    return corr


# =============================================================================
# I/O
# =============================================================================


def write_outputs(res: dict, outdir: Path, method: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    pre = outdir / f"garfte_{method}_"
    res["event"].to_csv(f"{pre}event.csv", index=False)
    res["units"].to_csv(f"{pre}units.csv", index=False)
    res["vcov"].to_csv(f"{pre}vcov.csv", index_label="name")
    with open(f"{pre}summary.csv", "w") as fh:
        fh.write("key,value\n")
        for k, v in res["summary"].items():
            if isinstance(v, (bool, np.bool_)):
                v = int(v)
            elif isinstance(v, (int, np.integer)):
                v = int(v)
            elif isinstance(v, (float, np.floating)):
                v = repr(float(v)) if np.isfinite(v) else "."
            fh.write(f"{k},{v}\n")
    import textwrap
    with open(f"{pre}messages.txt", "w") as fh:
        for m in res["messages"]:
            # wrapped for the Stata log; continuation lines start with two spaces
            for line in textwrap.wrap(m, width=76, subsequent_indent="  ", break_on_hyphens=False):
                fh.write(line + "\n")


def _clear(outdir: Path, method: str) -> None:
    for s in ("event.csv", "units.csv", "vcov.csv", "summary.csv", "messages.txt", "error.txt"):
        f = outdir / f"garfte_{method}_{s}"
        if f.exists():
            f.unlink()


def _print(res: dict) -> None:
    s = res["summary"]
    print(f"\nGARF estimation -- method {s['method']}, outcome {s['outcome']}, "
          f"{s['n_treated']} treated units, T0 = {s['t0']}")
    ev = res["event"]
    print(ev[["name", "period", "att", "se", "lb", "ub", "n_treated"]].to_string(index=False,
          float_format=lambda x: f"{x:9.4f}"))
    for m in res["messages"]:
        print(m)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="GARF estimation stage (did | gap)")
    ap.add_argument("--method", required=True, choices=METHODS)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--from", dest="from_dir", help="a garf output folder "
                     "(matched_panel.csv + diagnostics.json)")
    src.add_argument("--panel", help="a matched-panel CSV (settings given as options)")
    ap.add_argument("--t0", type=int); ap.add_argument("--npre", type=int)
    ap.add_argument("--includet0", type=int, default=0)
    ap.add_argument("--id"); ap.add_argument("--time"); ap.add_argument("--matched-outcome")
    ap.add_argument("--outcome", help="outcome to estimate on (default: the matched one)")
    ap.add_argument("--pre", type=int, default=-1, help="keep event times >= -pre (-1 = all)")
    ap.add_argument("--post", type=int, default=-1, help="keep event times <= post (-1 = all)")
    ap.add_argument("--level", type=float, default=95.0)
    ap.add_argument("--outdir", help="default: <from>/estimations")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    outdir = Path(a.outdir) if a.outdir else (
        Path(a.from_dir) / "estimations" if a.from_dir else Path("garf_output/estimations"))
    outdir.mkdir(parents=True, exist_ok=True)
    _clear(outdir, a.method)
    try:
        if a.from_dir:
            src_dir = Path(a.from_dir)
            mp, dj = src_dir / "matched_panel.csv", src_dir / "diagnostics.json"
            for f in (mp, dj):
                if not f.exists():
                    raise GarfteError(f"{f} not found: is '{src_dir}' a garf output folder "
                                      f"(saved with store / export() / outdir())?")
            cfg = json.loads(dj.read_text())["config"]
            settings = dict(id_col=cfg["id_col"], time_col=cfg["year_col"],
                            t0=int(cfg["treatment_year"]), npre=int(cfg["n_pre_years"]),
                            includet0=bool(cfg["include_treatment_year"]),
                            matched_outcome=cfg["outcome_col"])
            panel = pd.read_csv(mp)
        else:
            for opt in ("t0", "npre", "id", "time", "matched_outcome"):
                if getattr(a, opt) is None:
                    raise GarfteError(f"--panel needs --{opt.replace('_', '-')}")
            settings = dict(id_col=a.id, time_col=a.time, t0=a.t0, npre=a.npre,
                            includet0=bool(a.includet0), matched_outcome=a.matched_outcome)
            panel = pd.read_csv(a.panel)
        res = estimate(panel, method=a.method, outcome=a.outcome,
                       pre=a.pre if a.pre >= 0 else None, post=a.post if a.post >= 0 else None,
                       level=a.level, **settings)
    except GarfteError as exc:
        (outdir / f"garfte_{a.method}_error.txt").write_text(str(exc) + "\n")
        print(f"garfte error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # unexpected: keep the traceback for the user
        import traceback
        (outdir / f"garfte_{a.method}_error.txt").write_text(
            f"unexpected error in garf_estimate.py: {exc!r}\n" + traceback.format_exc())
        raise
    write_outputs(res, outdir, a.method)
    if not a.quiet:
        _print(res)
        print(f"\nresults written to {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
