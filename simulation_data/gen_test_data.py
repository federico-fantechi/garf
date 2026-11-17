"""
gen_test_data.py - simulated 50,000 x 15 panel with three treatment-assignment scenarios
(random, selection on observables incl. trajectories, selection on unobservables) for testing
GARF, PSM, DiD, etc.

"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats, optimize

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
rng = np.random.default_rng(42)
N, T, T0, NTREAT = 50_000, 15, 9, 500          # units, years, first treated year, treated units
ids = np.arange(1, N + 1)
years = np.arange(1, T + 1)
truth = []                                       # rows: (parameter, value)

# ---------------------------------------------------------------- fixed effects
unit_effect = rng.normal(0, 1, N)
year_effect = rng.normal(0, 0.5, T) + 0.05 * years
# DGP NOTES (v2) ----------------------------------------------------------------------
# g_i: unit-specific trend, centred on the middle of the selection window (year 6), so that
#      it changes the SLOPE of x_1 and y0 over years 4-8 but not their window MEANS.
# h_i: unobserved shock that enters y0 only from year 9 (never visible pre-treatment).
G_X1, G_Y, H_Y = 0.15, 0.08, 0.30                # loadings of g on x_1 and y0, of h on y0
g = rng.normal(0, 1, N)
h = rng.normal(0, 1, N)
WIN = 5                                          # selection window: years T0-WIN .. T0-1
tc = years - (T0 - (WIN + 1) / 2)                # centred on the window midpoint (year 6)
truth.append(("trend_g", f"g_i ~ N(0,1); x_1 latent += {G_X1}*g_i*(t-6); y0 += {G_Y}*g_i*(t-6)"))
truth.append(("shock_h", f"h_i ~ N(0,1); y0 += {H_Y}*h_i*1[t>=9] (unobserved)"))
for t in years:
    truth.append((f"year_effect_{t}", year_effect[t - 1]))

# ---------------------------------------------------------------- x_1..x_7 (time-varying)
def latent(rng):
    """unit + year + noise component with a variable-specific mix; returns N x T standard-normal-ish array."""
    w = rng.dirichlet([2, 1, 2])                 # shares: unit, year, noise
    u = rng.normal(0, 1, N)[:, None]
    v = rng.normal(0, 1, T)[None, :]
    e = rng.normal(0, 1, (N, T))
    return np.sqrt(w[0]) * u + np.sqrt(w[1]) * v + np.sqrt(w[2]) * e, w

transforms = {
    1: ("normal", lambda z: z),
    2: ("normal", lambda z: z),
    3: ("normal", lambda z: z),
    4: ("lognormal (right-skewed)", lambda z: np.exp(z)),
    5: ("chi-square 3 df (right-skewed)", lambda z: stats.chi2.ppf(stats.norm.cdf(z), df=3)),
    6: ("negative lognormal (left-skewed)", lambda z: -np.exp(z)),
    7: ("t 5 df (heavy-tailed)", lambda z: stats.t.ppf(stats.norm.cdf(z), df=5)),
}
X = {}
for k, (label, f) in transforms.items():
    z, w = latent(rng)
    if k == 1:
        z = z + G_X1 * g[:, None] * tc[None, :]   # trend in x_1
    x = f(z / z.std())
    X[k] = (x - x.mean()) / x.std()
    truth.append((f"x_{k}_distribution", label))
    truth.append((f"x_{k}_shares_unit_year_noise", ", ".join(f"{s:.3f}" for s in w)))

# ---------------------------------------------------------------- x_8..x_10 (time-invariant)
target = (1.40 - 0.02) / 2.68                    # Beta median that maps to 1.40
a8 = optimize.brentq(lambda a: stats.beta.median(a, 2.0) - target, 0.5, 10)
x8 = 0.02 + 2.68 * rng.beta(a8, 2.0, N)
truth.append(("x_8_distribution", f"0.02 + 2.68*Beta({a8:.4f}, 2)"))

p9 = rng.dirichlet(np.ones(60))
x9 = rng.choice(np.arange(1, 61), size=N, p=p9)
truth.append(("x_9_distribution", "categories 1..60, Dirichlet(1) proportions"))

p10 = 1 / np.arange(1, 21) ** 1.1
p10 /= p10.sum()
x10_idx = rng.choice(np.arange(1, 21), size=N, p=p10)
x10 = np.array([f"string_{i}" for i in x10_idx])
truth.append(("x_10_distribution", "string_1..string_20, Zipf-like p ~ 1/k^1.1"))

# ---------------------------------------------------------------- coefficients
b = rng.permutation(rng.uniform(0.02, 0.39, 8))  # b_1..b_8, randomly assigned
for k in range(8):
    truth.append((f"b_{k + 1}", b[k]))
c9, c10 = rng.uniform(0.02, 0.39, 2)
eff9 = rng.normal(0, 1, 60)
eff10 = rng.normal(0, 1, 20)
truth.append(("c_9 (multiplies x_9 category effects)", c9))
truth.append(("c_10 (multiplies x_10 category effects)", c10))
for j in range(60):
    truth.append((f"x_9_category_effect_{j + 1}", eff9[j]))
for j in range(20):
    truth.append((f"x_10_category_effect_string_{j + 1}", eff10[j]))

# ---------------------------------------------------------------- untreated outcome y0
rho = 0.3
e = np.empty((N, T))
e[:, 0] = rng.normal(0, 1, N)
for t in range(1, T):
    e[:, t] = rho * e[:, t - 1] + np.sqrt(1 - rho ** 2) * rng.normal(0, 1, N)
truth.append(("error", "AR(1) within unit, rho=0.3, stationary N(0,1)"))

y0 = (unit_effect[:, None] + year_effect[None, :] + e
      + sum(b[k - 1] * X[k] for k in range(1, 8))
      + b[7] * x8[:, None]
      + c9 * eff9[x9 - 1][:, None]
      + c10 * eff10[x10_idx - 1][:, None]
      + G_Y * g[:, None] * tc[None, :]
      + H_Y * h[:, None] * (years[None, :] >= T0))

# ---------------------------------------------------------------- treatment scenarios
pre = slice(T0 - 1 - WIN, T0 - 1)                # selection window: years 4..8
m1, m2, m3 = (X[k][:, pre].mean(axis=1) for k in (1, 2, 3))
tpre = np.arange(WIN) - (WIN - 1) / 2
s1 = (X[1][:, pre] * tpre).sum(axis=1) / (tpre ** 2).sum()      # OLS slope of x_1 in the window
s1 = (s1 - s1.mean()) / s1.std()
cat9 = rng.choice(np.arange(1, 61), 10, replace=False)          # random set of 10 categories
truth.append(("x_9_selected_categories", " ".join(map(str, sorted(cat9)))))
core_O = (0.8 * m1 + 0.5 * m2 ** 2 + 0.6 * (m3 > 0.5) + 0.4 * m1 * x8
          + 0.3 * np.isin(x9, cat9) + 0.8 * s1)
score_O = core_O + rng.normal(0, 0.5, N)
score_U = core_O + 1.0 * h + rng.normal(0, 0.5, N)
truth.append(("score_O", "0.8*m_1 + 0.5*m_2^2 + 0.6*1[m_3>0.5] + 0.4*m_1*x_8 + 0.3*1[x_9 in cat9] + 0.8*std_slope(x_1, years 4-8) + N(0, sd 0.5); top 500 treated"))
truth.append(("score_U", "same core as score_O + 1.0*h_i (post-only shock) + N(0, sd 0.5); top 500 treated"))
truth.append(("score_R", "500 units chosen at random"))

treated = {
    "r": np.isin(np.arange(N), rng.choice(N, NTREAT, replace=False)),
    "o": np.isin(np.arange(N), np.argsort(-score_O)[:NTREAT]),
    "u": np.isin(np.arange(N), np.argsort(-score_U)[:NTREAT]),
}

path = np.array([0] * (T0 - 1) + [0.30, 0.73, 1.17, 1.60, 1.40, 1.20, 1.00])
for t in years:
    truth.append((f"tau_year_{t}", path[t - 1]))

# ---------------------------------------------------------------- assemble long data
df = pd.DataFrame({
    "cod_id": np.repeat(ids, T),
    "year": np.tile(years, N),
})
for k in range(1, 8):
    df[f"x_{k}"] = X[k].ravel()
df["x_8"] = np.repeat(x8, T)
df["x_9"] = np.repeat(x9, T)
df["x_10"] = np.repeat(x10, T)
df["y0"] = y0.ravel()
df["post"] = (df["year"] >= T0).astype(int)
df["event_time"] = df["year"] - T0
tau_year = np.tile(path, N)
for s, tr in treated.items():
    trl = np.repeat(tr, T).astype(int)
    df[f"treated_{s}"] = trl
    df[f"D_{s}"] = trl * df["post"]
    df[f"tau_{s}"] = trl * tau_year
    df[f"y_{s}"] = df["y0"] + df[f"tau_{s}"]
    att = df[df[f"treated_{s}"] == 1].groupby("year")[f"tau_{s}"].mean()
    for t, v in att.items():
        truth.append((f"ATT_{s}_year_{t}", v))

os.makedirs(OUT, exist_ok=True)
data_path = os.path.join(OUT, "test_data.csv")
truth_path = os.path.join(OUT, "test_data_truth.csv")
df.to_csv(data_path, index=False, float_format="%.6f")
pd.DataFrame(truth, columns=["parameter", "value"]).to_csv(truth_path, index=False)

# ---------------------------------------------------------------- summary output
pd.set_option("display.width", 160, "display.max_columns", 40)
print("shape:", df.shape)
print("columns:", list(df.columns))
print(df.describe().T[["mean", "std", "min", "50%", "max"]].round(3))
print("\nmedian x_8:", round(float(np.median(x8)), 3))
print("top-5 x_9:\n", pd.Series(x9).value_counts().head(5).to_string())
print("top-5 x_10:\n", pd.Series(x10).value_counts().head(5).to_string())
print("\nSelection-window (years 4-8) balance, treated minus untreated means:")
bal = pd.DataFrame({"m_1": m1, "m_2": m2, "m_3": m3, "x_8": x8, "slope_x1": s1,
                    "unit_effect": unit_effect, "g": g, "h": h})
for s, tr in treated.items():
    diff = bal[tr].mean() - bal[~tr].mean()
    print(f"  {s.upper()} (n treated={tr.sum()}):", "  ".join(f"{c}={v:+.3f}" for c, v in diff.items()))
for p in (data_path, truth_path):
    print(f"{os.path.basename(p)}: {os.path.getsize(p) / 1e6:.1f} MB")
