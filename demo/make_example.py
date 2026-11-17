"""Generate a synthetic panel (input/panel.csv) for smoke-testing GARF (binary).

Selection rule is deliberately invisible to a level-based propensity score:
    treated if (pre-trend slope is low) AND (x1 above a threshold).
So a method that matches on the SHAPE of the pre-period (GARF) should beat one
that matches only on levels. Use this to check the install works end-to-end
before pointing the tool at your own data.
"""
import os
import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
N, T0, YEARS = 300, 2000, range(1990, 2011)

x1 = rng.normal(0, 1, N)
x2 = rng.normal(0, 1, N)
region = rng.integers(0, 4, N)
u = rng.normal(0, 1, N)                      # unit effect
slope = 0.3 + 0.25 * x2 - 0.4 * (x1 > 0.5) + rng.normal(0, 0.15, N)   # pre-trend slope

# selection on TRAJECTORY SHAPE + a threshold interaction
score = -2.0 * slope + 1.2 * (x1 > 0.3) + 0.8 * (x2 < 0) * (x1 > 0.3) + rng.normal(0, 0.4, N)
W = (score > np.quantile(score, 0.80)).astype(int)          # 20% treated

rows = []
for i in range(N):
    for t in YEARS:
        k = t - 1990
        y = 10 + u[i] + slope[i] * k + 0.5 * x1[i] + rng.normal(0, 0.3)
        if W[i] and t >= T0:
            y += 1.5 * (t - T0 + 1) * 0.4                    # true (post) effect
        rows.append({
            "id": f"U{i:03d}", "year": t,
            "W": int(W[i] and t >= T0),
            "Y": y,
            "x1": x1[i] + rng.normal(0, 0.05),               # time-varying
            "x2": x2[i] + 0.02 * k + rng.normal(0, 0.05),
            "region": f"R{region[i]}",                       # time-invariant, categorical
        })

os.makedirs("input", exist_ok=True)
pd.DataFrame(rows).to_csv("input/panel.csv", index=False)
print(f"input/panel.csv written: N={N}, treated={int(W.sum())}, T0={T0}")
