"""V3_null_calibration.py -- is the chi-square null that criterion C3 and the"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(BASE, "results/verification")
os.makedirs(OUT, exist_ok=True)

PAIRS = {
    "P4":  (f"{BASE}/results/applied_tierB/attrition_audit.csv",            ["ASD", "SCZ"]),
    "P1":  (f"{BASE}/results/panel/P1/attrition_audit.csv",  ["SCZ", "BIP"]),
    "P2a": (f"{BASE}/results/panel/P2a/attrition_audit.csv", ["MDD", "BIP"]),
    "P2b": (f"{BASE}/results/panel/P2b/attrition_audit.csv",
            ["MDD_noUKB", "BIP_noUKB"]),
    "P3":  (f"{BASE}/results/panel/P3/attrition_audit.csv",  ["ASD", "ADHD"]),
    "P5":  (f"{BASE}/results/panel/P5/attrition_audit.csv",  ["AN", "ASD"]),
}
ALPHA = 0.05

rows = []
for pair, (path, traits) in PAIRS.items():
    if not os.path.exists(path):
        print(f"MISSING {pair}")
        continue
    a = pd.read_csv(path)
    pro = a[a["processed"] == True].copy()
    K = pro["K"].to_numpy(float)
    q = stats.chi2.sf(K, K)

    for t in traits:
        col = f"p_{t}"
        if col not in pro.columns:
            continue
        m = pro[col].notna().to_numpy()
        p = pro.loc[m, col].to_numpy(float)
        qq = q[m]
        u = np.clip(p / qq, 0, 1)

        ks = stats.kstest(u, "uniform")
        rows.append(dict(
            pair=pair, trait=t, n_estimable=int(m.sum()),
            frac_estimable=float(m.mean()), null_frac_estimable=float(q.mean()),
            estimability_deficit=float(m.mean() - q.mean()),
            mean_u=float(u.mean()), median_u=float(np.median(u)),
            frac_u_below_0p1=float((u < 0.1).mean()),
            ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
            n_gate_pass=int((p < ALPHA).sum()),
            obs_gate_rate_cond=float((p < ALPHA).mean()),
            null_gate_rate_cond=float((ALPHA / qq).mean()),
        ))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "V3_null_calibration.csv"), index=False)

pd.set_option("display.width", 220, "display.max_columns", 60)
print("=== GOODNESS OF FIT OF THE COMPLETE LOCAL NULL ===")
print("u = p / P(chi2_K >= K) is Uniform(0,1) if the trait has no local signal.")
print("mean_u = 0.5 under the null; < 0.5 means real signal; > 0.5 means the")
print("chi-square model overstates the sampling variance.\n")
print(df[["pair", "trait", "n_estimable", "mean_u", "median_u",
          "frac_u_below_0p1", "ks_stat", "ks_p"]].to_string(index=False))

print("\n=== ESTIMABILITY: observed minus complete-null prediction ===")
print("A NEGATIVE deficit is impossible under the model: real signal can only")
print("raise P(omega.hat >= 0). Negative values indicate model mis-calibration.\n")
print(df[["pair", "trait", "frac_estimable", "null_frac_estimable",
          "estimability_deficit"]].to_string(index=False))

print("\n=== GATE PASS RATE, observed vs correct conditional null ===")
print(df[["pair", "trait", "n_gate_pass", "obs_gate_rate_cond",
          "null_gate_rate_cond"]].assign(
              enrichment=lambda d: d.obs_gate_rate_cond / d.null_gate_rate_cond
          ).to_string(index=False))

json.dump({r["pair"] + "_" + r["trait"]: r for r in rows},
          open(os.path.join(OUT, "V3_null_calibration.json"), "w"), indent=2)
print("\nwritten ->", OUT)
