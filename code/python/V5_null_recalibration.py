"""V5_null_recalibration.py -- close the one anomaly V3 left open, and check"""
import json
import os

import numpy as np
import pandas as pd
from scipy import optimize, stats

BASE = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(BASE, "results/verification")
os.makedirs(OUT, exist_ok=True)
ALPHA = 0.05

PAIRS = {
    "P4":  (f"{BASE}/results/applied_tierB/attrition_audit.csv",            ["ASD", "SCZ"]),
    "P1":  (f"{BASE}/results/panel/P1/attrition_audit.csv",  ["SCZ", "BIP"]),
    "P2a": (f"{BASE}/results/panel/P2a/attrition_audit.csv", ["MDD", "BIP"]),
    "P2b": (f"{BASE}/results/panel/P2b/attrition_audit.csv",
            ["MDD_noUKB", "BIP_noUKB"]),
    "P3":  (f"{BASE}/results/panel/P3/attrition_audit.csv",  ["ASD", "ADHD"]),
    "P5":  (f"{BASE}/results/panel/P5/attrition_audit.csv",  ["AN", "ASD"]),
}

rows = []
for pair, (path, traits) in PAIRS.items():
    if not os.path.exists(path):
        continue
    a = pd.read_csv(path)
    pro = a[a["processed"] == True].copy()
    K = pro["K"].to_numpy(float)

    for t in traits:
        col = f"p_{t}"
        if col not in pro.columns:
            continue
        est = pro[col].notna().to_numpy()
        obs_est = float(est.mean())

        def gap(c, K=K, target=obs_est):
            return float(np.mean(stats.chi2.sf(K / c, K))) - target

        try:
            c = optimize.brentq(gap, 0.5, 3.0, xtol=1e-10)
        except ValueError:
            c = np.nan

        p = pro.loc[est, col].to_numpy(float)
        Ke = K[est]
        q_re = stats.chi2.sf(Ke / c, Ke) if np.isfinite(c) else np.nan
        thr = stats.chi2.ppf(1 - ALPHA, Ke)
        g_re = stats.chi2.sf(thr / c, Ke) if np.isfinite(c) else np.nan

        u_obs = stats.chi2.isf(p, Ke)
        u_re = np.clip(stats.chi2.sf(u_obs / c, Ke) / q_re, 0, 1) \
            if np.isfinite(c) else np.full_like(p, np.nan)
        ks_re = stats.kstest(u_re, "uniform") if np.isfinite(c) else None

        q1 = stats.chi2.sf(Ke, Ke)
        u1 = np.clip(p / q1, 0, 1)
        ks1 = stats.kstest(u1, "uniform")

        obs_gate = float((p < ALPHA).mean())
        rows.append(dict(
            pair=pair, trait=t, n_estimable=int(est.sum()),
            obs_frac_estimable=obs_est,
            null_frac_estimable_c1=float(stats.chi2.sf(K, K).mean()),
            c_hat=float(c),
            ks_p_c1=float(ks1.pvalue), mean_u_c1=float(u1.mean()),
            ks_p_recal=float(ks_re.pvalue) if ks_re else np.nan,
            mean_u_recal=float(np.mean(u_re)) if np.isfinite(c) else np.nan,
            obs_gate_rate_cond=obs_gate,
            null_gate_rate_c1=float(np.mean(ALPHA / q1)),
            null_gate_rate_recal=float(np.mean(g_re / q_re))
            if np.isfinite(c) else np.nan,
        ))

df = pd.DataFrame(rows)
df["enrich_c1"] = df.obs_gate_rate_cond / df.null_gate_rate_c1
df["enrich_recal"] = df.obs_gate_rate_cond / df.null_gate_rate_recal
df.to_csv(os.path.join(OUT, "V5_null_recalibration.csv"), index=False)

pd.set_option("display.width", 250, "display.max_columns", 60)
print("=== VARIANCE CALIBRATION FACTOR c  (U ~ c * chi2_K under a complete null) ===")
print("c < 1 : the chi-square reference overstates the sampling variance.")
print("c > 1 : indistinguishable from the trait simply having real signal.\n")
print(df[["pair", "trait", "obs_frac_estimable", "null_frac_estimable_c1",
          "c_hat"]].to_string(index=False))

print("\n=== COMPLETE-NULL GOODNESS OF FIT, before and after recalibration ===")
print(df[["pair", "trait", "mean_u_c1", "ks_p_c1", "mean_u_recal",
          "ks_p_recal"]].to_string(index=False))

print("\n=== GATE ENRICHMENT over the conditional null (criterion C3) ===")
print("manuscript used a flat 0.05 baseline, giving 'enrich_manuscript'.\n")
print(df.assign(enrich_manuscript=lambda d: d.obs_gate_rate_cond / ALPHA)[
    ["pair", "trait", "obs_gate_rate_cond", "enrich_manuscript",
     "enrich_c1", "enrich_recal"]].to_string(index=False))

json.dump(df.to_dict("records"),
          open(os.path.join(OUT, "V5_null_recalibration.json"), "w"), indent=2)
print("\nwritten ->", OUT)
