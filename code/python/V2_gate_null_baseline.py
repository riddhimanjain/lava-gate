"""V2_gate_null_baseline.py -- adversarial check of the null baseline used by the"""
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
    "P4":  (f"{BASE}/results/applied_tierB/attrition_audit.csv",            "ASD", "SCZ"),
    "P1":  (f"{BASE}/results/panel/P1/attrition_audit.csv",  "SCZ", "BIP"),
    "P2a": (f"{BASE}/results/panel/P2a/attrition_audit.csv", "MDD", "BIP"),
    "P2b": (f"{BASE}/results/panel/P2b/attrition_audit.csv",
            "MDD_noUKB", "BIP_noUKB"),
    "P3":  (f"{BASE}/results/panel/P3/attrition_audit.csv",  "ASD", "ADHD"),
    "P5":  (f"{BASE}/results/panel/P5/attrition_audit.csv",  "AN",  "ASD"),
}
ALPHA = 0.05

rows = []
for pair, (path, t1, t2) in PAIRS.items():
    if not os.path.exists(path):
        print(f"MISSING {pair}: {path}")
        continue
    a = pd.read_csv(path)
    pro = a[a["processed"] == True].copy()
    K = pro["K"].to_numpy(float)

    q = stats.chi2.sf(K, K)
    p_gate_cond_null = ALPHA / q

    for weak, other in ((t1, t2), (t2, t1)):
        cw, co = f"p_{weak}", f"p_{other}"
        if cw not in pro.columns:
            print(f"  {pair}: no column {cw}")
            continue
        est_w = pro[cw].notna().to_numpy()
        both = est_w & pro[co].notna().to_numpy()

        obs_pass_cond = float((pro.loc[both, cw] < ALPHA).mean())
        n_pass = int((pro.loc[both, cw] < ALPHA).sum())
        n_both = int(both.sum())

        old_expected = n_both * ALPHA
        new_expected = float(p_gate_cond_null[both].sum())

        rows.append(dict(
            pair=pair, trait=weak,
            role="weaker" if weak == t1 else "stronger",
            n_processed=int(len(pro)),
            n_estimable=int(est_w.sum()),
            frac_estimable=float(est_w.mean()),
            null_frac_estimable=float(q.mean()),
            n_both_estimable=n_both,
            n_pass=n_pass,
            obs_pass_rate_cond=obs_pass_cond,
            null_pass_rate_cond=float(p_gate_cond_null[both].mean()),
            old_expected_null_passes=old_expected,
            new_expected_null_passes=new_expected,
            old_enrichment=obs_pass_cond / ALPHA,
            new_enrichment=n_pass / new_expected if new_expected > 0 else np.nan,
            old_frac_passes_null=old_expected / n_pass if n_pass else np.nan,
            new_frac_passes_null=new_expected / n_pass if n_pass else np.nan,
        ))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "V2_gate_null_baseline.csv"), index=False)

pd.set_option("display.width", 200, "display.max_columns", 50)
print("\n=== CHANNEL 2 SURVIVAL: observed vs pure-null prediction ===")
print("If a trait carried NO local signal anywhere, P(estimable) = P(chi2_K >= K).")
print(df[["pair", "trait", "role", "frac_estimable",
          "null_frac_estimable"]].to_string(index=False))

print("\n=== GATE PASS RATE among BOTH-ESTIMABLE loci ===")
print("old = flat 0.05 (what the manuscript used); new = 0.05 / P(chi2_K >= K)")
print(df[["pair", "trait", "role", "n_both_estimable", "n_pass",
          "obs_pass_rate_cond", "null_pass_rate_cond",
          "old_enrichment", "new_enrichment"]].to_string(index=False))

print("\n=== FRACTION OF GATE PASSES EXPECTED NULL (criterion C3) ===")
print(df[["pair", "trait", "role", "old_frac_passes_null",
          "new_frac_passes_null"]].to_string(index=False))

print("\n=== JOINT PASSES EXPECTED IF THE WEAKER TRAIT WERE ENTIRELY NULL ===")
for pair, (path, t1, t2) in PAIRS.items():
    if not os.path.exists(path):
        continue
    a = pd.read_csv(path)
    pro = a[a["processed"] == True].copy()
    K = pro["K"].to_numpy(float)
    pcn = ALPHA / stats.chi2.sf(K, K)
    strong_pass = (pro[f"p_{t2}"] < ALPHA).to_numpy()
    both = pro[f"p_{t1}"].notna().to_numpy() & pro[f"p_{t2}"].notna().to_numpy()
    sel = strong_pass & both
    n_joint = int(((pro[f"p_{t1}"] < ALPHA) & strong_pass).sum())
    print(f"  {pair:4s} weaker={t1:10s} observed joint={n_joint:5d}   "
          f"old expected={sel.sum() * ALPHA:8.1f}   "
          f"new expected={pcn[sel].sum():8.1f}")

print("\nwritten ->", os.path.join(OUT, "V2_gate_null_baseline.csv"))
