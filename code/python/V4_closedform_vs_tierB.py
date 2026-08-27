"""V4_closedform_vs_tierB.py -- test the closed form of Section 2.3 against Tier B"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gammaln

BASE = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
VER = os.path.join(BASE, "results/verification")
os.makedirs(VER, exist_ok=True)
ALPHA = 0.05

SRC = {
    "P4":  f"{BASE}/results/tierB_extended/tierB_extended.csv",
    "P1":  f"{BASE}/results/tierB_extended_multipair/P1/tierB_extended.csv",
    "P2a": f"{BASE}/results/tierB_extended_multipair/P2a/tierB_extended.csv",
    "P2b": f"{BASE}/results/tierB_extended_multipair/P2b/tierB_extended.csv",
    "P3":  f"{BASE}/results/tierB_extended_multipair/P3/tierB_extended.csv",
    "P5":  f"{BASE}/results/tierB_extended_multipair/P5/tierB_extended.csv",
}
SUMM = {k: v.replace("tierB_extended.csv", "tierB_extended_summary.json")
        for k, v in SRC.items()}


def _pois_w(lam, tol=1e-14):
    half = lam / 2.0
    jmax = int(max(60, half + 12.0 * np.sqrt(max(half, 1.0))))
    j = np.arange(jmax + 1)
    logw = -half + j * np.log(half + 1e-300) - gammaln(j + 1.0)
    w = np.exp(logw)
    if w.max() > 0:
        keep = w > tol * w.max()
        j, w = j[keep], w[keep]
    return j, w / w.sum()


def delta_i(K, theta, alpha=ALPHA):
    """Delta = E[S_ii | S_ii > t] - E[S_ii].  nref.scale == 1 in LAVA v0.1.5."""
    t = stats.chi2.ppf(1 - alpha, K)
    j, w = _pois_w(K * theta)
    m = K + 2.0 * j
    p_tail = float(np.sum(w * stats.chi2.sf(t, m)))
    e_trunc = float(np.sum(w * m * stats.chi2.sf(t, m + 2.0)))
    if p_tail <= 0:
        return np.nan, 0.0
    return e_trunc / p_tail - K * (1.0 + theta), p_tail


def slope_intercept(K, th1, th2, r_e):
    """S_locus and I_locus of Section 2.3."""
    D1, p1 = delta_i(K, th1)
    D2, p2 = delta_i(K, th2)
    if not np.isfinite(D1) or not np.isfinite(D2):
        return np.nan, np.nan, np.nan
    d1, d2 = th1 + D1 / K, th2 + D2 / K
    if d1 <= 0 or d2 <= 0:
        return np.nan, np.nan, np.nan
    den = np.sqrt(d1 * d2)
    S = np.sqrt(th1 * th2) * (1 + (D1 / (1 + 2 * th1) + D2 / (1 + 2 * th2)) / K) / den
    I = (r_e / K) * (D1 * (1 + th1) / (1 + 2 * th1)
                     + D2 * (1 + th2) / (1 + 2 * th2)) / den
    return S, I, p1 * p2


rows, per_locus = [], []
for pair, path in SRC.items():
    if not os.path.exists(path):
        print("MISSING", pair)
        continue
    tb = pd.read_csv(path)
    loci = tb.drop_duplicates("locus")[["locus", "K", "theta1", "theta2", "r_e"]]
    usable = set(tb.loc[tb["mean_rho_hat"].notna(), "locus"])
    loci = loci[loci["locus"].isin(usable)]

    S, I, PG = [], [], []
    for _, r in loci.iterrows():
        s, i, pg = slope_intercept(int(r.K), float(r.theta1), float(r.theta2),
                                   float(r.r_e))
        S.append(s); I.append(i); PG.append(pg)
        per_locus.append(dict(pair=pair, locus=int(r.locus), K=int(r.K),
                              theta1=r.theta1, theta2=r.theta2, r_e=r.r_e,
                              cf_slope=s, cf_intercept=i, p_gate_joint=pg))
    S, I = np.array(S, float), np.array(I, float)

    sm = json.load(open(SUMM[pair]))
    v1 = None
    p1f = os.path.join(VER, "V1_summary.csv")
    if os.path.exists(p1f):
        v = pd.read_csv(p1f)
        v = v[v.pair == pair]
        if len(v):
            v1 = v.iloc[0]

    rows.append(dict(
        pair=pair, n_loci=int(len(loci)),
        r_e=float(loci["r_e"].iloc[0]),
        cf_slope_mean=float(np.nanmean(S)),
        cf_slope_min=float(np.nanmin(S)), cf_slope_max=float(np.nanmax(S)),
        cf_intercept_mean=float(np.nanmean(I)),
        rule_0p68_intercept=0.68 * float(loci["r_e"].iloc[0]),
        tierB_slope_D=float(sm["slope_full"]),
        tierB_intercept_D=float(sm["intercept_full"]),
        v1_slope_C_notrunc=float(v1["slope_C_gated_raw"]) if v1 is not None else np.nan,
        v1_int_C_notrunc=float(v1["int_C"]) if v1 is not None else np.nan,
        v1_slope_D=float(v1["slope_D_gated_trunc"]) if v1 is not None else np.nan,
        v1_int_D=float(v1["int_D"]) if v1 is not None else np.nan,
    ))

df = pd.DataFrame(rows)
df["cf_minus_C_slope"] = df.cf_slope_mean - df.v1_slope_C_notrunc
df["cf_minus_C_int"] = df.cf_intercept_mean - df.v1_int_C_notrunc
df["rule_error_int"] = df.rule_0p68_intercept - df.tierB_intercept_D
df.to_csv(os.path.join(VER, "V4_closedform_vs_tierB.csv"), index=False)
pd.DataFrame(per_locus).to_csv(os.path.join(VER, "V4_per_locus_closedform.csv"),
                               index=False)

pd.set_option("display.width", 250, "display.max_columns", 60)
print("=== SLOPE: closed form vs Tier B (same loci) ===")
print("cf is compared to arm C (gate, NO truncation) -- what the theory models.")
print(df[["pair", "n_loci", "cf_slope_mean", "v1_slope_C_notrunc",
          "cf_minus_C_slope", "v1_slope_D", "tierB_slope_D"]].to_string(index=False))

print("\n=== INTERCEPT: the 0.68 x r_e rule vs the real closed form ===")
print(df[["pair", "r_e", "rule_0p68_intercept", "cf_intercept_mean",
          "v1_int_C_notrunc", "tierB_intercept_D",
          "rule_error_int"]].to_string(index=False))

print("\n=== closed-form slope heterogeneity across real loci, per pair ===")
print(df[["pair", "cf_slope_min", "cf_slope_max"]].to_string(index=False))
print("\nwritten ->", VER)
