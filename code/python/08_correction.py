"""08_correction.py -- Gate-aware correction of LAVA local genetic correlations."""
import json
import numpy as np
import pandas as pd
from scipy import stats
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()
R = machinery.SourceFileLoader("refined", "04_refined_theory.py").load_module()
S3 = machinery.SourceFileLoader("sim", "03_tierA_sim.py").load_module()

N_REF = 100_000
ALPHA = 0.05
THETA_MAX = 5.0
_GRID_CACHE = {}


def theta_obs_from_p(p, K, n_ref=N_REF, alpha=ALPHA):
    """Observed omega_hat_ii/Sigma_ii recovered from LAVA's univariate P-value."""
    nref = A.nref_scale(K, n_ref)
    stat = stats.chi2.isf(np.asarray(p, dtype=float), K)
    return stat / K - nref


def _grid(K, n_ref=N_REF, alpha=ALPHA, n=3000):
    """Cached (theta, g_K(theta), Delta(theta)) grid for one K. g_K is monotone."""
    key = (int(K), n_ref, alpha, n)
    if key in _GRID_CACHE:
        return _GRID_CACHE[key]
    th = np.concatenate([[0.0], np.geomspace(1e-7, THETA_MAX, n)])
    D = np.empty_like(th)
    for i, t in enumerate(th):
        D[i], _ = R.delta_i(K, max(t, 1e-12), n_ref=n_ref, alpha=alpha)
    g = th + D / K
    ok = np.isfinite(g) & np.isfinite(D)
    th, g, D = th[ok], g[ok], D[ok]
    order = np.argsort(g)
    out = (th[order], g[order], D[order])
    _GRID_CACHE[key] = out
    return out


def invert_theta(theta_obs, K, n_ref=N_REF, alpha=ALPHA):
    """Stage 1. Returns (theta_true, Delta, at_boundary). Vectorised over theta_obs."""
    th_g, g_g, D_g = _grid(K, n_ref, alpha)
    theta_obs = np.asarray(theta_obs, dtype=float)
    floor = g_g[0]
    th = np.interp(theta_obs, g_g, th_g, left=th_g[0], right=np.nan)
    at_bound = theta_obs <= floor
    th = np.where(at_bound, 0.0, th)
    D = np.interp(th, th_g, D_g, left=D_g[0], right=np.nan)
    return th, D, at_bound


def locus_slope_intercept(K, theta1, theta2, D1, D2, r_e=0.0):
    """SLOPE and INTERCEPT of E[rho_hat | G] = SLOPE*rho + INTERCEPT, per locus."""
    theta1 = np.asarray(theta1, float); theta2 = np.asarray(theta2, float)
    D1 = np.asarray(D1, float);         D2 = np.asarray(D2, float)
    c = np.sqrt(theta1 * theta2)
    d1 = theta1 + D1 / K
    d2 = theta2 + D2 / K
    denom = np.sqrt(d1 * d2)
    s = 1.0 + (D1 / (1.0 + 2.0 * theta1) + D2 / (1.0 + 2.0 * theta2)) / K
    slope = c * s / denom
    icept = (r_e / K) * (D1 * (1.0 + theta1) / (1.0 + 2.0 * theta1)
                         + D2 * (1.0 + theta2) / (1.0 + 2.0 * theta2)) / denom
    return slope, icept


def correct_from_theta_obs(rho_hat, K, theta1_obs, theta2_obs, r_e=0.0,
                           n_ref=N_REF, alpha=ALPHA):
    """Full two-stage correction. Returns a dict of arrays."""
    th1, D1, b1 = invert_theta(theta1_obs, K, n_ref, alpha)
    th2, D2, b2 = invert_theta(theta2_obs, K, n_ref, alpha)
    slope, icept = locus_slope_intercept(K, th1, th2, D1, D2, r_e)
    with np.errstate(invalid="ignore", divide="ignore"):
        rc = (np.asarray(rho_hat, float) - icept) / slope
    at_bound = (np.abs(rc) > 1.0) | b1 | b2
    rc_clipped = np.clip(rc, -1.0, 1.0)
    return dict(theta1=th1, theta2=th2, D1=D1, D2=D2,
                slope=slope, intercept=icept,
                rho_corrected_raw=rc, rho_corrected=rc_clipped,
                at_bound=at_bound)


def correct_locus(rho_hat, K, p1, p2, r_e=0.0, n_ref=N_REF, alpha=ALPHA):
    """Convenience wrapper taking LAVA's reported univariate P-values."""
    t1 = theta_obs_from_p(p1, K, n_ref, alpha)
    t2 = theta_obs_from_p(p2, K, n_ref, alpha)
    return correct_from_theta_obs(rho_hat, K, t1, t2, r_e, n_ref, alpha)


POOLED_SLOPE = 0.62364073
POOLED_INTERCEPT = 0.01817326

VAL_CELLS = {
    "V1_ASDlike_SCZlike": dict(K=300, theta1=0.030, theta2=0.090, r_e=0.02756),
    "V2_smallK":          dict(K=100, theta1=0.060, theta2=0.120, r_e=0.02756),
    "V3_largeK":          dict(K=500, theta1=0.025, theta2=0.070, r_e=0.02756),
    "V4_overlap_biobank": dict(K=250, theta1=0.050, theta2=0.100, r_e=0.20),
    "V5_no_overlap":      dict(K=200, theta1=0.045, theta2=0.150, r_e=0.00),
    "V6_symmetric":       dict(K=350, theta1=0.080, theta2=0.080, r_e=0.05),
}
VAL_RHOS = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]
N_VAL = 200_000


def validate():
    rng = np.random.default_rng(20260823)
    rows = []
    for name, cfg in VAL_CELLS.items():
        K, r_e = cfg["K"], cfg["r_e"]
        _grid(K)
        for rho in VAL_RHOS:
            res = S3.simulate_cell(K, cfg["theta1"], cfg["theta2"], rho, r_e,
                                   n_loci=N_VAL, rng=rng)
            passed = res["passed"]
            rh = res["rho_hat"][passed]
            o11, o22 = res["h2_1"][passed], res["h2_2"][passed]
            keep = np.isfinite(rh) & (np.abs(rh) <= S3.PARAM_LIM)
            rh = np.clip(rh[keep], -1.0, 1.0)
            o11, o22 = o11[keep], o22[keep]

            out = correct_from_theta_obs(rh, K, o11, o22, r_e)
            corr = out["rho_corrected"]
            ok = np.isfinite(corr)
            corr_ok = corr[ok]
            raw_corr = out["rho_corrected_raw"]
            raw_ok = raw_corr[np.isfinite(raw_corr)]

            th1b, D1b, _ = invert_theta(np.array([o11.mean()]), K)
            th2b, D2b, _ = invert_theta(np.array([o22.mean()]), K)
            sl_b, ic_b = locus_slope_intercept(K, th1b, th2b, D1b, D2b, r_e)
            agg = float((rh.mean() - ic_b[0]) / sl_b[0])

            pooled = np.clip((rh - POOLED_INTERCEPT) / POOLED_SLOPE, -1, 1)

            D1t, _ = R.delta_i(K, cfg["theta1"]); D2t, _ = R.delta_i(K, cfg["theta2"])
            sl_t, ic_t = locus_slope_intercept(K, cfg["theta1"], cfg["theta2"],
                                               D1t, D2t, r_e)
            oracle = np.clip((rh - ic_t) / sl_t, -1, 1)

            rows.append(dict(
                cell=name, K=K, theta1=cfg["theta1"], theta2=cfg["theta2"], r_e=r_e,
                rho_true=rho, n_reported=int(rh.size),
                mean_raw=float(rh.mean()),
                mean_corrected=float(corr_ok.mean()),
                sd_corrected=float(corr_ok.std(ddof=1)),
                se_corrected=float(corr_ok.std(ddof=1) / np.sqrt(corr_ok.size)),
                mean_pooled=float(pooled.mean()),
                mean_oracle=float(oracle.mean()),
                bias_raw=float(rh.mean() - rho),
                bias_corrected=float(corr_ok.mean() - rho),
                bias_pooled=float(pooled.mean() - rho),
                bias_oracle=float(oracle.mean() - rho),
                mean_corrected_unclipped=float(raw_ok.mean()),
                bias_corrected_unclipped=float(raw_ok.mean() - rho),
                sd_corrected_unclipped=float(raw_ok.std(ddof=1)),
                se_corrected_unclipped=float(raw_ok.std(ddof=1)/np.sqrt(raw_ok.size)),
                frac_outside_unit=float(np.mean(np.abs(raw_ok) > 1.0)),
                width_inflation=float(1.0/sl_t),
                mean_corrected_aggregate=float(agg),
                bias_corrected_aggregate=float(agg - rho),
                slope_aggregate=float(sl_b[0]),
                intercept_aggregate=float(ic_b[0]),
                theta1_agg=float(th1b[0]), theta2_agg=float(th2b[0]),
                mean_slope_perlocus=float(np.nanmean(out["slope"])),
                theory_slope=float(sl_t), theory_intercept=float(ic_t),
                mean_slope_est=float(np.nanmean(out["slope"][ok])),
                frac_theta_at_bound=float(np.mean(out["at_bound"])),
                frac_corr_pinned=float(np.mean(np.abs(corr_ok) >= 0.999)),
            ))
            print(f"{name:20s} rho={rho:+.2f}  raw={rh.mean():+.4f}  "
                  f"corr={corr_ok.mean():+.4f}  pooled={pooled.mean():+.4f}  "
                  f"agg={agg:+.4f}  slope_th={sl_t:.4f}  n={rh.size}",
                  flush=True)

    df = pd.DataFrame(rows)
    df.to_csv("out/correction_validation.csv", index=False)
    summ = {
        "n_cells": int(len(df)),
        "mean_abs_bias_raw": float(df["bias_raw"].abs().mean()),
        "mean_abs_bias_pooled": float(df["bias_pooled"].abs().mean()),
        "mean_abs_bias_corrected": float(df["bias_corrected"].abs().mean()),
        "mean_abs_bias_oracle": float(df["bias_oracle"].abs().mean()),
        "mean_abs_bias_corrected_unclipped": float(df["bias_corrected_unclipped"].abs().mean()),
        "max_abs_bias_corrected_unclipped": float(df["bias_corrected_unclipped"].abs().max()),
        "mean_abs_bias_corrected_aggregate": float(df["bias_corrected_aggregate"].abs().mean()),
        "max_abs_bias_corrected_aggregate": float(df["bias_corrected_aggregate"].abs().max()),
        "mean_frac_outside_unit": float(df["frac_outside_unit"].mean()),
        "mean_width_inflation": float(df["width_inflation"].mean()),
        "min_width_inflation": float(df["width_inflation"].min()),
        "max_width_inflation": float(df["width_inflation"].max()),
        "max_abs_bias_raw": float(df["bias_raw"].abs().max()),
        "max_abs_bias_pooled": float(df["bias_pooled"].abs().max()),
        "max_abs_bias_corrected": float(df["bias_corrected"].abs().max()),
        "max_abs_bias_oracle": float(df["bias_oracle"].abs().max()),
        "theory_slope_min": float(df["theory_slope"].min()),
        "theory_slope_max": float(df["theory_slope"].max()),
        "mean_sd_corrected": float(df["sd_corrected"].mean()),
        "pooled_slope_used": POOLED_SLOPE,
        "pooled_intercept_used": POOLED_INTERCEPT,
        "n_sim_loci_per_cell": N_VAL,
    }
    with open("out/correction_summary.json", "w") as f:
        json.dump(summ, f, indent=2)
    print("\n=== CORRECTION VALIDATION ===")
    for k, v in summ.items():
        print(f"  {k:28s} {v}")
    return df


if __name__ == "__main__":
    validate()
