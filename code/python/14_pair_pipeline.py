"""14_pair_pipeline.py -- the applied protocol, for any trait pair."""
import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
from importlib import machinery

_HERE = os.path.dirname(os.path.abspath(__file__))
_prev = os.getcwd()
os.chdir(_HERE)
A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()
R = machinery.SourceFileLoader("refined", "04_refined_theory.py").load_module()
C = machinery.SourceFileLoader("corr", "08_correction.py").load_module()
CV = machinery.SourceFileLoader("cov", "06_coverage.py").load_module()
os.chdir(_prev)

ALPHA = 0.05
SEED = 20260823
N_NULL = 40_000
N_SIGN = 6_000
N_DET = 250
N_BOOT_CI = 2_000

C1_MAX_SLOPE_FOLD = 1.25
C2_MIN_POOLED_Z = 5.0
C2_MAX_FRAC_BELOW = 0.40
C3_MAX_FRAC_NULL = 0.40


def load_reported(audit_path, t1, t2):
    """The loci this pair actually reports, with their two univariate P-values."""
    aud = pd.read_csv(audit_path)
    p1, p2 = f"p_{t1}", f"p_{t2}"
    for c in ("gate_pass", "rho_hat", "K", p1, p2):
        if c not in aud.columns:
            raise KeyError(f"{audit_path} has no column '{c}'. "
                           f"Columns present: {list(aud.columns)}")
    gp = aud["gate_pass"]
    if gp.dtype == object:
        gp = gp.astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    aud["gate_pass"] = gp.astype("boolean").fillna(False).astype(bool)
    aud["processed"] = (aud["processed"].astype(str).str.upper() == "TRUE"
                        if aud["processed"].dtype == object
                        else aud["processed"].astype(bool))
    rep = aud[aud["gate_pass"] & np.isfinite(aud["rho_hat"])].copy()
    rep["K"] = rep["K"].astype(int)
    rep["t1_obs"] = C.theta_obs_from_p(rep[p1].values, rep["K"].values)
    rep["t2_obs"] = C.theta_obs_from_p(rep[p2].values, rep["K"].values)
    return aud, rep


def apply_correction(rep, r_e, t1, t2, outdir):
    rows = []
    for r in rep.itertuples():
        K = int(r.K)
        out = C.correct_from_theta_obs(np.array([r.rho_hat]), K,
                                       np.array([r.t1_obs]), np.array([r.t2_obs]), r_e)
        sl = float(out["slope"][0])
        inv = (1.0 / sl) if sl > 0 else np.nan
        w = (r.rho_upper - r.rho_lower) if np.isfinite(getattr(r, "rho_upper", np.nan)) else np.nan
        rows.append(dict(
            locus=r.locus, chr=r.chr, K=K,
            p_1=getattr(r, f"p_{t1}"), p_2=getattr(r, f"p_{t2}"),
            theta_1=float(out["theta1"][0]), theta_2=float(out["theta2"][0]),
            slope=sl, intercept=float(out["intercept"][0]), floored=bool(sl <= 0),
            rho_hat=r.rho_hat,
            rho_corrected=float(out["rho_corrected"][0]),
            rho_corrected_unclipped=float(out["rho_corrected_raw"][0]),
            rho_lower=getattr(r, "rho_lower", np.nan),
            rho_upper=getattr(r, "rho_upper", np.nan),
            ci_width=w, ci_width_corrected=w * inv))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(outdir, "applied_correction.csv"), index=False)

    Kbar = int(round(rep["K"].mean()))
    th1b, D1b, _ = C.invert_theta(np.array([rep["t1_obs"].mean()]), Kbar)
    th2b, D2b, _ = C.invert_theta(np.array([rep["t2_obs"].mean()]), Kbar)
    sl_b, ic_b = C.locus_slope_intercept(Kbar, th1b, th2b, D1b, D2b, r_e)
    agg = float((rep["rho_hat"].mean() - ic_b[0]) / sl_b[0]) if sl_b[0] > 0 else np.nan

    live = df.loc[~df["floored"], "slope"]
    summ = dict(
        n_loci=int(len(df)), r_e=r_e, K_mean=Kbar,
        mean_rho_hat=float(rep["rho_hat"].mean()),
        n_floored_no_recoverable_signal=int(df["floored"].sum()),
        per_locus_slope_min=float(live.min()) if len(live) else None,
        per_locus_slope_max=float(live.max()) if len(live) else None,
        per_locus_slope_median=float(live.median()) if len(live) else None,
        per_locus_intercept_min=float(df["intercept"].min()),
        per_locus_intercept_max=float(df["intercept"].max()),
        per_locus_intercept_median=float(df["intercept"].median()),
        aggregate_slope=float(sl_b[0]), aggregate_intercept=float(ic_b[0]),
        mean_rho_corrected_aggregate=agg,
        n_corrected_outside_unit=int((df["rho_corrected_unclipped"].abs() > 1).sum()),
        mean_ci_width=float(np.nanmean(df["ci_width"])),
        mean_ci_width_corrected=float(np.nanmean(df["ci_width_corrected"])),
        frac_corrected_width_over_2=float((df["ci_width_corrected"].dropna() > 2).mean())
        if df["ci_width_corrected"].notna().any() else None,
    )
    with open(os.path.join(outdir, "applied_correction_summary.json"), "w") as f:
        json.dump(summ, f, indent=2)
    return df, summ


def strat_slope(d, n_strata, r_e):
    bins = (np.zeros(len(d), int) if n_strata == 1 else
            np.clip(np.digitize(d["K"], np.quantile(d["K"],
                    np.linspace(0, 1, n_strata + 1))[1:-1]), 0, n_strata - 1))
    sl_w = ic_w = 0.0
    for b in np.unique(bins):
        g = d[bins == b]
        Kb = int(round(g["K"].mean() / 25.0) * 25) or 25
        th1, D1, _ = C.invert_theta(np.array([g["t1_obs"].mean()]), Kb)
        th2, D2, _ = C.invert_theta(np.array([g["t2_obs"].mean()]), Kb)
        sl, ic = C.locus_slope_intercept(Kb, th1, th2, D1, D2, r_e)
        sl_w += float(sl[0]) * len(g); ic_w += float(ic[0]) * len(g)
    return sl_w / len(d), ic_w / len(d)


def diagnostic(rep, r_e, t1, t2, outdir):
    rng = np.random.default_rng(SEED)

    stab = {}
    for ns in (1, 2, 3, 4, 6, 8):
        sl, ic = strat_slope(rep, ns, r_e)
        cm = (rep["rho_hat"].mean() - ic) / sl if sl > 0 else np.nan
        stab[f"strata_{ns}"] = dict(slope=float(sl), intercept=float(ic),
                                    corrected_mean=float(cm))
    sl_vals = np.array([v["slope"] for v in stab.values()])
    cm_vals = np.array([v["corrected_mean"] for v in stab.values()])
    fold = float(sl_vals.max() / sl_vals.min()) if sl_vals.min() > 0 else np.inf

    kbin = (rep["K"] / 50).round() * 50
    rows = []
    for Kb in sorted({int(k) or 50 for k in kbin}):
        nref = A.nref_scale(Kb, C.N_REF)
        t = A.gate_threshold(Kb, ALPHA) / nref
        u = rng.chisquare(Kb, N_NULL)
        u = u[u > t]
        th_null = u / Kb - nref
        for trait, col in ((t1, "t1_obs"), (t2, "t2_obs")):
            sub = rep.loc[kbin == Kb, col]
            if not len(sub):
                continue
            rows.append(dict(
                K_bin=Kb, trait=trait, n_loci=int(len(sub)),
                observed_mean_theta=float(sub.mean()),
                null_mean_theta=float(th_null.mean()),
                null_sd_theta=float(th_null.std(ddof=1)),
                z_above_null=float((sub.mean() - th_null.mean())
                                   / (th_null.std(ddof=1) / np.sqrt(len(sub)))),
                pct_of_loci_below_null_mean=float((sub < th_null.mean()).mean())))
    fl = pd.DataFrame(rows)
    fl.to_csv(os.path.join(outdir, "floor_diagnostic.csv"), index=False)

    def pooled(trait):
        s = fl[fl.trait == trait]
        if not len(s):
            return dict(pooled_z=np.nan, frac_below=np.nan,
                        z_min=np.nan, z_max=np.nan)
        w = s["n_loci"].values
        return dict(pooled_z=float((s["z_above_null"] * w).sum() / w.sum()),
                    frac_below=float((s["pct_of_loci_below_null_mean"] * w).sum() / w.sum()),
                    z_min=float(s["z_above_null"].min()),
                    z_max=float(s["z_above_null"].max()))

    f1, f2 = pooled(t1), pooled(t2)
    summ = dict(
        n_loci=int(len(rep)), r_e=r_e,
        mean_rho_hat_observed=float(rep["rho_hat"].mean()),
        stability=stab,
        slope_range_over_binning=[float(sl_vals.min()), float(sl_vals.max())],
        slope_fold_change_over_binning=fold,
        corrected_mean_range_over_binning=[float(np.nanmin(cm_vals)),
                                           float(np.nanmax(cm_vals))],
        floor={t1: f1, t2: f2},
    )
    with open(os.path.join(outdir, "aggregate_correction.json"), "w") as f:
        json.dump(summ, f, indent=2)
    return summ


def signal_bound(aud, t1, t2, outdir):
    pro = aud[aud["processed"]]
    both = pro[np.isfinite(pro[f"p_{t1}"]) & np.isfinite(pro[f"p_{t2}"])]
    out = {"alpha": ALPHA, "n_processed": int(len(pro)),
           "n_both_estimable": int(len(both))}
    for trait in (t1, t2):
        p = both[f"p_{trait}"]
        n_pass = int((p < ALPHA).sum())
        obs = float((p < ALPHA).mean()) if len(both) else np.nan
        exp_null = ALPHA * len(both)
        out[trait] = dict(
            pass_rate=obs, n_pass=n_pass,
            expected_null_passes=float(exp_null),
            fraction_of_passes_expected_null=float(min(1.0, exp_null / max(n_pass, 1))),
            enrichment_over_alpha=float(obs / ALPHA) if np.isfinite(obs) else np.nan,
            pi_upper_bound_perfect_power=float(max(0.0, (obs - ALPHA) / (1 - ALPHA))))
    n_joint = int(((both[f"p_{t1}"] < ALPHA) & (both[f"p_{t2}"] < ALPHA)).sum())
    out["joint"] = dict(
        n_joint_pass=n_joint,
        expected_joint_if_t1_all_null=float(ALPHA * (both[f"p_{t2}"] < ALPHA).sum()),
        expected_joint_if_t2_all_null=float(ALPHA * (both[f"p_{t1}"] < ALPHA).sum()))
    with open(os.path.join(outdir, "signal_bound.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


def sim_gated(K, theta1, theta2, rho, r_e, n_want, rng, max_tries=60):
    """Draw until n_want gated loci are obtained. Returns (rho_hat, omega, n_drawn)."""
    nref = A.nref_scale(K, C.N_REF)
    t_crit = A.gate_threshold(K, ALPHA)
    Sigma = np.array([[1.0, r_e], [r_e, 1.0]])
    L_sig = np.linalg.cholesky(Sigma)
    Om = np.array([[theta1, rho * np.sqrt(theta1 * theta2)],
                   [rho * np.sqrt(theta1 * theta2), theta2]])
    M1 = np.linalg.cholesky(Om + 1e-15 * np.eye(2)).T * np.sqrt(K)
    m = K - 2
    rh_all, om_all, tries, drawn, have = [], [], 0, 0, 0
    batch = max(20000, n_want * 20)
    while have < n_want and tries < max_tries:
        n = batch
        Z1 = rng.standard_normal((n, 2, 2))
        a11 = np.sqrt(rng.chisquare(m, n)); a22 = np.sqrt(rng.chisquare(m - 1, n))
        a21 = rng.standard_normal(n)
        V = np.empty((n, 2, 2))
        V[:, 0, 0] = a11 ** 2
        V[:, 0, 1] = V[:, 1, 0] = a11 * a21
        V[:, 1, 1] = a21 ** 2 + a22 ** 2
        ZtZ = np.einsum("nki,nkj->nij", Z1, Z1) + V
        MtZ = np.einsum("ik,nkj->nij", M1.T, Z1)
        cross = MtZ @ L_sig.T
        S = (M1.T @ M1)[None] + cross + np.swapaxes(cross, 1, 2) + L_sig @ ZtZ @ L_sig.T
        passed = ((S[:, 0, 0] / Sigma[0, 0] * nref > t_crit) &
                  (S[:, 1, 1] / Sigma[1, 1] * nref > t_crit))
        om = (S / K - Sigma[None]) * nref
        om = om[passed]
        with np.errstate(invalid="ignore", divide="ignore"):
            rh = om[:, 0, 1] / np.sqrt(om[:, 0, 0] * om[:, 1, 1])
        keep = np.isfinite(rh) & (np.abs(rh) <= C.PARAM_LIM if hasattr(C, "PARAM_LIM")
                                  else np.abs(rh) <= 1.25)
        rh = np.clip(rh[keep], -1.0, 1.0); om = om[keep]
        rh_all.append(rh); om_all.append(om)
        have += rh.size; drawn += n; tries += 1
    rh = np.concatenate(rh_all)[:n_want] if rh_all else np.array([])
    om = np.concatenate(om_all)[:n_want] if om_all else np.zeros((0, 2, 2))
    return rh, om, drawn


def poisson_binomial_pmf(p):
    pmf = np.zeros(len(p) + 1); pmf[0] = 1.0
    for pi in p:
        pmf[1:] = pmf[1:] * (1 - pi) + pmf[:-1] * pi
        pmf[0] *= (1 - pi)
    return pmf


def pb_two_sided_p(k_obs, p):
    pmf = poisson_binomial_pmf(p)
    tol = 1e-9 * pmf.max()
    return float(pmf[pmf <= pmf[k_obs] + tol].sum())


def null_sign(rep, r_e, t1, t2, outdir, do_part_b=True):
    rng = np.random.default_rng(SEED)
    th1, th2 = [], []
    for r in rep.itertuples():
        K = int(r.K)
        a, _, _ = C.invert_theta(np.array([r.t1_obs]), K)
        b, _, _ = C.invert_theta(np.array([r.t2_obs]), K)
        th1.append(max(float(a[0]), 1e-6)); th2.append(max(float(b[0]), 1e-6))
    rep = rep.copy(); rep["theta_1"] = th1; rep["theta_2"] = th2

    n_obs = len(rep)
    n_pos_obs = int((rep["rho_hat"] > 0).sum())
    det_obs = (np.isfinite(rep["rho_lower"]) & np.isfinite(rep["rho_upper"]) &
               ((rep["rho_lower"] > 0) | (rep["rho_upper"] < 0)))
    n_det_obs = int(det_obs.sum())
    n_det_pos_obs = int((det_obs & (rep["rho_hat"] > 0)).sum())

    p_pos = []
    for r in rep.itertuples():
        rh, _, _ = sim_gated(int(r.K), r.theta_1, r.theta_2, 0.0, r_e, N_SIGN, rng)
        p_pos.append(float((rh > 0).mean()) if rh.size else np.nan)
    p_pos = np.array(p_pos)
    partA = dict(mean_p_positive=float(np.nanmean(p_pos)),
                 min_p_positive=float(np.nanmin(p_pos)),
                 max_p_positive=float(np.nanmax(p_pos)),
                 expected_n_positive=float(np.nansum(p_pos)),
                 observed_n_positive=n_pos_obs, n_loci=n_obs,
                 poisson_binomial_two_sided_p=pb_two_sided_p(n_pos_obs,
                                                             np.nan_to_num(p_pos, nan=0.5)),
                 naive_binomial_null=0.5)

    partB = None
    if do_part_b:
        Sigma = np.array([[1.0, r_e], [r_e, 1.0]])
        p_det, p_det_pos = [], []
        for r in rep.itertuples():
            rh, om, _ = sim_gated(int(r.K), r.theta_1, r.theta_2, 0.0, r_e, N_DET, rng)
            nd = npd = 0
            for j in range(rh.size):
                lo, hi = CV.lava_ci(om[j], Sigma, int(r.K), rng, n_boot=N_BOOT_CI)
                if not np.isfinite(lo):
                    continue
                if lo > 0 or hi < 0:
                    nd += 1
                    if rh[j] > 0:
                        npd += 1
            p_det.append(nd / max(rh.size, 1)); p_det_pos.append(npd / max(rh.size, 1))
        p_det = np.array(p_det); p_det_pos = np.array(p_det_pos)
        pooled_pos_given_det = (float(np.nansum(p_det_pos) / np.nansum(p_det))
                                if np.nansum(p_det) > 0 else np.nan)
        p_all_pos = pooled_pos_given_det ** n_det_obs
        p_all_neg = (1.0 - pooled_pos_given_det) ** n_det_obs
        partB = dict(
            null_mean_p_determinate=float(p_det.mean()),
            null_expected_n_determinate=float(p_det.sum()),
            observed_n_determinate=n_det_obs,
            poisson_binomial_two_sided_p_determinate=pb_two_sided_p(n_det_obs, p_det),
            null_p_positive_given_determinate=pooled_pos_given_det,
            observed_n_determinate_positive=n_det_pos_obs,
            direction_test_p=float(min(1.0, p_all_pos +
                                       (p_all_neg if p_all_neg <= p_all_pos else 0.0))))

    out = dict(observed=dict(n_loci=n_obs, n_positive=n_pos_obs,
                             n_determinate=n_det_obs,
                             n_determinate_positive=n_det_pos_obs),
               r_e=r_e, traits=[t1, t2], partA=partA, partB=partB)
    with open(os.path.join(outdir, "null_sign_summary.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out


def verdict(diag, sig, t1, t2):
    """Score the three frozen criteria. All three must pass."""
    c1_val = diag["slope_fold_change_over_binning"]
    c1 = bool(c1_val < C1_MAX_SLOPE_FOLD)

    f = diag["floor"]
    c2_parts = {}
    for t in (t1, t2):
        c2_parts[t] = bool(f[t]["pooled_z"] >= C2_MIN_POOLED_Z and
                           f[t]["frac_below"] <= C2_MAX_FRAC_BELOW)
    c2 = all(c2_parts.values())

    c3_parts = {t: bool(sig[t]["fraction_of_passes_expected_null"] <= C3_MAX_FRAC_NULL)
                for t in (t1, t2)}
    c3 = all(c3_parts.values())

    usable = c1 and c2 and c3
    fails = [n for n, ok in (("C1 stability", c1), ("C2 floor distance", c2),
                             ("C3 gate signal", c3)) if not ok]
    return dict(
        usable=usable,
        C1_stability=dict(passed=c1, slope_fold_change=c1_val,
                          threshold=C1_MAX_SLOPE_FOLD),
        C2_floor=dict(passed=c2, per_trait=c2_parts,
                      pooled_z={t: f[t]["pooled_z"] for t in (t1, t2)},
                      frac_below={t: f[t]["frac_below"] for t in (t1, t2)},
                      threshold_z=C2_MIN_POOLED_Z,
                      threshold_frac_below=C2_MAX_FRAC_BELOW),
        C3_signal=dict(passed=c3, per_trait=c3_parts,
                       frac_null={t: sig[t]["fraction_of_passes_expected_null"]
                                  for t in (t1, t2)},
                       threshold=C3_MAX_FRAC_NULL),
        failed_criteria=fails,
        statement=("A corrected magnitude MAY be reported for this pair."
                   if usable else
                   "NOT APPLICABLE. Report the uncorrected estimate with the "
                   "attenuation stated as a bound. Failed: " + ", ".join(fails)))


def run_pair(pair, t1, t2, r_e, audit, outdir, do_null_sign=False,
             do_part_b=True, quiet=False):
    os.makedirs(outdir, exist_ok=True)
    aud, rep = load_reported(audit, t1, t2)
    if not quiet:
        print(f"\n=== {pair}: {t1} x {t2}   r_e = {r_e:.5f} ===")
        print(f"  blocks {len(aud)}, processed {int(aud['processed'].sum())}, "
              f"reported {len(rep)}")
    if len(rep) == 0:
        summ = dict(pair=pair, traits=[t1, t2], r_e=r_e, n_reported=0,
                    note="No locus reported a rho for this pair; nothing to correct.")
        with open(os.path.join(outdir, "pair_summary.json"), "w") as f:
            json.dump(summ, f, indent=2)
        return summ

    _, corr = apply_correction(rep, r_e, t1, t2, outdir)
    diag = diagnostic(rep, r_e, t1, t2, outdir)
    sig = signal_bound(aud, t1, t2, outdir)
    ver = verdict(diag, sig, t1, t2)
    ns = null_sign(rep, r_e, t1, t2, outdir, do_part_b) if do_null_sign else None

    pro = int(aud["processed"].sum())
    summ = dict(
        pair=pair, traits=[t1, t2], r_e=r_e,
        predicted_gate_intercept=0.68 * r_e,
        fitted_aggregate_intercept=corr["aggregate_intercept"],
        n_blocks=int(len(aud)), n_processed=pro,
        n_reported=int(len(rep)),
        pct_reported=round(100 * len(rep) / len(aud), 2),
        mean_rho_hat=corr["mean_rho_hat"],
        aggregate_slope=corr["aggregate_slope"],
        mean_rho_corrected_aggregate=corr["mean_rho_corrected_aggregate"],
        n_floored=corr["n_floored_no_recoverable_signal"],
        mean_ci_width=corr["mean_ci_width"],
        mean_ci_width_corrected=corr["mean_ci_width_corrected"],
        diagnostic=ver, signal=sig, null_sign=ns)
    with open(os.path.join(outdir, "pair_summary.json"), "w") as f:
        json.dump(summ, f, indent=2)

    if not quiet:
        print(f"  mean rho_hat            {corr['mean_rho_hat']:+.4f}")
        print(f"  aggregate slope         {corr['aggregate_slope']:.4f}")
        print(f"  predicted intercept     {0.68*r_e:.4f}  "
              f"fitted {corr['aggregate_intercept']:.4f}")
        print(f"  C1 slope fold-change    {ver['C1_stability']['slope_fold_change']:.3f}"
              f"  (< {C1_MAX_SLOPE_FOLD})  {'PASS' if ver['C1_stability']['passed'] else 'FAIL'}")
        for t in (t1, t2):
            print(f"  C2 {t:<10s} pooled z {ver['C2_floor']['pooled_z'][t]:+7.2f}  "
                  f"below-null {ver['C2_floor']['frac_below'][t]:.2f}")
        for t in (t1, t2):
            print(f"  C3 {t:<10s} frac of passes expected null "
                  f"{ver['C3_signal']['frac_null'][t]:.3f}")
        print(f"  VERDICT: {'USABLE' if ver['usable'] else 'NOT APPLICABLE'} "
              f"{ver['failed_criteria']}")
    return summ


def selftest():
    """Reproduce the published P4 numbers through this generic code path."""
    base = os.path.abspath(os.path.join(_HERE, "..", ".."))
    audit = os.path.join(base, "results", "applied_tierB", "attrition_audit.csv")
    ref = json.load(open(os.path.join(base, "results", "correction",
                                      "aggregate_correction.json")))
    ref_app = json.load(open(os.path.join(base, "results", "correction",
                                          "applied_correction_summary.json")))
    ref_sig = json.load(open(os.path.join(base, "results", "correction",
                                          "signal_bound.json")))
    out = os.path.join(base, "results", "panel", "_selftest_P4")
    s = run_pair("P4", "ASD", "SCZ", 0.02756, audit, out, do_null_sign=False)

    checks = [
        ("n reported loci", s["n_reported"], 78, 0),
        ("mean rho_hat", s["mean_rho_hat"], ref["mean_rho_hat_observed"], 1e-9),
        ("aggregate slope", s["aggregate_slope"], ref_app["aggregate_slope"], 1e-9),
        ("aggregate intercept", s["fitted_aggregate_intercept"],
         ref_app["aggregate_intercept"], 1e-9),
        ("slope fold-change over binning",
         s["diagnostic"]["C1_stability"]["slope_fold_change"],
         ref["slope_fold_change_over_binning"], 1e-9),
        ("ASD frac of passes expected null",
         s["signal"]["ASD"]["fraction_of_passes_expected_null"],
         ref_sig["ASD"]["fraction_of_passes_expected_null"], 1e-12),
        ("mean CI width", s["mean_ci_width"], ref_app["mean_ci_width"], 1e-9),
    ]
    print("\n=== SELFTEST: generic pipeline against published P4 ===")
    ok = True
    for name, got, want, tol in checks:
        good = abs(got - want) <= tol
        ok &= good
        print(f"  [{'ok ' if good else 'FAIL'}] {name:36s} got {got!r:>22}  want {want!r}")
    print(f"  [{'ok ' if not s['diagnostic']['usable'] else 'FAIL'}] "
          f"{'P4 verdict is NOT APPLICABLE':36s} "
          f"failed={s['diagnostic']['failed_criteria']}")
    ok &= not s["diagnostic"]["usable"]
    print("SELFTEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--pair"); ap.add_argument("--trait1"); ap.add_argument("--trait2")
    ap.add_argument("--r-e", type=float, dest="r_e")
    ap.add_argument("--audit"); ap.add_argument("--outdir")
    ap.add_argument("--null-sign", action="store_true")
    ap.add_argument("--no-part-b", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    for req in ("pair", "trait1", "trait2", "r_e", "audit", "outdir"):
        if getattr(a, req) is None:
            ap.error(f"--{req.replace('_','-')} is required")
    run_pair(a.pair, a.trait1, a.trait2, a.r_e, a.audit, a.outdir,
             a.null_sign, not a.no_part_b)
