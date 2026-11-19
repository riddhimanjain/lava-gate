"""09_null_sign.py -- The filtered null for the applied sign test."""
import json
import numpy as np
import pandas as pd
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()
C = machinery.SourceFileLoader("corr", "08_correction.py").load_module()
CV = machinery.SourceFileLoader("cov", "06_coverage.py").load_module()

N_REF = 100_000
ALPHA = 0.05
PARAM_LIM = 1.25
R_E = 0.02756
SEED = 20260823

N_SIGN = 6_000
N_DET = 250
N_BOOT_CI = 2_000


def sim_gated(K, theta1, theta2, rho, r_e, n_want, rng, max_tries=60):
    """Draw until n_want loci have passed the gate. Returns (rho_hat, omega)."""
    nref = A.nref_scale(K, N_REF)
    t_crit = A.gate_threshold(K, ALPHA)
    Sigma = np.array([[1.0, r_e], [r_e, 1.0]])
    L_sig = np.linalg.cholesky(Sigma)
    Om = np.array([[theta1, rho * np.sqrt(theta1 * theta2)],
                   [rho * np.sqrt(theta1 * theta2), theta2]])
    M1 = np.linalg.cholesky(Om + 1e-15 * np.eye(2)).T * np.sqrt(K)
    m = K - 2

    rh_all, om_all, tries, drawn = [], [], 0, 0
    have = 0
    batch = max(20000, n_want * 20)
    while have < n_want and tries < max_tries:
        n = batch
        Z1 = rng.standard_normal((n, 2, 2))
        a11 = np.sqrt(rng.chisquare(m, n))
        a22 = np.sqrt(rng.chisquare(m - 1, n))
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
        keep = np.isfinite(rh) & (np.abs(rh) <= PARAM_LIM)
        rh = np.clip(rh[keep], -1.0, 1.0)
        om = om[keep]
        rh_all.append(rh); om_all.append(om)
        have += rh.size
        drawn += n
        tries += 1
    rh = np.concatenate(rh_all)[:n_want] if rh_all else np.array([])
    om = np.concatenate(om_all)[:n_want] if om_all else np.zeros((0, 2, 2))
    return rh, om, drawn


def poisson_binomial_pmf(p):
    """Exact PMF of sum of independent Bernoulli(p_i), by direct convolution."""
    pmf = np.zeros(len(p) + 1)
    pmf[0] = 1.0
    for pi in p:
        pmf[1:] = pmf[1:] * (1 - pi) + pmf[:-1] * pi
        pmf[0] *= (1 - pi)
    return pmf


def pb_two_sided_p(k_obs, p):
    """Two-sided P-value by the method of small probabilities."""
    pmf = poisson_binomial_pmf(p)
    tol = 1e-9 * pmf.max()
    return float(pmf[pmf <= pmf[k_obs] + tol].sum())


def main():
    rng = np.random.default_rng(SEED)
    aud = pd.read_csv("out_R/attrition_audit.csv")
    rep = aud[aud["gate_pass"] == True].copy()
    rep = rep[np.isfinite(rep["rho_hat"])]
    print(f"Reported loci read from attrition_audit.csv: {len(rep)}")

    th1, th2, Ks = [], [], []
    for _, r in rep.iterrows():
        K = int(r["K"])
        t1o = C.theta_obs_from_p(r["p_ASD"], K)
        t2o = C.theta_obs_from_p(r["p_SCZ"], K)
        a, _, _ = C.invert_theta(np.array([t1o]), K)
        b, _, _ = C.invert_theta(np.array([t2o]), K)
        th1.append(max(float(a[0]), 1e-6)); th2.append(max(float(b[0]), 1e-6))
        Ks.append(K)
    rep["K_int"] = Ks; rep["theta_ASD"] = th1; rep["theta_SCZ"] = th2

    n_obs = len(rep)
    n_pos_obs = int((rep["rho_hat"] > 0).sum())
    det_obs = np.isfinite(rep["rho_lower"]) & np.isfinite(rep["rho_upper"]) & \
              ((rep["rho_lower"] > 0) | (rep["rho_upper"] < 0))
    n_det_obs = int(det_obs.sum())
    n_det_pos_obs = int((det_obs & (rep["rho_hat"] > 0)).sum())
    print(f"OBSERVED: {n_pos_obs}/{n_obs} positive; "
          f"{n_det_obs} sign-determinate, {n_det_pos_obs} of those positive")

    RHO_GRID = [0.0, -0.4, 0.37]
    partA = {}
    for rho0 in RHO_GRID:
        p_pos = []
        for i, r in enumerate(rep.itertuples()):
            rh, _, _ = sim_gated(r.K_int, r.theta_ASD, r.theta_SCZ, rho0, R_E,
                                 N_SIGN, rng)
            p_pos.append(float((rh > 0).mean()) if rh.size else np.nan)
        p_pos = np.array(p_pos)
        pval = pb_two_sided_p(n_pos_obs, p_pos)
        partA[f"rho_{rho0}"] = dict(
            mean_p_positive=float(p_pos.mean()),
            min_p_positive=float(p_pos.min()), max_p_positive=float(p_pos.max()),
            expected_n_positive=float(p_pos.sum()),
            observed_n_positive=n_pos_obs, n_loci=n_obs,
            poisson_binomial_two_sided_p=pval,
        )
        print(f"  [A] rho0={rho0:+.2f}: null P(pos)={p_pos.mean():.4f} "
              f"expected={p_pos.sum():.1f}/{n_obs} observed={n_pos_obs} "
              f"PB P={pval:.3e}", flush=True)

    Sigma = np.array([[1.0, R_E], [R_E, 1.0]])
    p_det, p_det_pos, p_pos_given_det = [], [], []
    for i, r in enumerate(rep.itertuples()):
        rh, om, _ = sim_gated(r.K_int, r.theta_ASD, r.theta_SCZ, 0.0, R_E,
                              N_DET, rng)
        nd = npd = 0
        for j in range(rh.size):
            lo, hi = CV.lava_ci(om[j], Sigma, r.K_int, rng, n_boot=N_BOOT_CI)
            if not np.isfinite(lo):
                continue
            if lo > 0 or hi < 0:
                nd += 1
                if rh[j] > 0:
                    npd += 1
        p_det.append(nd / max(rh.size, 1))
        p_det_pos.append(npd / max(rh.size, 1))
        p_pos_given_det.append(npd / nd if nd else np.nan)
        if (i + 1) % 10 == 0:
            print(f"    [B] {i+1}/{n_obs} loci done", flush=True)

    p_det = np.array(p_det); p_det_pos = np.array(p_det_pos)
    pgd = np.array(p_pos_given_det)
    pooled_pos_given_det = float(np.nansum(p_det_pos) / np.nansum(p_det)) \
        if np.nansum(p_det) > 0 else np.nan

    pval_det = pb_two_sided_p(n_det_obs, p_det)
    p_dir = float(pooled_pos_given_det)
    p_all_pos = p_dir ** n_det_obs
    p_all_neg = (1.0 - p_dir) ** n_det_obs
    pval_dir = float(min(1.0, p_all_pos + (p_all_neg if p_all_neg <= p_all_pos else 0.0)))

    partB = dict(
        n_sim_per_locus=N_DET, n_boot_ci=N_BOOT_CI,
        null_mean_p_determinate=float(p_det.mean()),
        null_expected_n_determinate=float(p_det.sum()),
        observed_n_determinate=n_det_obs,
        poisson_binomial_two_sided_p_determinate=pval_det,
        null_p_positive_given_determinate=pooled_pos_given_det,
        observed_n_determinate_positive=n_det_pos_obs,
        direction_test_p=pval_dir,
        note=("p_positive_given_determinate is pooled as sum(P(det & pos))/sum(P(det)) "
              "across loci, which is the correct weighting when per-locus "
              "determinacy probabilities differ."),
    )
    print("\n=== PART B ===")
    for k, v in partB.items():
        print(f"  {k:44s} {v}")

    rep_out = rep[["locus", "chr", "K_int", "theta_ASD", "theta_SCZ",
                   "p_ASD", "p_SCZ", "rho_hat", "rho_lower", "rho_upper"]].copy()
    rep_out["null_p_positive_rho0"] = np.nan
    rep_out["null_p_determinate_rho0"] = p_det
    rep_out.to_csv("out/null_sign_perlocus.csv", index=False)

    with open("out/null_sign_summary.json", "w") as f:
        json.dump(dict(observed=dict(n_loci=n_obs, n_positive=n_pos_obs,
                                     n_determinate=n_det_obs,
                                     n_determinate_positive=n_det_pos_obs),
                       r_e=R_E, partA=partA, partB=partB), f, indent=2)
    print("\nwrote out/null_sign_summary.json")


if __name__ == "__main__":
    main()
