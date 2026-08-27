"""06_coverage.py — Deliverable D: coverage AND width of LAVA's own intervals,"""
import json
import numpy as np
import pandas as pd
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()

N_REF = 100_000
ALPHA = 0.05
PARAM_LIM = 1.25
N_BOOT = 10_000
SEED = 606


def _wishart2(nu, Sigma, Theta, n, rng):
    """Sample n draws from noncentral Wishart_2(nu, Sigma, Theta)."""
    L = np.linalg.cholesky(Sigma)
    ev = np.linalg.eigvalsh(Theta)
    Th = Theta if ev.min() > 0 else Theta + (abs(ev.min()) + 1e-12) * np.eye(2)
    M1 = np.linalg.cholesky(Th).T

    m = nu - 2
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
    cross = MtZ @ L.T
    return (M1.T @ M1)[None] + cross + np.swapaxes(cross, 1, 2) + L @ ZtZ @ L.T


def lava_ci(omega_hat, Sigma, K, rng, n_boot=N_BOOT):
    """Replicate ci.bivariate() for a single locus. Returns (lower, upper)."""
    d = np.diag(omega_hat).copy()
    if np.any(d <= 0) or not np.all(np.isfinite(omega_hat)):
        return np.nan, np.nan
    s = np.sqrt(d)
    corr = omega_hat / np.outer(s, s)
    corr = np.clip(corr, -0.99999, 0.99999)
    np.fill_diagonal(corr, 1.0)
    om = np.outer(s, s) * corr

    draws = _wishart2(K, Sigma / K, om, n_boot, rng)
    o = draws - Sigma[None]
    den = o[:, 0, 0] * o[:, 1, 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        r = o[:, 0, 1] / np.sqrt(den)
    r = r[np.isfinite(r)]
    if r.size < 100:
        return np.nan, np.nan
    lo, hi = np.quantile(r, [0.025, 0.975])
    return float(np.clip(lo, -1, 1)), float(np.clip(hi, -1, 1))


def run_cell(K, theta1, theta2, rho, r_e, n_loci, rng, max_ci=1500):
    nref = A.nref_scale(K, N_REF)
    t_crit = A.gate_threshold(K, ALPHA)
    Sigma = np.array([[1.0, r_e], [r_e, 1.0]])
    L_sig = np.linalg.cholesky(Sigma)
    Om = np.array([[theta1, rho * np.sqrt(theta1 * theta2)],
                   [rho * np.sqrt(theta1 * theta2), theta2]])
    M1 = np.linalg.cholesky(Om + 1e-15 * np.eye(2)).T * np.sqrt(K)

    m = K - 2
    Z1 = rng.standard_normal((n_loci, 2, 2))
    a11 = np.sqrt(rng.chisquare(m, n_loci))
    a22 = np.sqrt(rng.chisquare(m - 1, n_loci))
    a21 = rng.standard_normal(n_loci)
    V = np.empty((n_loci, 2, 2))
    V[:, 0, 0] = a11 ** 2
    V[:, 0, 1] = V[:, 1, 0] = a11 * a21
    V[:, 1, 1] = a21 ** 2 + a22 ** 2
    ZtZ = np.einsum("nki,nkj->nij", Z1, Z1) + V
    MtZ = np.einsum("ik,nkj->nij", M1.T, Z1)
    cross = MtZ @ L_sig.T
    S = (M1.T @ M1)[None] + cross + np.swapaxes(cross, 1, 2) + L_sig @ ZtZ @ L_sig.T

    passed = ((S[:, 0, 0] / Sigma[0, 0] * nref > t_crit) &
              (S[:, 1, 1] / Sigma[1, 1] * nref > t_crit))
    omega = (S / K - Sigma[None]) * nref

    out = {}
    for label, idx in (("post_gate", np.flatnonzero(passed)),
                       ("pre_gate", np.arange(n_loci))):
        if idx.size == 0:
            continue
        take = idx if idx.size <= max_ci else rng.choice(idx, max_ci, replace=False)
        cov, wid, vac, nan_ci, rep = 0, [], 0, 0, 0
        for j in take:
            oh = omega[j]
            rh = oh[0, 1] / np.sqrt(oh[0, 0] * oh[1, 1]) if oh[0, 0] > 0 and oh[1, 1] > 0 else np.nan
            if not np.isfinite(rh) or abs(rh) > PARAM_LIM:
                continue
            rep += 1
            lo, hi = lava_ci(oh, Sigma, K, rng)
            if not np.isfinite(lo):
                nan_ci += 1
                continue
            wid.append(hi - lo)
            if lo <= rho <= hi:
                cov += 1
            if (hi - lo) > 1.0:
                vac += 1
        n_ok = len(wid)
        out[label] = dict(
            n_attempted=int(take.size), n_reported=int(rep), n_ci_ok=int(n_ok),
            n_ci_failed=int(nan_ci),
            coverage=float(cov / n_ok) if n_ok else np.nan,
            mean_width=float(np.mean(wid)) if n_ok else np.nan,
            median_width=float(np.median(wid)) if n_ok else np.nan,
            frac_width_gt_1=float(vac / n_ok) if n_ok else np.nan,
        )
    out["gate_rate"] = float(passed.mean())
    return out


def main():
    rng = np.random.default_rng(SEED)
    SCEN = {
        "ASDlike_x_SCZlike": dict(K=300, theta1=0.065, theta2=0.220, r_e=0.0),
        "matched_moderate":  dict(K=300, theta1=0.100, theta2=0.100, r_e=0.0),
        "with_overlap_r10":  dict(K=300, theta1=0.065, theta2=0.220, r_e=0.10),
    }
    RHOS = [0.0, 0.3, 0.6]
    rows = []
    for name, cfg in SCEN.items():
        for rho in RHOS:
            r = run_cell(rho=rho, n_loci=30_000, rng=rng, **cfg)
            for label in ("pre_gate", "post_gate"):
                if label in r:
                    rows.append(dict(scenario=name, rho=rho, gate=label,
                                     gate_rate=r["gate_rate"], **cfg, **r[label]))
            print(f"{name:20s} rho={rho:.1f}  "
                  f"pre: cov={r['pre_gate']['coverage']:.3f} w={r['pre_gate']['mean_width']:.3f} | "
                  f"post: cov={r['post_gate']['coverage']:.3f} w={r['post_gate']['mean_width']:.3f}",
                  flush=True)
    df = pd.DataFrame(rows)
    df.to_csv("out/coverage.csv", index=False)

    summ = {}
    for label in ("pre_gate", "post_gate"):
        d = df[df.gate == label]
        summ[label] = dict(
            mean_coverage=float(d.coverage.mean()),
            min_coverage=float(d.coverage.min()),
            mean_width=float(d.mean_width.mean()),
            mean_frac_width_gt_1=float(d.frac_width_gt_1.mean()),
        )
    with open("out/coverage_summary.json", "w") as f:
        json.dump(summ, f, indent=2)
    print("\n=== SUMMARY ===")
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
