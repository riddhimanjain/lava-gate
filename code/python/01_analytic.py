"""01_analytic.py — Closed-form machinery for gate-induced bias in LAVA's local rg."""

import numpy as np
from scipy import stats
from scipy.special import gammaln

__all__ = [
    "nref_scale", "gate_threshold", "ncx2_upper_moments",
    "h2_inflation", "predict_rho_given_gate", "gate_pass_prob",
]


def nref_scale(K, n_ref):
    """LAVA input_processing.R:134 — (N_ref - K - 1)/(N_ref - 1)."""
    return (n_ref - K - 1.0) / (n_ref - 1.0)


def gate_threshold(K, alpha=0.05):
    """Critical value of the univariate chi2_K test at level alpha (binary traits)."""
    return stats.chi2.ppf(1.0 - alpha, K)


def _poisson_weights(lam, tol=1e-14, jmax=None):
    """Poisson(lam/2) mixture weights for the noncentral chi-square representation."""
    half = lam / 2.0
    if jmax is None:
        jmax = int(max(60, half + 12.0 * np.sqrt(max(half, 1.0))))
    j = np.arange(jmax + 1)
    logw = -half + j * np.log(half + 1e-300) - gammaln(j + 1.0)
    w = np.exp(logw)
    if w.max() > 0:
        keep = w > tol * w.max()
        j, w = j[keep], w[keep]
    return j, w / w.sum()


def ncx2_upper_moments(K, lam, t):
    """Exact P(U>t) and E[U*1{U>t}] for U ~ noncentral chi2_K(lam)."""
    j, w = _poisson_weights(lam)
    m = K + 2.0 * j
    p_tail = np.sum(w * stats.chi2.sf(t, m))
    e_trunc = np.sum(w * m * stats.chi2.sf(t, m + 2.0))
    return p_tail, e_trunc


def gate_pass_prob(K, theta1, theta2, n_ref=None, alpha=0.05, r_e=0.0):
    """Marginal per-trait gate pass probabilities."""
    nref = 1.0 if n_ref is None else nref_scale(K, n_ref)
    t = gate_threshold(K, alpha) / nref
    p1, _ = ncx2_upper_moments(K, K * theta1, t)
    p2, _ = ncx2_upper_moments(K, K * theta2, t)
    return p1, p2, p1 * p2


def h2_inflation(K, theta, n_ref=None, alpha=0.05):
    """a = E[omega_hat_ii | gate_i] / Omega_true[i,i], the conditional h2 inflation."""
    nref = 1.0 if n_ref is None else nref_scale(K, n_ref)
    t = gate_threshold(K, alpha) / nref
    lam = K * theta
    p_tail, e_trunc = ncx2_upper_moments(K, lam, t)
    if p_tail <= 0:
        return np.nan, np.nan
    e_cond = e_trunc / p_tail
    a = (e_cond / K - 1.0) / theta
    if np.isfinite(a) and a < 1.0 - 1e-8:
        raise AssertionError(
            f"h2_inflation returned a={a:.6f} < 1 for K={K}, theta={theta}. "
            "Conditioning on U>t must inflate E[U]; this indicates a numerical bug."
        )
    return a, p_tail


def predict_rho_given_gate(rho, K, theta1, theta2, n_ref=None, alpha=0.05, r_e=0.0):
    """Leading-order closed form for E[rho_hat | G]."""
    a1, p1 = h2_inflation(K, theta1, n_ref, alpha)
    a2, p2 = h2_inflation(K, theta2, n_ref, alpha)
    atten = 1.0 / np.sqrt(a1 * a2)

    coupling = 0.0
    for (a, th) in ((a1, theta1), (a2, theta2)):
        coupling += 0.5 * (a - 1.0) / (1.0 + K * th / 2.0)
    num_mult = 1.0 + coupling

    denom_scale = np.sqrt(a1 * theta1 * a2 * theta2)
    intercept = 0.0
    if r_e != 0.0 and denom_scale > 0:
        excess = 0.5 * ((a1 - 1.0) * theta1 + (a2 - 1.0) * theta2)
        intercept = r_e * excess / denom_scale / K * 2.0

    pred = rho * atten * num_mult + intercept
    return {
        "rho": rho, "K": K, "theta1": theta1, "theta2": theta2, "r_e": r_e,
        "a1": a1, "a2": a2,
        "p_gate1": p1, "p_gate2": p2, "p_gate_joint_indep": p1 * p2,
        "attenuation": atten, "num_mult": num_mult, "intercept": intercept,
        "E_rho_given_G": pred,
    }


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for m, t in ((10, 5.0), (50, 60.0), (300, 341.0)):
        mc = rng.chisquare(m, 4_000_000)
        emp = np.mean(mc * (mc > t))
        exact = m * stats.chi2.sf(t, m + 2)
        print(f"m={m:4d} t={t:7.2f}  MC={emp:12.5f}  exact={exact:12.5f}  "
              f"reldiff={abs(emp-exact)/exact:.2e}")
