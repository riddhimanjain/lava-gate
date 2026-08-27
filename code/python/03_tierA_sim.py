"""03_tierA_sim.py — Tier A Monte Carlo validation of the gate-bias closed form."""
import json
import numpy as np
import pandas as pd
from scipy import stats
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()

N_REF = 100_000
ALPHA = 0.05
PARAM_LIM = 1.25
N_LOCI = 200_000
BATCH = 20_000
SEED = 20260814


def simulate_cell(K, theta1, theta2, rho, r_e, n_loci=N_LOCI, rng=None):
    """Draw n_loci independent loci and return the estimator outputs."""
    rng = rng or np.random.default_rng(SEED)
    nref = A.nref_scale(K, N_REF)
    t_crit = A.gate_threshold(K, ALPHA)

    Sigma = np.array([[1.0, r_e], [r_e, 1.0]])
    L_sig = np.linalg.cholesky(Sigma)

    Om = np.array([[theta1, rho * np.sqrt(theta1 * theta2)],
                   [rho * np.sqrt(theta1 * theta2), theta2]])
    L_om = np.linalg.cholesky(Om + 1e-15 * np.eye(2))
    M1 = L_om.T * np.sqrt(K)
    assert np.allclose(M1.T @ M1, K * Om, atol=1e-8), "mean matrix construction failed"

    m = K - 2
    out = {k: [] for k in ("rho_hat", "passed", "h2_1", "h2_2")}
    done = 0
    while done < n_loci:
        n = min(BATCH, n_loci - done)
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
        S = (M1.T @ M1)[None, :, :] + cross + np.swapaxes(cross, 1, 2) \
            + L_sig @ ZtZ @ L_sig.T

        stat1 = S[:, 0, 0] / Sigma[0, 0] * nref
        stat2 = S[:, 1, 1] / Sigma[1, 1] * nref
        passed = (stat1 > t_crit) & (stat2 > t_crit)

        om = (S / K - Sigma[None, :, :]) * nref
        o11, o22, o12 = om[:, 0, 0], om[:, 1, 1], om[:, 0, 1]
        with np.errstate(invalid="ignore", divide="ignore"):
            rh = o12 / np.sqrt(o11 * o22)

        out["rho_hat"].append(rh)
        out["passed"].append(passed)
        out["h2_1"].append(o11)
        out["h2_2"].append(o22)
        done += n

    res = {k: np.concatenate(v) for k, v in out.items() if k != "dropped"}
    return res


def summarize(res, rho):
    """Apply LAVA's reporting rule and summarise, mirroring run.bivar()."""
    rh, passed = res["rho_hat"], res["passed"]
    gated = rh[passed]
    n_gate = gated.size
    finite = np.isfinite(gated)
    keep = finite & (np.abs(gated) <= PARAM_LIM)
    reported = np.clip(gated[keep], -1.0, 1.0)
    return dict(
        rho=rho,
        n_gate_pass=int(n_gate),
        gate_rate=float(passed.mean()),
        n_nonfinite=int((~finite).sum()),
        n_paramlim_dropped=int((finite & (np.abs(gated) > PARAM_LIM)).sum()),
        paramlim_drop_frac=float((finite & (np.abs(gated) > PARAM_LIM)).sum() / max(n_gate, 1)),
        n_reported=int(reported.size),
        mean_rho_hat=float(reported.mean()) if reported.size else np.nan,
        se_mean=float(reported.std(ddof=1) / np.sqrt(reported.size)) if reported.size > 1 else np.nan,
        median_rho_hat=float(np.median(reported)) if reported.size else np.nan,
        mean_rho_hat_ungated=float(np.nanmean(rh[np.isfinite(rh) & (np.abs(rh) <= PARAM_LIM)])),
        mean_h2_1_gated=float(res["h2_1"][passed].mean()),
        mean_h2_2_gated=float(res["h2_2"][passed].mean()),
    )


def main():
    preds = pd.read_csv("out/predictions.csv")
    scen = (preds[["scenario", "K", "theta1", "theta2", "r_e"]]
            .drop_duplicates().set_index("scenario"))
    rho_grid = sorted(preds["rho"].unique())

    rng = np.random.default_rng(SEED)
    rows = []
    for name, cfg in scen.iterrows():
        for rho in rho_grid:
            res = simulate_cell(int(cfg.K), cfg.theta1, cfg.theta2, rho, cfg.r_e, rng=rng)
            s = summarize(res, rho)
            s["scenario"] = name
            rows.append(s)
            print(f"{name:24s} rho={rho:.1f}  gate={s['gate_rate']:.4f}  "
                  f"mean_rho_hat={s['mean_rho_hat']:.4f}  "
                  f"paramlim_drop={s['paramlim_drop_frac']:.4f}", flush=True)
    sim = pd.DataFrame(rows)
    sim.to_csv("out/tierA_sim.csv", index=False)

    cmp = sim.merge(preds[["scenario", "rho", "E_rho_given_G", "attenuation",
                           "a1", "a2", "p_gate1", "p_gate2", "p_gate_joint_indep"]],
                    on=["scenario", "rho"])
    cmp["resid"] = cmp["mean_rho_hat"] - cmp["E_rho_given_G"]
    cmp["z"] = cmp["resid"] / cmp["se_mean"]
    cmp.to_csv("out/tierA_vs_predictions.csv", index=False)

    fits = {}
    for name in scen.index:
        d = cmp[cmp.scenario == name]
        sl, ic = np.polyfit(d["rho"], d["mean_rho_hat"], 1)
        resid = d["mean_rho_hat"] - (sl * d["rho"] + ic)
        pred_sl = float(np.polyfit(d["rho"], d["E_rho_given_G"], 1)[0])
        fits[name] = dict(
            sim_slope=float(sl), sim_intercept=float(ic),
            sim_max_abs_resid_linear=float(np.abs(resid).max()),
            pred_slope=pred_sl,
            attenuation_1_over_sqrt_a1a2=float(d["attenuation"].iloc[0]),
            slope_vs_pred_reldiff=float(abs(sl - pred_sl) / pred_sl),
            slope_vs_attenuation_reldiff=float(
                abs(sl - d["attenuation"].iloc[0]) / d["attenuation"].iloc[0]),
            max_abs_resid_vs_pred=float(d["resid"].abs().max()),
            max_abs_z=float(d["z"].abs().max()),
        )
    with open("out/tierA_fits.json", "w") as f:
        json.dump(fits, f, indent=2)
    print("\n=== SLOPE COMPARISON ===")
    for k, v in fits.items():
        print(f"{k:24s} sim={v['sim_slope']:.4f} pred={v['pred_slope']:.4f} "
              f"1/sqrt(a1a2)={v['attenuation_1_over_sqrt_a1a2']:.4f} "
              f"reldiff_pred={v['slope_vs_pred_reldiff']:.3%} "
              f"linres={v['sim_max_abs_resid_linear']:.5f}")


if __name__ == "__main__":
    main()
