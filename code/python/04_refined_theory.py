"""04_refined_theory.py — Second-order closed form, after Tier A falsified the"""
import numpy as np
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()

N_REF = 100_000
ALPHA = 0.05


def delta_i(K, theta, n_ref=N_REF, alpha=ALPHA):
    """Delta = E[S_ii | S_ii > t] - E[S_ii], exact via noncentral chi-square."""
    nref = A.nref_scale(K, n_ref)
    t = A.gate_threshold(K, alpha) / nref
    lam = K * theta
    p_tail, e_trunc = A.ncx2_upper_moments(K, lam, t)
    e_cond = e_trunc / p_tail
    return e_cond - K * (1.0 + theta), p_tail


def predict_v2(rho, K, theta1, theta2, r_e=0.0, n_ref=N_REF, alpha=ALPHA):
    """Second-order closed form for E[rho_hat | G]."""
    Om12 = rho * np.sqrt(theta1 * theta2)
    D1, p1 = delta_i(K, theta1, n_ref, alpha)
    D2, p2 = delta_i(K, theta2, n_ref, alpha)

    b1 = (r_e + Om12 + r_e * theta1) / (1.0 + 2.0 * theta1)
    b2 = (r_e + Om12 + r_e * theta2) / (1.0 + 2.0 * theta2)

    num = Om12 + (b1 * D1 + b2 * D2) / K
    d1 = theta1 + D1 / K
    d2 = theta2 + D2 / K
    pred = num / np.sqrt(d1 * d2)

    return {
        "rho": rho, "K": K, "theta1": theta1, "theta2": theta2, "r_e": r_e,
        "Delta1": D1, "Delta2": D2, "b1": b1, "b2": b2,
        "E_omega12_given_G": num, "E_omega11_given_G": d1, "E_omega22_given_G": d2,
        "p_gate1": p1, "p_gate2": p2, "p_gate_joint": p1 * p2,
        "E_rho_given_G_v2": pred,
    }


if __name__ == "__main__":
    import pandas as pd
    import json
    from importlib import machinery as mm
    S3 = mm.SourceFileLoader("sim", "03_tierA_sim.py").load_module()

    HELD_OUT = {
        "H1_lowpower_both":   dict(K=200, theta1=0.040, theta2=0.070, r_e=0.00),
        "H2_asym_extreme":    dict(K=400, theta1=0.030, theta2=0.350, r_e=0.00),
        "H3_overlap_high":    dict(K=250, theta1=0.080, theta2=0.180, r_e=0.20),
        "H4_overlap_mid_bigK":dict(K=600, theta1=0.120, theta2=0.120, r_e=0.08),
        "H5_tinyK":           dict(K=60,  theta1=0.150, theta2=0.250, r_e=0.00),
        "H6_negative_rho":    dict(K=300, theta1=0.065, theta2=0.220, r_e=0.00),
    }
    RHOS = [-0.6, -0.3, 0.0, 0.2, 0.45, 0.7, 0.9]

    rng = np.random.default_rng(99)
    rows = []
    for name, cfg in HELD_OUT.items():
        for rho in RHOS:
            res = S3.simulate_cell(cfg["K"], cfg["theta1"], cfg["theta2"],
                                   rho, cfg["r_e"], n_loci=300_000, rng=rng)
            s = S3.summarize(res, rho)
            v1 = A.predict_rho_given_gate(rho=rho, K=cfg["K"], theta1=cfg["theta1"],
                                          theta2=cfg["theta2"], r_e=cfg["r_e"],
                                          n_ref=N_REF, alpha=ALPHA)
            v2 = predict_v2(rho=rho, **cfg)
            rows.append(dict(
                scenario=name, rho=rho, **{k: cfg[k] for k in cfg},
                sim=s["mean_rho_hat"], se=s["se_mean"],
                sim_gate=s["gate_rate"], pred_gate=v2["p_gate_joint"],
                v1=v1["E_rho_given_G"], v2=v2["E_rho_given_G_v2"],
                paramlim_drop=s["paramlim_drop_frac"],
            ))
            print(f"{name:22s} rho={rho:+.2f} sim={s['mean_rho_hat']:+.4f} "
                  f"v1={v1['E_rho_given_G']:+.4f} v2={v2['E_rho_given_G_v2']:+.4f}",
                  flush=True)

    df = pd.DataFrame(rows)
    df["err_v1"] = df["sim"] - df["v1"]
    df["err_v2"] = df["sim"] - df["v2"]
    df["z_v2"] = df["err_v2"] / df["se"]
    df.to_csv("out/heldout_validation.csv", index=False)

    summ = {
        "n_cells": int(len(df)),
        "v1_mean_abs_err": float(df["err_v1"].abs().mean()),
        "v1_max_abs_err": float(df["err_v1"].abs().max()),
        "v2_mean_abs_err": float(df["err_v2"].abs().mean()),
        "v2_max_abs_err": float(df["err_v2"].abs().max()),
        "v2_max_abs_z": float(df["z_v2"].abs().max()),
        "v2_frac_within_0.01": float((df["err_v2"].abs() < 0.01).mean()),
        "gate_max_reldiff": float(((df.sim_gate - df.pred_gate).abs()
                                   / df.pred_gate).max()),
    }
    with open("out/heldout_summary.json", "w") as f:
        json.dump(summ, f, indent=2)
    print("\n=== HELD-OUT SUMMARY ===")
    for k, v in summ.items():
        print(f"  {k:24s} {v}")
