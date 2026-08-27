"""02_predictions.py — Generate PRE-REGISTERED predictions from the closed form."""
import json
import numpy as np
import pandas as pd
from importlib import machinery
from datetime import datetime, timezone

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()

N_REF = 100_000
ALPHA = 0.05

SCENARIOS = {
    "S1_asdlike_x_sczlike": dict(K=300, theta1=0.065, theta2=0.220, r_e=0.0),
    "S2_matched_moderate":  dict(K=300, theta1=0.100, theta2=0.100, r_e=0.0),
    "S3_overlap_re05":      dict(K=300, theta1=0.065, theta2=0.220, r_e=0.05),
    "S4_overlap_re10":      dict(K=300, theta1=0.065, theta2=0.220, r_e=0.10),
    "S5_smallK":            dict(K=100, theta1=0.065, theta2=0.220, r_e=0.0),
    "S6_largeK":            dict(K=500, theta1=0.065, theta2=0.220, r_e=0.0),
}
RHO_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def main():
    rows = []
    for name, cfg in SCENARIOS.items():
        for rho in RHO_GRID:
            out = A.predict_rho_given_gate(rho=rho, n_ref=N_REF, alpha=ALPHA, **cfg)
            out["scenario"] = name
            rows.append(out)
    df = pd.DataFrame(rows)
    df.to_csv("out/predictions.csv", index=False)

    fits = {}
    for name in SCENARIOS:
        d = df[df.scenario == name]
        sl, ic = np.polyfit(d["rho"], d["E_rho_given_G"], 1)
        resid = d["E_rho_given_G"] - (sl * d["rho"] + ic)
        fits[name] = dict(slope=float(sl), intercept=float(ic),
                          max_abs_resid=float(np.abs(resid).max()))
    with open("out/predictions_fits.json", "w") as f:
        json.dump(fits, f, indent=2)

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    L = []
    L.append("# PRE-REGISTERED PREDICTIONS — gate-induced bias in LAVA local rg\n")
    L.append(f"**Generated:** {ts}  \n**Source:** closed form in `01_analytic.py`, "
             "no simulation involved.\n")
    L.append("> Frozen before `03_tierA_sim.py` was written or run. If the simulation "
             "disagrees with these numbers, the disagreement is the result and gets "
             "reported as such. This file is not to be edited afterwards.\n")

    L.append("\n## Parameters\n")
    L.append(f"- Reference panel N = {N_REF:,} (from LD panel `.info` NOBS ~ 99,339)")
    L.append(f"- Univariate gate: chi-square, alpha = {ALPHA} (both traits must pass)")
    L.append("- theta_i = Omega_ii/Sigma_ii, the local signal-to-noise ratio")
    L.append("- param.lim = 1.25, cap.estimates = TRUE (applied in simulation, not here)\n")

    L.append("\n## Falsifiable claims\n")
    L.append("**P1 — Direction.** E[rho_hat | G] is attenuated toward zero: "
             "|E[rho_hat|G]| < |rho| for every scenario with r_e = 0.\n")
    L.append("**P2 — Mechanism is the denominator.** Attenuation equals "
             "1/sqrt(a_1*a_2) to leading order, where a_i = E[omega_hat_ii|G]/Omega_ii "
             "is the conditional inflation of local h2. The numerator plays only a "
             "second-order role.\n")
    L.append("**P3 — Near-multiplicativity.** Because a_i depends on Omega_ii but not "
             "on Omega_12, the attenuation factor is (to leading order) independent of "
             "rho. So E[rho_hat|G] should be close to LINEAR through the origin in rho, "
             "with residuals from a straight-line fit below 0.01 in absolute value when "
             "r_e = 0.\n")
    L.append("**P4 — Sample overlap creates an intercept.** With r_e > 0, "
             "E[rho_hat|G] != 0 at rho = 0. The intercept grows with r_e.\n")
    L.append("**P5 — K dependence.** Attenuation worsens (slope falls) as K grows at "
             "fixed theta, because the gate threshold tightens relative to the signal.\n")

    L.append("\n## Predicted slope and intercept of E[rho_hat|G] against rho\n")
    L.append("| Scenario | K | theta1 | theta2 | r_e | a1 | a2 | 1/sqrt(a1 a2) | "
             "pred slope | pred intercept | max resid |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for name, cfg in SCENARIOS.items():
        d = df[df.scenario == name].iloc[0]
        f = fits[name]
        L.append(f"| {name} | {cfg['K']} | {cfg['theta1']} | {cfg['theta2']} | "
                 f"{cfg['r_e']} | {d['a1']:.3f} | {d['a2']:.3f} | "
                 f"{d['attenuation']:.4f} | {f['slope']:.4f} | {f['intercept']:.4f} | "
                 f"{f['max_abs_resid']:.5f} |")

    L.append("\n## Predicted gate pass rates\n")
    L.append("| Scenario | P(trait1 passes) | P(trait2 passes) | P(both, indep approx) |")
    L.append("|---|---|---|---|")
    for name in SCENARIOS:
        d = df[df.scenario == name].iloc[0]
        L.append(f"| {name} | {d['p_gate1']:.4f} | {d['p_gate2']:.4f} | "
                 f"{d['p_gate_joint_indep']:.4f} |")

    L.append("\n## Full predicted curve\n")
    L.append("| Scenario | rho | predicted E[rho_hat \\| G] |")
    L.append("|---|---|---|")
    for name in SCENARIOS:
        for _, r in df[df.scenario == name].iterrows():
            L.append(f"| {name} | {r['rho']:.1f} | {r['E_rho_given_G']:.4f} |")

    L.append("\n## What would falsify this\n")
    L.append("- Simulated E[rho_hat|G] *above* rho (over- rather than under-estimation) "
             "at r_e = 0 kills P1 and the paper's framing with it.")
    L.append("- Straight-line residuals above 0.01 kill P3, and the bias would then "
             "have to be modelled as genuinely rho-dependent.")
    L.append("- Simulated slope differing from 1/sqrt(a1 a2) by more than ~10% relative "
             "means the denominator-inflation account (P2) is incomplete.")

    with open("out/PREDICTIONS_PREREGISTERED.md", "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    print("\n".join(L[:40]))
    print("\n... written to out/PREDICTIONS_PREREGISTERED.md")


if __name__ == "__main__":
    main()
