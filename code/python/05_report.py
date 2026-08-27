"""05_report.py — Build RESULTS_TIERA.md and figures from the saved CSVs."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime, timezone

OUT = "out"
sim = pd.read_csv(f"{OUT}/tierA_sim.csv")
cmp_ = pd.read_csv(f"{OUT}/tierA_vs_predictions.csv")
fits = json.load(open(f"{OUT}/tierA_fits.json"))
ho = pd.read_csv(f"{OUT}/heldout_validation.csv")
hos = json.load(open(f"{OUT}/heldout_summary.json"))
pl = pd.read_csv(f"{OUT}/paramlim_isolation.csv")

_p = pd.read_csv(f"{OUT}/predictions.csv")[
    ["scenario", "K", "theta1", "theta2", "r_e"]].drop_duplicates()
cmp_ = cmp_.merge(_p, on="scenario", how="left")
KMAP = _p.set_index("scenario")["K"].to_dict()

f = lambda x: f"{x:.4f}"

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

a = ax[0]
for name in ["S5_smallK", "S1_asdlike_x_sczlike", "S6_largeK"]:
    d = cmp_[cmp_.scenario == name].sort_values("rho")
    a.plot(d.rho, d.mean_rho_hat, "o-", label=f"K={KMAP[name]}", ms=4)
a.plot([0, 0.9], [0, 0.9], "k--", lw=1, label="unbiased")
a.set_xlabel(r"true $\rho$"); a.set_ylabel(r"E[$\hat\rho$ | gate]")
a.set_title("Gate attenuation, by locus size K"); a.legend(fontsize=8); a.grid(alpha=.3)

a = ax[1]
for name, c in [("S1_asdlike_x_sczlike", "C0"), ("S3_overlap_re05", "C1"),
                ("S4_overlap_re10", "C2")]:
    d = cmp_[cmp_.scenario == name].sort_values("rho")
    re = {"S1_asdlike_x_sczlike": 0.0, "S3_overlap_re05": 0.05,
          "S4_overlap_re10": 0.10}[name]
    a.plot(d.rho, d.mean_rho_hat, "o-", color=c, ms=4, label=f"$r_e$={re}")
a.axhline(0, color="k", lw=.8)
a.set_xlabel(r"true $\rho$"); a.set_ylabel(r"E[$\hat\rho$ | gate]")
a.set_title("Sample overlap induces a nonzero intercept")
a.legend(fontsize=8); a.grid(alpha=.3)

a = ax[2]
a.scatter(ho.v1, ho.sim, s=22, alpha=.75, label="leading order (v1)")
a.scatter(ho.v2, ho.sim, s=22, alpha=.75, label="second order (v2)")
lim = [-0.7, 1.0]
a.plot(lim, lim, "k--", lw=1)
a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("predicted"); a.set_ylabel("simulated")
a.set_title("Held-out validation (42 cells)"); a.legend(fontsize=8); a.grid(alpha=.3)

plt.tight_layout()
plt.savefig(f"{OUT}/FigureA_gate_bias.png", dpi=160)
print("wrote", f"{OUT}/FigureA_gate_bias.png")

ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
L = []
L.append("# Tier A Results — gate-induced bias in LAVA local genetic correlation\n")
L.append(f"**Generated:** {ts} from `out/*.csv`. Predictions were frozen in "
         "`out/PREDICTIONS_PREREGISTERED.md` before any simulation ran.\n")
L.append("**Scope.** Pure sampling-model simulation: validates the derivation against "
         "its own assumptions. It does **not** validate the assumptions. That is Tier B, "
         "which needs the real LD panel and the real LAVA code path.\n")
L.append("> ### ⚠️ Superseded as the headline — read this first\n"
         "> Every slope below is **Tier A's** answer, and Tier A has since been superseded on "
         "the headline number. Tier B (real loci, real per-locus *K*, real Σ, LAVA's own code "
         "path) gives an attenuation slope of **0.5512** (95% CI 0.514–0.589), and Tier A's "
         "**0.6064 falls outside that interval**. Per the rule committed to before either ran, "
         "Tier B wins.\n"
         ">\n"
         "> **The paper's headline attenuation is 0.551, not 0.606.** The numbers here remain "
         "valid as Tier A's own output and are what the Tier A vs Tier B comparison is built "
         "from — but do not quote them as the result. Tier B: `out_R/tierB.txt`.\n"
         ">\n"
         "> One further caveat on the `param.lim` section below: on real data `param.lim` "
         "discarded **zero** loci. Its effect here is a statement about the simulated "
         "configurations, not about ASD × SCZ.\n")

L.append("\n## Scorecard against the pre-registered predictions\n")
L.append("| # | Claim | Verdict |")
L.append("|---|---|---|")
L.append("| P1 | Attenuation toward zero | **Confirmed** |")
L.append("| P2 | Attenuation = 1/sqrt(a1 a2), denominator only | **Falsified** |")
L.append("| P3 | Near-multiplicative in rho | **Confirmed (5/6)** |")
L.append("| P4 | Sample overlap creates an intercept | **Direction confirmed, magnitude badly wrong** |")
L.append("| P5 | Attenuation worsens as K grows | **Falsified — and the prose contradicted my own formula** |")

L.append("\n## P1 — Attenuation is real and large\n")
L.append("| Scenario | K | fitted slope | interpretation |")
L.append("|---|---|---|---|")
for name, v in fits.items():
    L.append(f"| {name} | {KMAP[name]} | {f(v['sim_slope'])} | true rho=0.5 reports as "
             f"{0.5*v['sim_slope']:.3f} |")
L.append("\nAt the ASD-like x SCZ-like configuration the slope is "
         f"**{f(fits['S1_asdlike_x_sczlike']['sim_slope'])}** — roughly "
         f"{(1-fits['S1_asdlike_x_sczlike']['sim_slope'])*100:.0f}% of the true signal "
         "is lost before it is ever reported.\n")

L.append("\n## P2 — Falsified: the denominator account is incomplete\n")
L.append("| Scenario | sim slope | 1/sqrt(a1 a2) | rel. error |")
L.append("|---|---|---|---|")
for name, v in fits.items():
    L.append(f"| {name} | {f(v['sim_slope'])} | "
             f"{f(v['attenuation_1_over_sqrt_a1a2'])} | "
             f"{v['slope_vs_attenuation_reldiff']:.1%} |")
L.append("\nErrors of 8–33% mean the gate inflates the **numerator** too, partially "
         "offsetting the denominator inflation. The true attenuation is *less* severe "
         "than the leading-order account predicts. `04_refined_theory.py` derives the "
         "second-order form from the exact noncentral Wishart covariance.\n")

L.append("\n## Refined theory, validated on held-out parameters\n")
L.append("Six parameter settings never used in building the refinement, 7 values of rho "
         "each, 300,000 simulated loci per cell.\n")
L.append(f"- leading order (v1): mean abs error **{hos['v1_mean_abs_err']:.4f}**, "
         f"max {hos['v1_max_abs_err']:.4f}")
L.append(f"- second order (v2): mean abs error **{hos['v2_mean_abs_err']:.4f}**, "
         f"max {hos['v2_max_abs_err']:.4f}")
L.append(f"- {hos['v2_frac_within_0.01']:.0%} of held-out cells within 0.01\n")
L.append("v2 remains statistically distinguishable from the simulation "
         f"(max |z| = {hos['v2_max_abs_z']:.0f} at these Monte Carlo sizes). It is a "
         "second-order expansion, not an identity, and should be reported as such.\n")

L.append("\n## P4 — Sample overlap manufactures correlation from nothing\n")
L.append("| r_e | true rho | simulated E[rho_hat \\| G] | v1 predicted | v2 predicted |")
L.append("|---|---|---|---|---|")
for _, r in ho[(ho.rho == 0) & (ho.r_e > 0)].iterrows():
    L.append(f"| {r.r_e} | 0 | {f(r['sim'])} | ~0.0002 | {f(r['v2'])} |")
for name in ["S3_overlap_re05", "S4_overlap_re10"]:
    d = cmp_[(cmp_.scenario == name) & (cmp_.rho == 0)].iloc[0]
    re = 0.05 if "05" in name else 0.10
    L.append(f"| {re} | 0 | {f(d['mean_rho_hat'])} | ~0.0002 | — |")
L.append("\n**This is the most consequential single result.** Two traits with *zero* true "
         "local genetic correlation, measured in overlapping samples, report a positive "
         "local correlation after the gate. At r_e = 0.20 the spurious value is "
         "~0.156. The leading-order formula under-predicted this by a factor of ~280; "
         "the second-order form captures it.\n")

L.append("\n## P5 — Falsified, and the error was mine\n")
L.append("The pre-registration prose asserted attenuation *worsens* as K grows. The closed "
         "form it was supposedly describing predicted the opposite, and the simulation "
         "agrees with the formula:\n")
L.append("| K | sim slope |")
L.append("|---|---|")
for name, K in [("S5_smallK", 100), ("S1_asdlike_x_sczlike", 300), ("S6_largeK", 500)]:
    L.append(f"| {K} | {f(fits[name]['sim_slope'])} |")
L.append("\nAt fixed theta, larger K means larger noncentrality (lambda = K*theta), so the "
         "gate selects less aggressively. Recorded as a falsified prediction rather than "
         "silently corrected, because the prose and the mathematics disagreed and only "
         "the simulation revealed it.\n")

L.append("\n## Channel 3: param.lim, tested separately\n")
L.append("**Hypothesis tested and rejected:** that v2's residual error was caused by "
         "`param.lim` truncation. Re-running with the truncation disabled barely moved "
         f"the error (mean |err| {pl.err_vs_gateonly.abs().mean():.5f} without vs "
         f"{pl.err_vs_LAVA.abs().mean():.5f} with). The residual is delta-method "
         "curvature, not channel 3.\n")
L.append("`param.lim` is nonetheless a real and separate distortion:\n")
L.append("| \\|rho\\| | fraction dropped | shift in reported mean |")
L.append("|---|---|---|")
g = pl.assign(absrho=pl.rho.abs()).groupby("absrho").agg(
    drop=("drop_frac", "mean"), eff=("paramlim_effect", "mean"))
for idx, r in g.iterrows():
    L.append(f"| {idx:.2f} | {r['drop']:.4f} | {r['eff']:+.4f} |")
L.append("\nNote the reporting rule is discontinuous: rho_hat in (1, 1.25] is **capped to "
         "1**, while rho_hat > 1.25 becomes **NA**. Two nearly identical estimates receive "
         "opposite treatment.\n")

L.append("\n## What Tier B must test\n")
L.append("Everything above assumes the sampling model is correct. Tier B checks that "
         "assumption against reality:\n")
L.append("1. Real per-locus K from the LD panel — not the assumed values swept here.")
L.append("2. Real per-locus Sigma and the actual ASD/SCZ sample overlap.")
L.append("3. The real LAVA code path, including `process.locus` upstream attrition "
         "(channel 1), which this simulation does not model at all.")
L.append("4. Whether the empirical distribution of theta across the 2,495 blocks puts "
         "ASD x SCZ in the severe-attenuation regime.\n")

with open(f"{OUT}/RESULTS_TIERA.md", "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
print("wrote", f"{OUT}/RESULTS_TIERA.md")
