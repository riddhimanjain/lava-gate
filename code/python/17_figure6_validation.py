"""17_figure6_validation.py -- Figure 6: the closed form against measurement, and"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import _figstyle as F

HERE = os.path.dirname(os.path.abspath(__file__))
SUB = os.path.abspath(os.path.join(HERE, "..", ".."))
VER = os.path.join(SUB, "results", "verification")
OUT = os.path.join(SUB, "figures", "Figure6_validation")

F.apply()
v4 = pd.read_csv(os.path.join(VER, "V4_closedform_vs_tierB.csv"))
v1 = pd.read_csv(os.path.join(VER, "V1_summary.csv"))
ORDER = ["P4", "P1", "P2a", "P2b", "P3", "P5"]
v4 = v4.set_index("pair").loc[ORDER].reset_index()
v1 = v1.set_index("pair").loc[ORDER].reset_index()

fig, axes = plt.subplots(1, 3, figsize=(3 * F.PANEL_W, F.PANEL_H))

ax = axes[0]
lo, hi = 0.54, 0.96
F.refline(ax, 0, axis="y", kind="identity")
ax.lines[-1].remove()
ax.plot([lo, hi], [lo, hi], color=F.BLACK, ls="--", lw=F.LW_REF,
        zorder=2, label="exact agreement")
ax.scatter(v4.cf_slope_mean, v4.tierB_slope_D, s=60, color=F.C_BLUE,
           alpha=0.9, zorder=3, edgecolor="white", linewidth=0.8)
for _, r in v4.iterrows():
    ax.annotate(r["pair"], (r.cf_slope_mean, r.tierB_slope_D),
                textcoords="offset points", xytext=(7, -3), fontsize=9)
worst = (v4.cf_slope_mean - v4.tierB_slope_D).abs().max()
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("Closed-form slope $S_{locus}$ (theory)")
ax.set_ylabel("Measured slope (Tier B, LAVA's code path)")
ax.set_title("Closed form against measurement")
ax.legend(loc="upper left")
ax.text(0.97, 0.05, f"max |error| = {worst:.3f}\nno parameter fitted to the panel",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=F.C_GREY, alpha=0.9))

ax = axes[1]
x = np.arange(len(ORDER))
ax.plot(x, v1.slope_B_estimable_trunc, "o--", color=F.C_GREY, ms=F.MS_SERIES,
        alpha=F.A_QUIET, label="no filter (estimable only)")
ax.plot(x, v1.slope_C_gated_raw, "o-", color=F.C_BLUE, ms=F.MS_SERIES,
        lw=F.LW_EMPH, label="filter, no truncation")
ax.plot(x, v1.slope_D_gated_trunc, "o-", color=F.C_ORANGE, ms=F.MS_SERIES,
        label="filter + `cap.estimates` (LAVA)")
for i in x:
    ax.annotate("", xy=(i, v1.slope_D_gated_trunc[i]),
                xytext=(i, v1.slope_C_gated_raw[i]),
                arrowprops=dict(arrowstyle="-", color=F.C_ORANGE, lw=0.8,
                                alpha=0.5))
F.refline(ax, 1.0, axis="y", kind="identity")
ax.set_xticks(x); ax.set_xticklabels(ORDER)
ax.set_ylim(0.5, 1.05)
ax.set_ylabel("Slope of $E[\\hat{\\rho}]$ on $\\rho$")
ax.set_xlabel("Trait pair, ordered as in Table 3")
ax.set_title("What the attenuation is caused by")
ax.legend(loc="lower right")
gap = (v1.slope_C_gated_raw - v1.slope_D_gated_trunc).abs()
ax.text(0.03, 0.05,
        f"truncation accounts for\n{gap.min():.3f}–{gap.max():.3f} of the slope",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=F.C_GREY, alpha=0.9))

ax = axes[2]
ax.scatter(v4.r_e, v4.rule_0p68_intercept, s=F.S_SCATTER * 2, color=F.C_RED,
           marker="^", alpha=F.A_SCATTER, zorder=3, label="rule: 0.68 × $r_e$")
ax.scatter(v4.r_e, v4.cf_intercept_mean, s=F.S_SCATTER * 2, color=F.C_BLUE,
           alpha=F.A_SCATTER, zorder=3, label="closed form $I_{locus}$")
ax.scatter(v4.r_e, v4.tierB_intercept_D, s=F.S_SCATTER * 2.4, color=F.BLACK,
           marker="x", zorder=4, label="measured (Tier B)")
OFF = {"P4": (7, -4), "P1": (-26, -4), "P2a": (4, -13), "P2b": (7, 2),
       "P3": (-28, -4), "P5": (6, -12)}
for _, r in v4.iterrows():
    ax.annotate(r["pair"], (r.r_e, r.tierB_intercept_D),
                textcoords="offset points", xytext=OFF[r["pair"]], fontsize=9)
    ax.plot([r.r_e, r.r_e], [r.tierB_intercept_D, r.rule_0p68_intercept],
            color=F.C_RED, lw=0.7, alpha=0.4, zorder=1)
ax.set_xlabel("Sample-overlap correlation $r_e$")
ax.set_ylabel("Manufactured intercept at $\\rho = 0$")
ax.set_title("Why a constant rule fails")
ax.legend(loc="upper left")
ax.set_ylim(0.0, 0.185)
mae_rule = (v4.rule_0p68_intercept - v4.tierB_intercept_D).abs().mean()
mae_cf = (v4.cf_intercept_mean - v4.tierB_intercept_D).abs().mean()
ax.text(0.50, 0.97, f"MAE   rule {mae_rule:.4f}\n         closed form {mae_cf:.4f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=F.C_GREY, alpha=0.9))

F.save(fig, OUT)
print(f"  panel a: max closed-form error {worst:.4f}")
print(f"  panel b: truncation share {gap.min():.4f}-{gap.max():.4f}")
print(f"  panel c: MAE rule {mae_rule:.4f} vs closed form {mae_cf:.4f}")
