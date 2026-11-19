"""13_figure5.py -- Figure 5: the correction, and the diagnostic that says when it"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from importlib import machinery

FS = machinery.SourceFileLoader("figstyle", "_figstyle.py").load_module()
FS.apply()

val = pd.read_csv("out/correction_validation.csv")
csum = json.load(open("out/correction_summary.json"))
agg = json.load(open("out/aggregate_correction.json"))
flr = pd.read_csv("out/floor_diagnostic.csv")

fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))

ax = axes[0]
labels = ["Uncorrected", "Single pooled\nfactor", "Gate-aware,\nper locus",
          "Gate-aware,\naggregate"]
vals = [csum["mean_abs_bias_raw"], csum["mean_abs_bias_pooled"],
        csum["mean_abs_bias_corrected"], csum["mean_abs_bias_corrected_aggregate"]]
cols = [FS.C_GREY, FS.C_ORANGE, FS.C_BLUE, FS.C_GREEN]
bars = ax.bar(range(4), vals, color=cols, edgecolor="black", linewidth=0.6, width=0.62)
for i, (b, v) in enumerate(zip(bars, vals)):
    ax.text(b.get_x() + b.get_width()/2, v + 0.006, f"{v:.3f}",
            ha="center", va="bottom", fontsize=9)
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel(r"mean absolute bias in $\rho$")
ax.set_ylim(0, max(vals) * 1.18)
ax.set_title(f"Bias by estimator ({csum['n_cells']} simulation cells)")
ax.grid(alpha=FS.A_GRID, axis="y")

ax = axes[1]
cell = "V1_ASDlike_SCZlike"
d = val[val.cell == cell].sort_values("rho_true")
ax.plot(d.rho_true, d.mean_raw, "o-", ms=FS.MS_SERIES, color=FS.C_GREY,
        label="uncorrected")
ax.plot(d.rho_true, d.mean_pooled, "o-", ms=FS.MS_SERIES, color=FS.C_ORANGE,
        label="single pooled factor")
ax.plot(d.rho_true, d.mean_corrected_aggregate, "o-", ms=FS.MS_SERIES,
        color=FS.C_GREEN, label="gate-aware, aggregate")
lim = [-0.75, 0.95]
ax.plot(lim, lim, "k--", lw=1, label="identity")
ax.set_xlim(*lim); ax.set_ylim(*lim)
ax.set_xlabel(r"true $\rho$"); ax.set_ylabel(r"recovered $\rho$")
ax.set_title(f"Recovery, K={int(d.K.iloc[0])}, "
             r"$r_e$=" + f"{d.r_e.iloc[0]:.4f}")
ax.legend(fontsize=FS.FS_LEGEND, frameon=True, loc="upper left")
ax.grid(alpha=FS.A_GRID)

ax = axes[2]
for trait, col, mk in (("ASD", FS.C_RED, "o"), ("SCZ", FS.C_BLUE, "s")):
    g = flr[flr.trait == trait].sort_values("K_bin")
    ax.scatter(g.K_bin, g.z_above_null, s=FS.S_SCATTER * 1.6, color=col,
               marker=mk, alpha=FS.A_SCATTER, edgecolor="black", linewidth=0.4,
               label=trait, zorder=3)
ax.axhline(0, color="black", lw=0.8, zorder=2)
ax.axhspan(-2, 2, color=FS.BAND, zorder=1)
ax.text(320, 0.35, "indistinguishable from a locus\nwith no true local signal",
        fontsize=7.5, ha="center", color="#444444")
ax.set_yscale("symlog", linthresh=5)
ax.set_xlabel("K (binned)")
ax.set_ylabel(r"$z$ of observed $\hat{\theta}$ above the gated null")
ax.set_title("Applicability diagnostic, 78 reported loci")
ax.legend(fontsize=FS.FS_LEGEND, frameon=True, loc="upper right")
ax.grid(alpha=FS.A_GRID)

fig.tight_layout()
FS.save(fig, "out/figures/Figure5_correction")
print("wrote Figure5_correction")
print("  panel a values:", dict(zip(labels, vals)))
print("  panel c ASD z range:", flr[flr.trait=='ASD'].z_above_null.min(),
      flr[flr.trait=='ASD'].z_above_null.max())
