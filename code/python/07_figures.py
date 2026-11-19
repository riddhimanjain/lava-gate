"""07_figures.py -- the paper's figures."""
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import _figstyle as S

S.apply()

OUT = "out/figures"
os.makedirs(OUT, exist_ok=True)

_SUB = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
audit = pd.read_csv("out_R/attrition_audit.csv")
tierB = pd.read_csv("out_R/tierB.csv")
tierB_ext = pd.read_csv(os.path.join(_SUB, "results", "tierB_extended", "tierB_extended.csv"))
tierB_ext_slopes = pd.read_csv(os.path.join(_SUB, "results", "tierB_extended", "tierB_per_locus_slopes.csv"))
tierB_ext_summary = json.load(open(os.path.join(_SUB, "results", "tierB_extended", "tierB_extended_summary.json")))
cov = pd.read_csv("out/coverage.csv")
sim = pd.read_csv("out/tierA_sim.csv")
signdet = pd.DataFrame(json.load(open("out/sign_determinacy.json")))

proc = audit[audit.processed == True].copy()
rep = audit[audit.rho_hat.notna()].copy()
rep["width"] = rep.rho_upper - rep.rho_lower
rep["determinate"] = np.sign(rep.rho_lower) == np.sign(rep.rho_upper)

pooled_slope = tierB_ext_summary["slope_full"]
pooled_int = tierB_ext_summary["intercept_full"]
CI_LO, CI_HI = tierB_ext_summary["slope_full_lo"], tierB_ext_summary["slope_full_hi"]

per_locus = tierB_ext_slopes[["locus", "K", "slope"]].sort_values("slope").reset_index(drop=True)

loc_theta = tierB_ext.drop_duplicates("locus")[["locus", "K", "theta1", "theta2"]].copy()
loc_theta["weaker_theta"] = loc_theta[["theta1", "theta2"]].min(axis=1)
loc_theta["noncentrality"] = loc_theta["K"] * loc_theta["weaker_theta"]
per_locus = per_locus.merge(loc_theta[["locus", "noncentrality"]], on="locus", how="left")

N_PROC = len(proc)
BONF = 0.05 / N_PROC


def genome_axis(df):
    """Cumulative genomic position for a locus-level Manhattan, plus chromosome"""
    d = df.sort_values(["chr", "start"]).copy()
    offsets, centres, edges, run = {}, [], [], 0.0
    for c in sorted(d.chr.unique()):
        g = d[d.chr == c]
        span = float(g.stop.max())
        offsets[c] = run
        centres.append(run + span / 2.0)
        edges.append((run, run + span))
        run += span
    d["gpos"] = d.chr.map(offsets) + (d.start + d.stop) / 2.0
    return d, centres, edges, run


gw, CHR_CENTRES, CHR_EDGES, GENOME_LEN = genome_axis(proc)
CHR_ORDER = sorted(proc.chr.unique())
CHR_LABELS = [str(c) if (c <= 18 or c % 2 == 0) else "" for c in CHR_ORDER]


def _chrom_bands(ax, ymin, ymax):
    """Alternating chromosome bands as recessive background."""
    for i, (a, b) in enumerate(CHR_EDGES):
        if i % 2 == 1:
            ax.axvspan(a, b, color=S.BAND, lw=0, zorder=0)
    ax.set_xlim(0, GENOME_LEN)
    ax.set_ylim(ymin, ymax)


def figure1_genomewide():
    """The gate, genome-wide. Three stacked panels on a shared genomic x-axis."""
    fig, axes = plt.subplots(3, 1, figsize=(S.PANEL_W * 2.9, S.PANEL_H * 1.5),
                             sharex=True)

    band = gw.chr.map(lambda c: S.CHROM_B if CHR_ORDER.index(c) % 2 == 0
                      else S.CHROM_A)

    for ax, trait, ymax in ((axes[0], "ASD", 5.4), (axes[1], "SCZ", 42.0)):
        p = gw[f"p_{trait}"]
        ok = p.notna()
        yv = -np.log10(p[ok])
        rug = ymax * 0.055

        _chrom_bands(ax, -rug * 1.8, ymax)
        ax.scatter(gw.gpos[ok], yv, s=S.S_DENSE, c=band[ok], linewidths=0,
                   zorder=3)
        ax.vlines(gw.gpos[~ok], -rug * 1.65, -rug * 0.4, color=S.C_GREY,
                  lw=S.LW_RUG, alpha=S.A_RUG, zorder=2)

        ax.axhline(-np.log10(0.05), color=S.C_ORANGE, lw=S.LW_REF, zorder=4)
        ax.axhline(-np.log10(BONF), color=S.C_RED, ls="--", lw=S.LW_REF,
                   zorder=4)
        S.refline(ax, 0, kind="zero", zorder=4)

        n_pass = int((p < 0.05).sum())
        ax.set_title(f"{trait}: {n_pass:,} of {N_PROC:,} blocks clear "
                     f"$P<0.05$ ({100*n_pass/N_PROC:.1f}%); "
                     f"{int((~ok).sum()):,} not estimable")
        ax.set_ylabel(r"$-\log_{10}P$  (local $h^2$)")
        ax.grid(axis="x", visible=False)
        ax.tick_params(axis="x", length=0)

    axes[0].set_yticks([0, 1, 2, 3, 4, 5])
    axes[1].set_yticks([0, 10, 20, 30, 40])

    ax = axes[2]
    r = gw[gw.rho_hat.notna()]
    _chrom_bands(ax, -1.12, 1.12)
    S.refline(ax, 0, kind="zero", zorder=4)

    for _, row in r.iterrows():
        c = S.C_POS if row.rho_hat > 0 else S.C_NEG
        ax.plot([row.gpos, row.gpos], [row.rho_lower, row.rho_upper], color=c,
                lw=S.LW_INTERVAL, alpha=S.A_INTERVAL, solid_capstyle="butt",
                zorder=3)
    det = (np.sign(r.rho_lower) == np.sign(r.rho_upper)).values
    cols = np.where(r.rho_hat > 0, S.C_POS, S.C_NEG)
    ax.scatter(r.gpos[~det], r.rho_hat[~det], s=S.S_HOLLOW, facecolors="none",
               edgecolors=cols[~det], linewidths=S.MEW, zorder=4)
    ax.scatter(r.gpos[det], r.rho_hat[det], s=S.S_FILLED, c=cols[det],
               alpha=S.A_SCATTER, linewidths=0, zorder=5)

    n_pos = int((r.rho_hat > 0).sum())
    ax.set_title(f"Reported local genetic correlation at {len(r)} loci: "
                 f"{n_pos} positive ({100*n_pos/len(r):.1f}%), "
                 f"{int(det.sum())} excluding zero")
    ax.set_ylabel(r"$\hat\rho$  (95% CI)")
    ax.set_xlabel("chromosome")
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xticks(CHR_CENTRES)
    ax.set_xticklabels(CHR_LABELS)
    ax.grid(axis="x", visible=False)
    ax.tick_params(axis="x", length=0)

    fig.legend(handles=[
        Line2D([], [], color=S.C_ORANGE, lw=S.LW_REF,
               label="$P$ = 0.05 (the filter)"),
        Line2D([], [], color=S.C_RED, ls="--", lw=S.LW_REF,
               label=f"$P$ = 0.05/{N_PROC:,} (Bonferroni)"),
        Line2D([], [], color=S.C_GREY, lw=S.LW_REF,
               label="phenotype not estimable"),
        Line2D([], [], color=S.C_POS, marker="o", lw=0, alpha=S.A_SCATTER,
               label=r"$\hat\rho$ excludes zero"),
        Line2D([], [], color="none", marker="o", lw=0, markeredgecolor=S.C_GREY,
               markeredgewidth=S.MEW, label=r"$\hat\rho$ crosses zero"),
    ], loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=5,
        fontsize=S.FS_LEGEND)

    S.save(fig, f"{OUT}/Figure1_genomewide")


def figure2_attenuation():
    fig, ax = plt.subplots(1, 3, figsize=(S.PANEL_W * 3, S.PANEL_H))

    a = ax[0]
    rho = np.linspace(0, 0.9, 200)
    for key, K in [("S5_smallK", 100), ("S1_asdlike_x_sczlike", 300),
                   ("S6_largeK", 500)]:
        d = sim[sim.scenario == key].sort_values("rho")
        a.plot(d.rho, d.mean_rho_hat, "o-", ms=S.MS_SERIES,
               label=f"Tier A, K={K}")
    a.fill_between(rho, CI_LO * rho + pooled_int, CI_HI * rho + pooled_int,
                   color=S.C_RED, alpha=S.A_BAND, lw=0)
    a.plot(rho, pooled_slope * rho + pooled_int, color=S.C_RED, lw=S.LW_EMPH,
           label=f"Tier B, real loci (slope {pooled_slope:.4f})")
    tb = tierB_ext.groupby("rho_true").mean_rho_hat.mean()
    a.plot(tb.index, tb.values, "o", color=S.C_RED, ms=S.MS_SERIES + 2)
    a.plot([0, 0.9], [0, 0.9], "k--", lw=S.LW_REF, label="unbiased")
    a.set_xlabel(r"true $\rho$")
    a.set_ylabel(r"E[$\hat\rho$ | gate]")
    a.set_xlim(0, 0.95)
    a.set_ylim(0, 0.95)
    a.set_title("Tier A understates the attenuation Tier B measures")
    a.legend(fontsize=S.FS_LEGEND, loc="upper left")

    a = ax[1]
    lo, hi = per_locus.iloc[0], per_locus.iloc[-1]
    a.axhspan(CI_LO, CI_HI, color=S.C_RED, alpha=S.A_BAND, lw=0)
    a.axhline(pooled_slope, color=S.C_RED, lw=S.LW_EMPH,
              label=f"pooled slope {pooled_slope:.4f}")
    a.scatter(per_locus.K, per_locus.slope, s=S.S_SCATTER * 1.6,
              color=S.C_BLUE, alpha=S.A_SCATTER, label="one locus")
    for row, dy, va in ((lo, -12, "top"), (hi, 12, "bottom")):
        a.annotate(f"{row.slope:.3f}", xy=(row.K, row.slope), xytext=(0, dy),
                   textcoords="offset points", ha="center", va=va, fontsize=9)
    a.set_xlabel("locus size K")
    a.set_ylabel("attenuation slope at that locus")
    a.set_ylim(0.28, 1.02)
    a.set_title(f"Per-locus slope spans {lo.slope:.3f}"
                f"-{hi.slope:.3f} ({hi.slope/lo.slope:.1f}-fold)")
    a.legend(fontsize=S.FS_LEGEND, loc="upper left")

    a = ax[2]
    mech_slope = tierB_ext_summary["mech_slope"]
    mech_p = tierB_ext_summary["mech_p"]
    x = np.log(per_locus.noncentrality)
    b, c0 = np.polyfit(x, per_locus.slope, 1)
    xx = np.linspace(x.min(), x.max(), 100)
    a.scatter(per_locus.noncentrality, per_locus.slope, s=S.S_SCATTER * 1.6,
              color=S.C_BLUE, alpha=S.A_SCATTER, label="one locus (extended Tier B)")
    a.plot(np.exp(xx), b * xx + c0, color=S.C_RED, lw=S.LW_EMPH,
           label=f"slope on log noncentrality {mech_slope:.3f} ($P$={mech_p:.1e})")
    a.set_xscale("log")
    a.set_xlabel(r"weaker trait's noncentrality $K\theta$")
    a.set_ylabel("attenuation slope at that locus")
    a.set_title("Channel 2 selects on power: survivors skew\ntoward less-attenuated, higher-noncentrality loci")
    a.legend(fontsize=S.FS_LEGEND, loc="upper left")

    S.save(fig, f"{OUT}/Figure2_attenuation")


def figure3_intervals():
    post = cov[cov.gate == "post_gate"]
    w = rep.width.values
    fig, ax = plt.subplots(1, 2, figsize=(S.PANEL_W * 2, S.PANEL_H))

    a = ax[0]
    a.hist(w, bins=np.arange(0.55, 2.05, 0.10), color=S.C_BLUE,
           alpha=S.A_SCATTER, edgecolor=S.C_BLUE, zorder=3)
    a.axvline(w.mean(), color=S.C_RED, lw=S.LW_EMPH,
              label=f"observed mean {w.mean():.3f}")
    a.axvline(post.mean_width.mean(), color=S.C_GREEN, lw=S.LW_EMPH,
              label=f"simulation predicted {post.mean_width.mean():.3f}")
    a.axvline(2.0, color="k", ls="--", lw=S.LW_REF,
              label=r"width 2 = full $[-1,1]$")
    a.axvline(1.0, color="k", lw=S.LW_ZERO, alpha=S.A_QUIET, zorder=4,
              label=f"{100*(w>1).mean():.1f}% wider than 1")
    a.set_xlabel("95% confidence-interval width")
    a.set_ylabel(f"loci (of {len(w)})")
    a.set_xlim(0.5, 2.1)
    a.set_ylim(0, 14.5)
    a.set_title(f"Intervals are calibrated ({post.coverage.mean():.3f}) "
                "and nearly vacuous")
    a.legend(fontsize=S.FS_LEGEND, loc="upper left")

    a = ax[1]
    g0 = signdet[signdet.scenario == "ASDlike_x_SCZlike"].sort_values("rho")
    g1 = signdet[signdet.scenario == "with_overlap_r10"].sort_values("rho")
    obs = float(rep.determinate.mean())

    a.plot(g0.rho, 1 - g0.spans_zero, "o-", ms=S.MS_SERIES, color=S.C_BLUE,
           label="$r_e$=0.0")
    a.plot(g1.rho, 1 - g1.spans_zero, "o-", ms=S.MS_SERIES, color=S.C_ORANGE,
           label="$r_e$=0.10")
    a.axhline(obs, color=S.C_RED, lw=S.LW_EMPH,
              label=f"observed in real data {obs:.3f}")
    for r_, dy, va in ((-0.4, 10, "bottom"), (0.4, -11, "top")):
        v = float(1 - g0[np.isclose(g0.rho, r_)].spans_zero.iloc[0])
        a.annotate(f"{v:.3f}", xy=(r_, v), xytext=(0, dy),
                   textcoords="offset points", ha="center", va=va, fontsize=9,
                   color=S.C_BLUE)
    a.set_xlabel(r"true $\rho$")
    a.set_ylabel("fraction of intervals excluding zero")
    a.set_xlim(-0.55, 0.7)
    a.set_ylim(0, 0.33)
    a.set_title(r"Determinacy is near-symmetric in $\rho$,"
                "\nso its count cannot discriminate")
    a.legend(fontsize=S.FS_LEGEND, loc="upper left")

    S.save(fig, f"{OUT}/Figure3_intervals")


def figure4_the_78():
    d = rep.sort_values("rho_hat").reset_index(drop=True)
    y = np.arange(len(d))
    det = d.determinate.values
    n_pos = int((d.rho_hat > 0).sum())
    n_neg = int((d.rho_hat < 0).sum())
    n_det = int(det.sum())
    mean_r = d.rho_hat.mean()

    fig, ax = plt.subplots(figsize=(S.PANEL_W * 1.2, S.PANEL_H * 2.0))

    for yi, row in d.iterrows():
        c = S.C_POS if row.rho_hat > 0 else S.C_NEG
        ax.plot([row.rho_lower, row.rho_upper], [yi, yi], color=c,
                lw=S.LW_INTERVAL, alpha=S.A_INTERVAL, solid_capstyle="butt",
                zorder=3)
    cols = np.where(d.rho_hat > 0, S.C_POS, S.C_NEG)
    ax.scatter(d.rho_hat[~det], y[~det], s=S.S_HOLLOW, facecolors="none",
               edgecolors=cols[~det], linewidths=S.MEW, zorder=4,
               label=f"crosses zero ({len(d)-n_det})")
    ax.scatter(d.rho_hat[det], y[det], s=S.S_FILLED, c=cols[det],
               alpha=S.A_SCATTER, linewidths=0, zorder=5,
               label=f"excludes zero ({n_det})")

    S.refline(ax, 0, axis="x", kind="zero", zorder=4)
    ax.axvline(mean_r, color="k", ls="--", lw=S.LW_REF, zorder=3,
               label=f"mean {mean_r:.3f}")

    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.5, len(d) + 0.5)
    ax.set_yticks([])
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel(r"$\hat\rho$  (95% CI)")
    ax.set_ylabel(f"the {len(d)} reported loci, sorted")
    ax.grid(axis="y", visible=False)
    ax.set_title(f"{n_pos} positive ({100*n_pos/len(d):.1f}%), {n_neg} negative;"
                 f"\n{len(d)-n_det} of {len(d)} intervals cross zero")
    ax.legend(fontsize=S.FS_LEGEND, loc="lower right")

    S.save(fig, f"{OUT}/Figure4_the_78_loci")


if __name__ == "__main__":
    print("building figures ->", OUT)
    figure1_genomewide()
    figure2_attenuation()
    figure3_intervals()
    figure4_the_78()
    print("done.")
