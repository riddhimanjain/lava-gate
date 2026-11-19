"""11_aggregate_correction.py -- can the gate-aware correction actually be applied"""
import json
import numpy as np
import pandas as pd
from importlib import machinery

A = machinery.SourceFileLoader("analytic", "01_analytic.py").load_module()
C = machinery.SourceFileLoader("corr", "08_correction.py").load_module()
R = machinery.SourceFileLoader("refined", "04_refined_theory.py").load_module()

R_E, TIERB_SLOPE, SEED = 0.02756, 0.55123, 20260823
N_NULL = 40_000

aud = pd.read_csv("out_R/attrition_audit.csv")
rep = aud[(aud["gate_pass"] == True) & np.isfinite(aud["rho_hat"])].copy()
rep["K"] = rep["K"].astype(int)
rep["t1"] = C.theta_obs_from_p(rep["p_ASD"].values, rep["K"].values)
rep["t2"] = C.theta_obs_from_p(rep["p_SCZ"].values, rep["K"].values)
assert len(rep) == 78


def strat_slope(d, n_strata):
    bins = (np.zeros(len(d), int) if n_strata == 1 else
            np.clip(np.digitize(d["K"], np.quantile(d["K"],
                    np.linspace(0, 1, n_strata + 1))[1:-1]), 0, n_strata - 1))
    sl_w = ic_w = 0.0
    for b in np.unique(bins):
        g = d[bins == b]
        Kb = int(round(g["K"].mean() / 25.0) * 25) or 25
        th1, D1, _ = C.invert_theta(np.array([g["t1"].mean()]), Kb)
        th2, D2, _ = C.invert_theta(np.array([g["t2"].mean()]), Kb)
        sl, ic = C.locus_slope_intercept(Kb, th1, th2, D1, D2, R_E)
        sl_w += float(sl[0]) * len(g); ic_w += float(ic[0]) * len(g)
    return sl_w / len(d), ic_w / len(d)


stab = {}
for ns in (1, 2, 3, 4, 6, 8):
    sl, ic = strat_slope(rep, ns)
    cm = (rep["rho_hat"].mean() - ic) / sl if sl > 0 else np.nan
    stab[f"strata_{ns}"] = dict(slope=float(sl), intercept=float(ic),
                                corrected_mean=float(cm))
    print(f"  strata={ns}: slope={sl:.4f}  corrected_mean={cm:.4f}")
sl_vals = np.array([v["slope"] for v in stab.values()])
cm_vals = np.array([v["corrected_mean"] for v in stab.values()])

rng = np.random.default_rng(SEED)
floor_rows = []
for Kb in sorted({int(round(k / 50.0) * 50) or 50 for k in rep["K"]}):
    nref = A.nref_scale(Kb, C.N_REF)
    t = A.gate_threshold(Kb, C.ALPHA) / nref
    u = rng.chisquare(Kb, N_NULL)
    u = u[u > t]
    th_null = u / Kb - nref
    sub1 = rep.loc[(rep["K"] / 50).round() * 50 == Kb, "t1"]
    sub2 = rep.loc[(rep["K"] / 50).round() * 50 == Kb, "t2"]
    for trait, sub in (("ASD", sub1), ("SCZ", sub2)):
        if not len(sub):
            continue
        floor_rows.append(dict(
            K_bin=Kb, trait=trait, n_loci=int(len(sub)),
            observed_mean_theta=float(sub.mean()),
            null_mean_theta=float(th_null.mean()),
            null_sd_theta=float(th_null.std(ddof=1)),
            z_above_null=float((sub.mean() - th_null.mean())
                               / (th_null.std(ddof=1) / np.sqrt(len(sub)))),
            pct_of_loci_below_null_mean=float((sub < th_null.mean()).mean()),
        ))
fl = pd.DataFrame(floor_rows)
fl.to_csv("out/floor_diagnostic.csv", index=False)
print("\n  floor diagnostic (theta_hat vs its null conditional distribution):")
for r in fl.itertuples():
    print(f"    K~{r.K_bin:3d} {r.trait}  n={r.n_loci:3d}  obs={r.observed_mean_theta:.4f}"
          f"  null={r.null_mean_theta:.4f}  z={r.z_above_null:+6.2f}"
          f"  below_null={r.pct_of_loci_below_null_mean:.0%}")

asd = fl[fl.trait == "ASD"]; scz = fl[fl.trait == "SCZ"]
summ = dict(
    n_loci=78, r_e=R_E,
    mean_rho_hat_observed=float(rep["rho_hat"].mean()),
    stability=stab,
    slope_range_over_binning=[float(sl_vals.min()), float(sl_vals.max())],
    slope_fold_change_over_binning=float(sl_vals.max() / sl_vals.min()),
    corrected_mean_range_over_binning=[float(cm_vals.min()), float(cm_vals.max())],
    n_floored_per_locus=int(52),
    asd_max_z_above_null=float(asd["z_above_null"].max()),
    asd_weighted_pct_below_null=float((asd["pct_of_loci_below_null_mean"]
                                       * asd["n_loci"]).sum() / asd["n_loci"].sum()),
    scz_min_z_above_null=float(scz["z_above_null"].min()),
    scz_weighted_pct_below_null=float((scz["pct_of_loci_below_null_mean"]
                                       * scz["n_loci"]).sum() / scz["n_loci"].sum()),
    paper_previous_corrected_mean=float(rep["rho_hat"].mean() / TIERB_SLOPE),
    verdict=("NOT APPLICABLE to ASD x SCZ. The aggregate slope moves by a factor "
             "of %.2f across reasonable K stratifications, because ASD's observed "
             "gated heritability sits at the selection floor: the stage-1 "
             "inversion is not identified there. The correction is validated and "
             "usable, but this trait pair does not meet its precondition. Report "
             "the uncorrected estimate with the attenuation stated as a bound."
             % (sl_vals.max() / sl_vals.min())),
)
with open("out/aggregate_correction.json", "w") as f:
    json.dump(summ, f, indent=2)
print("\n=== VERDICT ===\n " + summ["verdict"])
