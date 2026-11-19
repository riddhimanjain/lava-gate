"""10_apply_correction.py -- the gate-aware correction applied to ASD x SCZ."""
import json
import numpy as np
import pandas as pd
from importlib import machinery

C = machinery.SourceFileLoader("corr", "08_correction.py").load_module()

R_E = 0.02756
POOLED_SLOPE, POOLED_INTERCEPT = 0.55123, 0.028293

aud = pd.read_csv("out_R/attrition_audit.csv")
rep = aud[(aud["gate_pass"] == True) & np.isfinite(aud["rho_hat"])].copy()
assert len(rep) == 78, f"expected 78 reported loci, got {len(rep)}"

rows = []
for r in rep.itertuples():
    K = int(r.K)
    t1o = C.theta_obs_from_p(r.p_ASD, K)
    t2o = C.theta_obs_from_p(r.p_SCZ, K)
    out = C.correct_from_theta_obs(np.array([r.rho_hat]), K,
                                   np.array([t1o]), np.array([t2o]), R_E)
    sl = float(out["slope"][0])
    inv = (1.0 / sl) if sl > 0 else np.nan
    rows.append(dict(
        locus=r.locus, chr=r.chr, K=K,
        p_ASD=r.p_ASD, p_SCZ=r.p_SCZ,
        theta_ASD=float(out["theta1"][0]), theta_SCZ=float(out["theta2"][0]),
        slope=sl, intercept=float(out["intercept"][0]),
        floored=bool(sl <= 0),
        rho_hat=r.rho_hat,
        rho_corrected=float(out["rho_corrected"][0]),
        rho_corrected_unclipped=float(out["rho_corrected_raw"][0]),
        rho_lower=r.rho_lower, rho_upper=r.rho_upper,
        ci_width=r.rho_upper - r.rho_lower,
        ci_width_corrected=(r.rho_upper - r.rho_lower) * inv,
    ))
df = pd.DataFrame(rows)
df.to_csv("out/applied_correction.csv", index=False)

th1b, D1b, _ = C.invert_theta(np.array([C.theta_obs_from_p(rep["p_ASD"], rep["K"]).mean()]),
                              int(round(rep["K"].mean())))
th2b, D2b, _ = C.invert_theta(np.array([C.theta_obs_from_p(rep["p_SCZ"], rep["K"]).mean()]),
                              int(round(rep["K"].mean())))
Kbar = int(round(rep["K"].mean()))
sl_b, ic_b = C.locus_slope_intercept(Kbar, th1b, th2b, D1b, D2b, R_E)
agg = float((rep["rho_hat"].mean() - ic_b[0]) / sl_b[0])

summ = dict(
    n_loci=int(len(df)), r_e=R_E, K_mean=Kbar,
    mean_rho_hat=float(rep["rho_hat"].mean()),
    n_floored_no_recoverable_signal=int(df["floored"].sum()),
    per_locus_slope_min=float(df.loc[~df["floored"], "slope"].min()),
    per_locus_slope_max=float(df.loc[~df["floored"], "slope"].max()),
    per_locus_slope_median=float(df.loc[~df["floored"], "slope"].median()),
    per_locus_slope_iqr=[float(df.loc[~df["floored"], "slope"].quantile(.25)),
                         float(df.loc[~df["floored"], "slope"].quantile(.75))],
    per_locus_intercept_min=float(df["intercept"].min()),
    per_locus_intercept_max=float(df["intercept"].max()),
    per_locus_intercept_median=float(df["intercept"].median()),
    aggregate_slope=float(sl_b[0]), aggregate_intercept=float(ic_b[0]),
    mean_rho_corrected_aggregate=agg,
    mean_rho_corrected_pooledfactor=float(((rep["rho_hat"] - POOLED_INTERCEPT)
                                           / POOLED_SLOPE).mean()),
    mean_rho_corrected_paper_naive=float(rep["rho_hat"].mean() / POOLED_SLOPE),
    n_corrected_outside_unit=int((df["rho_corrected_unclipped"].abs() > 1).sum()),
    mean_ci_width=float(df["ci_width"].mean()),
    mean_ci_width_corrected=float(df["ci_width_corrected"].mean(skipna=True)),
    frac_corrected_width_over_2=float((df["ci_width_corrected"].dropna() > 2).mean()),
)
with open("out/applied_correction_summary.json", "w") as f:
    json.dump(summ, f, indent=2)
for k, v in summ.items():
    print(f"  {k:36s} {v}")
