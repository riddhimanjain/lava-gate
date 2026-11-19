"""16_audit_numbers.py -- every numeric claim in the manuscript, checked against"""
import argparse
import json
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MS = os.path.join(ROOT, "manuscript", "manuscript.md")
RES = os.path.join(ROOT, "results")

_cache = {}


def J(rel):
    p = os.path.join(RES, rel)
    if p not in _cache:
        _cache[p] = json.load(open(p)) if os.path.exists(p) else None
    return _cache[p]


def D(rel):
    p = os.path.join(RES, rel)
    if p not in _cache:
        _cache[p] = pd.read_csv(p) if os.path.exists(p) else None
    return _cache[p]


def txt(rel):
    p = os.path.join(RES, rel)
    if p not in _cache:
        _cache[p] = open(p).read() if os.path.exists(p) else None
    return _cache[p]


def num_from(pattern, source, group=1):
    """Pull a number out of a text report (tierB.txt and friends)."""
    if source is None:
        return None
    m = re.search(pattern, source)
    return float(m.group(group)) if m else None


def claims():
    cf = D("correction/closedform_vs_tierA.csv")
    ho = J("theory_tierA/heldout_summary.json")
    cov = J("theory_tierA/coverage_summary.json")
    cs = J("correction/correction_summary.json")
    ac = J("correction/applied_correction_summary.json")
    ag = J("correction/aggregate_correction.json")
    sb = J("correction/signal_bound.json")
    ns = J("correction/null_sign_summary.json")
    tb = txt("applied_tierB/tierB.txt")
    at = txt("applied_tierB/attrition_audit.txt")

    C = []

    def add(cid, sec, printed, tol, fn):
        C.append((cid, sec, printed, tol, fn))

    if cf is not None:
        for _, r in cf.iterrows():
            add(f"T3_{r.scenario}_simslope", "3.1", f"{r.sim_slope:.4f}", 5e-5,
                lambda r=r: r.sim_slope)
            add(f"T3_{r.scenario}_cfslope", "3.1", f"{r.cf_slope:.4f}", 5e-5,
                lambda r=r: r.cf_slope)
        add("3.1_mean_reldiff", "3.1", "1.14%", 5e-5,
            lambda: cf.reldiff.mean() * 100 / 100)
        add("3.1_max_reldiff", "3.1", "2.06%", 5e-5,
            lambda: cf.reldiff.max() * 100 / 100)
        add("3.1_max_abs_intercept_err", "3.1", "0.0052", 5e-5,
            lambda: (cf.sim_intercept - cf.cf_intercept).abs().max())
        add("3.1_slope_K100", "3.1", "0.4555", 5e-5,
            lambda: float(cf.loc[cf.K == 100, "cf_slope"].iloc[0]))
        add("3.1_slope_K500", "3.1", "0.6897", 5e-5,
            lambda: float(cf.loc[cf.K == 500, "cf_slope"].iloc[0]))
        add("Abs_cfslope_min", "Abs", "0.455", 5e-4, lambda: cf.cf_slope.min())
        add("Abs_cfslope_max", "Abs", "0.690", 5e-4, lambda: cf.cf_slope.max())
    if ho:
        add("3.1_heldout_cells", "3.1", "42", 0, lambda: ho["n_cells"])
        add("3.1_heldout_mae", "3.1", "0.0066", 5e-5, lambda: ho["v2_mean_abs_err"])
        add("3.1_heldout_within01", "3.1", "88%", 0.005,
            lambda: ho["v2_frac_within_0.01"])

    xv = J("crossvalidation/r_vs_python.json")
    if xv:
        add("2.4_xval_ncombos", "2.4", "500", 0,
            lambda: xv["n_param_combinations"])
        add("2.4_xval_K_min", "2.4", "50", 0, lambda: xv["K_min"])
        add("2.4_xval_K_max", "2.4", "500", 0, lambda: xv["K_max"])
        add("2.4_xval_theta_min", "2.4", "0.005", 1e-9, lambda: xv["theta_min"])
        add("2.4_xval_theta_max", "2.4", "0.30", 1e-9, lambda: xv["theta_max"])
        add("2.4_xval_re_max", "2.4", "0.20", 1e-9, lambda: xv["r_e_max"])
        C.append(("2.4_xval_delta_bound", "2.4", "2 × 10⁻¹²", "BOUND",
                  lambda: xv["max_diff_delta"]))
        C.append(("2.4_xval_slope_bound", "2.4", "2 × 10⁻¹⁵", "BOUND",
                  lambda: xv["max_diff_slope_intercept"]))

    if tb:
        add("3.2_pooled_slope", "3.2", "0.5512", 5e-5,
            lambda: num_from(r"Pooled slope\s*:\s*([0-9.]+)", tb))
        add("3.2_pooled_intercept", "3.2", "0.0283", 5e-5,
            lambda: num_from(r"Pooled intercept\s*:\s*([0-9.]+)", tb))

    tx = J("tierB_extended/tierB_extended_summary.json")
    if tx:
        add("3.2x_candidates", "3.2", "119", 0, lambda: tx["n_selected"])
        add("3.2x_usable", "3.2", "46", 0, lambda: tx["n_usable"])
        add("3.2x_ch2_loss", "3.2", "61.3%", 0.0005,
            lambda: tx["survivorship_loss_pct"] / 100)
        add("3.2x_slope", "3.2", "0.6236", 5e-5, lambda: tx["slope_full"])
        add("3.2x_slope_lo", "3.2", "0.6102", 5e-5, lambda: tx["slope_full_lo"])
        add("3.2x_slope_hi", "3.2", "0.6371", 5e-5, lambda: tx["slope_full_hi"])
        add("3.2x_intercept", "3.2", "0.0182", 5e-5, lambda: tx["intercept_full"])
        add("3.2x_r2", "3.2", "0.9527", 5e-5, lambda: tx["r2_full"])
        add("3.2x_quad_p", "3.2", "0.936", 5e-4, lambda: tx["quad_term_p"])
        add("3.2x_mech_slope", "3.2", "+0.178", 5e-4, lambda: tx["mech_slope"])
        add("3.2x_perlocus_min", "3.2", "0.371", 5e-4,
            lambda: tx["per_locus_slope_min"])
        add("3.2x_perlocus_max", "3.2", "0.971", 5e-4,
            lambda: tx["per_locus_slope_max"])
        add("3.2x_origgrid_slope", "3.2", "0.6202", 5e-5,
            lambda: tx["slope_origgrid"])
        add("3.2x_origgrid_lo", "3.2", "0.5912", 5e-5,
            lambda: tx["slope_origgrid_lo"])
        add("3.2x_origgrid_hi", "3.2", "0.6493", 5e-5,
            lambda: tx["slope_origgrid_hi"])
        add("3.2x_attenuation_pct", "3.2", "38%", 0.005,
            lambda: 1 - tx["slope_full"])

    if cov:
        add("3.4_cov_pre", "3.4", "0.974", 5e-4,
            lambda: cov["pre_gate"]["mean_coverage"])
        add("3.4_cov_post", "3.4", "0.949", 5e-4,
            lambda: cov["post_gate"]["mean_coverage"])
        add("3.4_width_pre", "3.4", "1.651", 5e-4,
            lambda: cov["pre_gate"]["mean_width"])
        add("3.4_width_post", "3.4", "1.386", 5e-4,
            lambda: cov["post_gate"]["mean_width"])
        add("3.4_fracgt1_pre", "3.4", "0.936", 5e-4,
            lambda: cov["pre_gate"]["mean_frac_width_gt_1"])
        add("3.4_fracgt1_post", "3.4", "0.869", 5e-4,
            lambda: cov["post_gate"]["mean_frac_width_gt_1"])
    if ac:
        add("3.4_width_real", "3.4", "1.365", 5e-4, lambda: ac["mean_ci_width"])

    if cs:
        add("T7_bias_raw", "3.5", "0.247", 5e-4, lambda: cs["mean_abs_bias_raw"])
        add("T7_bias_pooled", "3.5", "0.195", 5e-4, lambda: cs["mean_abs_bias_pooled"])
        add("T7_bias_perlocus", "3.5", "0.106", 5e-4,
            lambda: cs["mean_abs_bias_corrected"])
        add("T7_bias_aggregate", "3.5", "0.026", 5e-4,
            lambda: cs["mean_abs_bias_corrected_aggregate"])
        add("T7_max_raw", "3.5", "0.602", 5e-4, lambda: cs["max_abs_bias_raw"])
        add("T7_max_pooled", "3.5", "0.579", 5e-4, lambda: cs["max_abs_bias_pooled"])
        add("T7_max_perlocus", "3.5", "0.304", 5e-4,
            lambda: cs["max_abs_bias_corrected"])
        add("T7_max_aggregate", "3.5", "0.107", 5e-4,
            lambda: cs["max_abs_bias_corrected_aggregate"])
        add("3.5_pooled_slope_used", "3.5", "0.6236", 5e-5,
            lambda: cs["pooled_slope_used"])
        add("3.5_infl_min", "3.5", "1.95", 5e-3, lambda: cs["min_width_inflation"])
        add("3.5_infl_max", "3.5", "3.01", 5e-3, lambda: cs["max_width_inflation"])
        add("3.5_infl_mean", "3.5", "2.58", 5e-3, lambda: cs["mean_width_inflation"])
    if ac:
        add("3.5_corrected_width", "3.5", "2.27", 5e-3,
            lambda: ac["mean_ci_width_corrected"])
        add("3.5_frac_width_gt2", "3.5", "42%", 0.005,
            lambda: ac["frac_corrected_width_over_2"])

    if at:
        add("3.6_processed", "3.6", "2,400", 0,
            lambda: num_from(r"processed by process\.locus:\s*(\d+)", at))
    ad = D("applied_tierB/attrition_audit.csv")
    if ad is not None:
        pro = ad[ad["processed"] == True]
        both = pro[pro["p_ASD"].notna() & pro["p_SCZ"].notna()]
        add("3.6_both_estimable", "3.6", "906", 0, lambda: len(both))
        add("3.6_reported", "3.6", "78", 0, lambda: int(ad["rho_hat"].notna().sum()))
        add("3.6_ch2_asd_n", "3.6", "1,429", 0,
            lambda: int((pro["p_ASD"].isna()).sum()))
        add("3.6_ch2_scz_n", "3.6", "65", 0,
            lambda: int((pro["p_SCZ"].isna()).sum()))
        add("3.6_ch2_asd_pct", "3.6", "59.5%", 0.0005,
            lambda: (pro["p_ASD"].isna()).mean())
        add("3.6_ch2_scz_pct", "3.6", "2.7%", 0.0005,
            lambda: (pro["p_SCZ"].isna()).mean())
        add("3.6_cond_pass", "3.6", "10.26%", 0.0005,
            lambda: (both["p_ASD"] < 0.05).mean())
        add("3.6_uncond_pass", "3.6", "4.25%", 0.0005,
            lambda: (pro["p_ASD"] < 0.05).mean())

    if ag:
        add("3.7_slope_fold", "3.7", "2.66", 5e-3,
            lambda: ag["slope_fold_change_over_binning"])
        add("3.7_slope_hi", "3.7", "0.769", 5e-4,
            lambda: max(ag["slope_range_over_binning"]))
        add("3.7_slope_lo", "3.7", "0.288", 5e-4,
            lambda: min(ag["slope_range_over_binning"]))
    fd = D("correction/floor_diagnostic.csv")
    if fd is not None:
        a = fd[fd.trait == "ASD"]; s = fd[fd.trait == "SCZ"]
        add("3.7_asd_zmin", "3.7", "−1.67", 5e-3, lambda: a.z_above_null.min())
        add("3.7_asd_zmax", "3.7", "+3.05", 5e-3, lambda: a.z_above_null.max())
        add("3.7_scz_zmin", "3.7", "+1.95", 5e-3, lambda: s.z_above_null.min())
        add("3.7_scz_zmax", "3.7", "+48.3", 5e-2, lambda: s.z_above_null.max())
    if sb:
        add("3.7_asd_enrich", "3.7", "2.05", 5e-3,
            lambda: sb["ASD"]["enrichment_over_alpha"])
        add("3.7_asd_fracnull", "3.7", "49%", 0.005,
            lambda: sb["ASD"]["fraction_of_passes_expected_null"])
        add("3.7_joint_expected", "3.7", "31.8", 0.05,
            lambda: sb["joint"]["expected_joint_if_ASD_all_null"])

    if ns:
        o, A_, B_ = ns["observed"], ns["partA"]["rho_0.0"], ns["partB"]
        add("3.8_n_positive", "3.8", "57", 0, lambda: o["n_positive"])
        add("3.8_n_determinate", "3.8", "10", 0, lambda: o["n_determinate"])
        add("3.8_null_ppos", "3.8", "0.531", 5e-4,
            lambda: A_["mean_p_positive"])
        add("3.8_expected_pos", "3.8", "41.4", 0.05,
            lambda: A_["expected_n_positive"])
        add("3.8_null_pdet", "3.8", "0.058", 5e-4,
            lambda: B_["null_mean_p_determinate"])
        add("3.8_expected_det", "3.8", "4.5", 0.05,
            lambda: B_["null_expected_n_determinate"])
        add("3.8_null_pos_given_det", "3.8", "0.601", 5e-4,
            lambda: B_["null_p_positive_given_determinate"])

    PANEL = {
        "P4":  ("panel/PANEL_SUMMARY_P4.csv", "ASD", "SCZ"),
        "P1":  ("panel/PANEL_SUMMARY_P1.csv", "SCZ", "BIP"),
        "P2a": ("panel/PANEL_SUMMARY_P2.csv", "MDD", "BIP"),
        "P2b": ("panel/PANEL_SUMMARY_P2.csv", "MDD_noUKB", "BIP_noUKB"),
        "P3":  ("panel/PANEL_SUMMARY_P3.csv", "ASD", "ADHD"),
        "P5":  ("panel/PANEL_SUMMARY_P5.csv", "ASD", "AN"),
    }
    T8 = {
        "P4":  ("2,400", "971", "906", "78", "0", "59.5%", "2.7%", "0.0276"),
        "P1":  ("2,487", "2,369", "2,305", "1,369", "98", "4.7%", "2.6%", "0.2139"),
        "P2a": ("2,483", "2,153", "2,082", "816", "20", "13.3%", "3.1%", "0.1090"),
        "P2b": ("2,483", "1,848", "1,810", "545", "8", "25.6%", "2.0%", "0.1132"),
        "P3":  ("2,424", "1,017", "977", "84", "0", "58.0%", "1.7%", "0.2372"),
        "P5":  ("2,458", "984", "961", "73", "0", "60.0%", "0.9%", "0.0945"),
    }
    AUDITS = {
        "P4": "applied_tierB/attrition_audit.csv",
        "P1": "panel/P1/attrition_audit.csv", "P2a": "panel/P2a/attrition_audit.csv",
        "P2b": "panel/P2b/attrition_audit.csv", "P3": "panel/P3/attrition_audit.csv",
        "P5": "panel/P5/attrition_audit.csv",
    }
    for pair, (csvp, weak, strong) in PANEL.items():
        _s = D(csvp)
        row = None if _s is None else _s[_s["pair"] == pair]
        pr, we, be, rp, bf, cw, cs_, re_ = T8[pair]
        if row is not None and len(row):
            r0 = row.iloc[0]
            add(f"T8_{pair}_processed", "3.6", pr, 0, lambda r0=r0: r0["processed"])
            add(f"T8_{pair}_both", "3.6", be, 0,
                lambda r0=r0: r0["both_estimable"])
            add(f"T8_{pair}_reported", "3.6", rp, 0, lambda r0=r0: r0["reported"])
            add(f"T8_{pair}_bonf", "3.6", bf, 0,
                lambda r0=r0: r0["bonferroni_survivors"])
            add(f"T8_{pair}_weakest", "3.6", we, 0,
                lambda r0=r0, w=weak: r0["est_1"]
                if str(r0["weaker_trait"]) == w and r0["est_1"] <= r0["est_2"]
                else min(r0["est_1"], r0["est_2"]))
        _ad = D(AUDITS[pair])
        if _ad is not None and f"p_{weak}" in _ad.columns:
            _pro = _ad[_ad["processed"] == True]
            add(f"T8_{pair}_ch2_weak", "3.6", cw, 0.0005,
                lambda pro=_pro, w=weak: pro[f"p_{w}"].isna().mean())
            add(f"T8_{pair}_ch2_strong", "3.6", cs_, 0.0005,
                lambda pro=_pro, s=strong: pro[f"p_{s}"].isna().mean())

    DIAG = {
        "P4": ("2.66", "−0.16", "0.49"), "P1": ("1.02", "+38.8", "0.07"),
        "P2a": ("1.03", "+18.0", "0.10"), "P2b": ("1.08", "+8.1", "0.14"),
        "P3": ("2.87", "−0.2", "0.44"), "P5": ("4.70", "−0.2", "0.53"),
    }
    for pair, (c1, z, fn) in DIAG.items():
        ps = J(f"panel/{pair}/pair_summary.json")
        if not ps:
            continue
        _dg = ps["diagnostic"]
        zz = _dg["C2_floor"]["pooled_z"]
        wk = min(zz, key=zz.get)
        add(f"D_{pair}_C1", "3.6", c1, 5e-3,
            lambda d=_dg: d["C1_stability"]["slope_fold_change"])
        add(f"D_{pair}_C2z", "3.6", z, 0.05, lambda zz=zz, wk=wk: zz[wk])
        add(f"D_{pair}_C3", "3.6", fn, 5e-3,
            lambda d=_dg, wk=wk: d["C3_signal"]["frac_null"][wk])

    v1 = D("verification/V1_summary.csv")
    if v1 is not None:
        ARMS = {
            "P4": ("1.228", "0.689", "0.647", "0.626"),
            "P1": ("1.163", "0.951", "0.939", "0.907"),
            "P2a": ("1.163", "0.965", "0.957", "0.923"),
            "P2b": ("1.203", "0.935", "0.915", "0.881"),
            "P3": ("1.210", "0.686", "0.635", "0.615"),
            "P5": ("1.217", "0.637", "0.601", "0.581"),
        }
        for pair, (_a, _b, _c, _d) in ARMS.items():
            r = v1[v1["pair"] == pair]
            if not len(r):
                continue
            r0 = r.iloc[0]
            for lbl, printed, col in (
                    ("A", _a, "slope_A_estimable_raw"),
                    ("B", _b, "slope_B_estimable_trunc"),
                    ("C", _c, "slope_C_gated_raw"),
                    ("D", _d, "slope_D_gated_trunc")):
                add(f"T6b_{pair}_{lbl}", "3.2", printed, 5e-4,
                    lambda r0=r0, col=col: r0[col])
        C.append(("3.2_trunc_min", "3.2", "0.020", "BOUND",
                  lambda: 0.020 - v1["trunc_effect_at_gate"].abs().min()))
        C.append(("3.2_trunc_max", "3.2", "0.034", "BOUND",
                  lambda: v1["trunc_effect_at_gate"].abs().max() - 0.034 + 1e-9))
        _tot = 1.0 - v1["slope_D_gated_trunc"]
        _shr = (v1["slope_C_gated_raw"] - v1["slope_D_gated_trunc"]).abs() / _tot
        _sh = dict(zip(v1["pair"], _shr * 100))
        _tt = dict(zip(v1["pair"], _tot * 100))
        for _p, _pc in (("P5", "4.7%"), ("P3", "5.3%"), ("P4", "5.6%"),
                        ("P2b", "28%"), ("P1", "34%"), ("P2a", "44%")):
            add(f"3.2_truncshare_{_p}", "3.2", _pc, 0.005,
                lambda _p=_p: _sh[_p] / 100.0)
        for _p, _pc in (("P5", "41.9%"), ("P3", "38.5%"), ("P4", "37.4%"),
                        ("P2a", "7.7%"), ("P1", "9.3%"), ("P2b", "11.9%")):
            add(f"3.2_totalatten_{_p}", "3.2", _pc, 0.0005,
                lambda _p=_p: _tt[_p] / 100.0)

    v4 = D("verification/V4_closedform_vs_tierB.csv")
    if v4 is not None:
        CF = {
            "P4": ("0.620", "0.624", "0.0174", "0.0182"),
            "P1": ("0.902", "0.908", "0.0360", "0.0283"),
            "P2a": ("0.917", "0.923", "0.0166", "0.0125"),
            "P2b": ("0.877", "0.882", "0.0249", "0.0205"),
            "P3": ("0.619", "0.612", "0.1476", "0.1389"),
            "P5": ("0.578", "0.579", "0.0696", "0.0652"),
        }
        for pair, (cs2, ms, ci, mi) in CF.items():
            r = v4[v4["pair"] == pair]
            if not len(r):
                continue
            r0 = r.iloc[0]
            add(f"T8_{pair}_cfslope", "3.6", cs2, 5e-4,
                lambda r0=r0: r0["cf_slope_mean"])
            add(f"T8_{pair}_tbslope", "3.6", ms, 5e-4,
                lambda r0=r0: r0["tierB_slope_D"])
            add(f"T4b_{pair}_cfint", "3.3", ci, 5e-5,
                lambda r0=r0: r0["cf_intercept_mean"])
            add(f"T4b_{pair}_tbint", "3.3", mi, 5e-5,
                lambda r0=r0: r0["tierB_intercept_D"])
            add(f"T4b_{pair}_rule", "3.3",
                f"{0.68 * float(r0['r_e']):.4f}", 5e-5,
                lambda r0=r0: r0["rule_0p68_intercept"])
        add("3.3_mae_rule", "3.3", "0.0432", 5e-5,
            lambda: (v4.rule_0p68_intercept - v4.tierB_intercept_D).abs().mean())
        add("3.3_mae_cf", "3.3", "0.005", 5e-4,
            lambda: (v4.cf_intercept_mean - v4.tierB_intercept_D).abs().mean())
        add("3.6_cf_worst", "3.6", "0.007", 5e-4,
            lambda: (v4.cf_slope_mean - v4.tierB_slope_D).abs().max())
        add("3.6_cf_mean", "3.6", "0.0048", 5e-5,
            lambda: (v4.cf_slope_mean - v4.tierB_slope_D).abs().mean())
        add("Abs_panel_slope_lo", "Abs", "0.579", 5e-4,
            lambda: v4.tierB_slope_D.min())
        add("Abs_panel_slope_hi", "Abs", "0.923", 5e-4,
            lambda: v4.tierB_slope_D.max())
        add("KP_panel_slope_lo", "KP", "0.58", 5e-3,
            lambda: v4.tierB_slope_D.min())
        add("KP_panel_slope_hi", "KP", "0.92", 5e-3,
            lambda: v4.tierB_slope_D.max())
        add("3.2_armC_lo", "3.2", "0.601", 5e-4,
            lambda: D("verification/V1_summary.csv").slope_C_gated_raw.min())
        add("3.2_armC_hi", "3.2", "0.957", 5e-4,
            lambda: D("verification/V1_summary.csv").slope_C_gated_raw.max())

    v3 = D("verification/V3_null_calibration.csv")
    if v3 is not None:
        asd = v3[(v3.pair == "P4") & (v3.trait == "ASD")]
        if len(asd):
            a0 = asd.iloc[0]
            add("3.7_ks_p", "3.7", "0.47", 5e-3, lambda a0=a0: a0["ks_p"])
            add("3.7_mean_u", "3.7", "0.491", 5e-4, lambda a0=a0: a0["mean_u"])
            add("3.7_n_estimable", "3.7", "971", 0,
                lambda a0=a0: a0["n_estimable"])
        for pr, printed in (("P3", "0.078"), ("P5", "0.86")):
            row = v3[(v3.pair == pr) & (v3.trait == "ASD")]
            if len(row):
                add(f"3.7_ks_{pr}", "3.7", printed, 5e-3,
                    lambda row=row: row.iloc[0]["ks_p"])
        for pr, printed in (("P4", "40.5%"), ("P3", "42.0%"), ("P5", "40.0%")):
            row = v3[(v3.pair == pr) & (v3.trait == "ASD")]
            if len(row):
                add(f"4.4_estimable_{pr}", "4.4", printed, 0.0005,
                    lambda row=row: row.iloc[0]["frac_estimable"])
        add("4.4_null_estimable", "4.4", "48.7%", 0.0005,
            lambda: v3["null_frac_estimable"].mean())
        _oth = v3[v3.trait != "ASD"]
        C.append(("3.7_others_reject", "3.7", "10⁻¹²⁰", "BOUND",
                  lambda: _oth["ks_p"].max()))
        add("3.7_weakest_other", "3.7", "2 × 10⁻¹²¹", 5e-122,
            lambda: _oth["ks_p"].max())
        C.append(("3.7_n_underflow", "3.7", "8", "BOUND",
                  lambda: float(abs(int((_oth["ks_p"] == 0).sum()) - 8))))
        C.append(("Abs_others_reject", "3.7", "10^−100^", "BOUND",
                  lambda: _oth["ks_p"].max()))

    v2 = D("verification/V2_gate_null_baseline.csv")
    if v2 is not None:
        w = v2[(v2.pair == "P4") & (v2.trait == "ASD")]
        if len(w):
            w0 = w.iloc[0]
            add("3.7_enrich_corrected", "3.7", "1.00", 5e-3,
                lambda w0=w0: w0["new_enrichment"])
            add("3.7_null_cond_rate", "3.7", "0.103", 5e-4,
                lambda w0=w0: w0["null_pass_rate_cond"])

    p0 = os.path.join(RES, "verification", "V0_PASS.txt")
    if os.path.exists(p0):
        C.append(("2.6_V0_exact", "2.6", "0", "BOUND",
                  lambda: 0.0 if open(p0).read().strip() == "TRUE" else 1.0))

    info = os.path.join(ROOT, "results", "ldsc_panel", "multipair.info.txt")
    if os.path.exists(info):
        inf = pd.read_csv(info, sep="\t").set_index("phenotype")
        for ph, cs2, ct in (("MDD", "58,854", "106,796"),
                            ("MDD_noUKB", "45,396", "97,250"),
                            ("BIP", "41,917", "371,549"),
                            ("BIP_noUKB", "40,463", "313,436"),
                            ("ASD", "18,381", "27,969"),
                            ("SCZ", "53,386", "77,258")):
            if ph in inf.index:
                add(f"prov_{ph}_cases", "3.3", cs2, 0,
                    lambda ph=ph: inf.loc[ph, "cases"])
                add(f"prov_{ph}_controls", "3.3", ct, 0,
                    lambda ph=ph: inf.loc[ph, "controls"])
        for ph, diff in (("MDD", "23,004"), ("BIP", "59,567")):
            if ph in inf.index and f"{ph}_noUKB" in inf.index:
                add(f"prov_{ph}_diff", "3.3", diff, 0,
                    lambda ph=ph: (inf.loc[ph, "cases"] + inf.loc[ph, "controls"])
                    - (inf.loc[f"{ph}_noUKB", "cases"]
                       + inf.loc[f"{ph}_noUKB", "controls"]))

    add("2.2_p3_drops", "2.2", "1,407", 0,
        lambda: int(D("panel/P3/attrition_audit.csv")
                    .query("processed")["p_ASD"].isna().sum()))
    add("2.2_p5_drops", "2.2", "1,474", 0,
        lambda: int(D("panel/P5/attrition_audit.csv")
                    .query("processed")["p_ASD"].isna().sum()))

    ad4 = D("applied_tierB/attrition_audit.csv")
    if ad4 is not None:
        from scipy import stats as _st
        pro4 = ad4[ad4["processed"] == True]
        b4 = pro4[pro4["p_ASD"].notna() & pro4["p_SCZ"].notna()]
        tab = pd.crosstab(b4["p_ASD"] < 0.05, b4["p_SCZ"] < 0.05).values
        add("3.6_fisher_or", "3.6", "2.37", 5e-3,
            lambda tab=tab: _st.contingency.odds_ratio(tab).statistic)
        add("3.6_fisher_ci_lo", "3.6", "1.32", 5e-3,
            lambda tab=tab: _st.contingency.odds_ratio(tab)
            .confidence_interval(0.95).low)
        add("3.6_fisher_ci_hi", "3.6", "4.53", 5e-3,
            lambda tab=tab: _st.contingency.odds_ratio(tab)
            .confidence_interval(0.95).high)
        add("3.6_fisher_p", "3.6", "0.0018", 5e-5,
            lambda tab=tab: _st.fisher_exact(tab)[1])
    we = J("worked_example/p1_worked_example.json")
    if we:
        add("3.5.1_n_reporting", "3.5.1", "1,369", 0, lambda: we["n_reporting"])
        add("3.5.1_n_floor", "3.5.1", "353", 0, lambda: we["n_at_floor"])
        add("3.5.1_pct_floor", "3.5.1", "25.8%", 0.0005,
            lambda: we["pct_at_floor"] / 100.0)
        WE = {
            "A": [("locus", "874", 0), ("K", "357", 0),
                  ("p_weak", "0.01138", 5e-6),
                  ("rho_hat", "0.2588", 5e-5), ("theta1", "0.0800", 5e-5),
                  ("theta2", "0.5660", 5e-5), ("slope", "0.7217", 5e-5),
                  ("intercept", "0.0632", 5e-5),
                  ("rho_corrected", "0.2711", 5e-5),
                  ("ci_width", "0.9926", 5e-5),
                  ("ci_width_corrected", "1.3754", 5e-5),
                  ("naive_divide_only", "0.359", 5e-4)],
            "B": [("locus", "38", 0), ("K", "285", 0),
                  ("rho_hat", "0.4922", 5e-5), ("theta1", "0.3192", 5e-5),
                  ("theta2", "0.3913", 5e-5), ("slope", "0.9881", 5e-5),
                  ("intercept", "0.0064", 5e-5),
                  ("rho_corrected", "0.4916", 5e-5),
                  ("ci_width", "0.8553", 5e-5),
                  ("ci_width_corrected", "0.8656", 5e-5)],
            "C": [("locus", "168", 0), ("K", "345", 0),
                  ("p_weak", "0.03097", 5e-6),
                  ("rho_hat", "0.6064", 5e-5), ("theta1", "0.0000", 5e-5),
                  ("theta2", "0.3509", 5e-5)],
        }
        for _lab, _rows in WE.items():
            for _fld, _pr, _tol in _rows:
                add(f"WE_{_lab}_{_fld}", "3.5.1", _pr, _tol,
                    lambda _lab=_lab, _fld=_fld: abs(we[_lab][_fld]))

    return C


_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def parse_printed(s):
    """'1.14%' -> 0.0114 ; '2,400' -> 2400 ; '−1.67' -> -1.67 ; '2 × 10⁻¹²' -> 2e-12"""
    s = s.strip().replace(",", "").replace("−", "-")
    m = re.fullmatch(r"([0-9.]+)\s*×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", s)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2).translate(_SUP))
    m = re.fullmatch(r"10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", s)
    if m:
        return 10.0 ** int(m.group(1).translate(_SUP))
    m = re.fullmatch(r"10\^(-?\d+)\^", s)
    if m:
        return 10.0 ** int(m.group(1))
    s = s.replace("×", "")
    pct = s.endswith("%")
    if pct:
        s = s[:-1]
    v = float(s)
    return v / 100 if pct else v


def present(printed, text):
    """Is `printed` in the manuscript as a WHOLE number, not a fragment?"""
    return re.search(r"(?<![\d.])" + re.escape(printed) + r"(?!\d)",
                     text) is not None


def main(verbose=False):
    if not os.path.exists(MS):
        print(f"manuscript not found: {MS}")
        return 2
    text = open(MS, encoding="utf-8").read()

    rows = []
    for cid, sec, printed, tol, fn in claims():
        try:
            got = fn()
        except Exception as e:
            got = None
            if verbose:
                print(f"  ({cid}: getter raised {e!r})")
        if got is None:
            rows.append((cid, sec, printed, None, "PENDING", ""))
            continue
        want = parse_printed(printed)
        num_ok = (float(got) <= want) if tol == "BOUND" \
            else abs(float(got) - want) <= tol
        in_ms = present(printed, text)
        status = "ok" if (num_ok and in_ms) else "FAIL"
        why = "" if status == "ok" else (
            ("value differs; " if not num_ok else "") +
            ("printed string not found in manuscript" if not in_ms else ""))
        rows.append((cid, sec, printed, got, status, why))

    n_ok = sum(r[4] == "ok" for r in rows)
    n_fail = sum(r[4] == "FAIL" for r in rows)
    n_pend = sum(r[4] == "PENDING" for r in rows)

    for cid, sec, printed, got, status, why in rows:
        if status != "ok" or verbose:
            g = "-" if got is None else f"{float(got):.6g}"
            print(f"[{status:7s}] {sec:4s} {cid:36s} printed={printed:>10s} "
                  f"computed={g:>12s} {why}")

    print(f"\n{n_ok} passed, {n_fail} failed, {n_pend} pending "
          f"({len(rows)} claims checked).")

    ms_pending = len(re.findall(r"PENDING", text))
    print(f"manuscript contains {ms_pending} 'PENDING' marker(s).")
    if n_pend and ms_pending == 0:
        print("INCONSISTENT: claims are unavailable but the manuscript marks nothing "
              "as PENDING. A number may have been filled in by hand.")
        return 1
    return 1 if n_fail else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    sys.exit(main(ap.parse_args().verbose))
