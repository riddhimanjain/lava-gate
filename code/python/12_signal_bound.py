"""12_signal_bound.py -- an upper bound on the fraction of loci carrying real local"""
import json
import numpy as np
import pandas as pd

ALPHA = 0.05
aud = pd.read_csv("out_R/attrition_audit.csv")
pro = aud[aud["processed"] == True]
both = pro[np.isfinite(pro["p_ASD"]) & np.isfinite(pro["p_SCZ"])]

out = {"alpha": ALPHA, "n_processed": int(len(pro)), "n_both_estimable": int(len(both))}
for trait in ("ASD", "SCZ"):
    p = both[f"p_{trait}"]
    obs = float((p < ALPHA).mean())
    exp_null_passes = ALPHA * len(both)
    frac_passes_null = min(1.0, exp_null_passes / max((p < ALPHA).sum(), 1))
    out[trait] = dict(
        pass_rate=obs,
        n_pass=int((p < ALPHA).sum()),
        expected_null_passes=float(exp_null_passes),
        fraction_of_passes_expected_null=float(frac_passes_null),
        enrichment_over_alpha=float(obs / ALPHA),
        pi_upper_bound_perfect_power=float(max(0.0, (obs - ALPHA) / (1 - ALPHA))),
    )
    print(f"  {trait}: pass={obs:.4f} ({int((p<ALPHA).sum())} loci), "
          f"{frac_passes_null:.0%} of passes expected null under a complete null, "
          f"enrichment={obs/ALPHA:.2f}x")

n_joint = int(((both["p_ASD"] < ALPHA) & (both["p_SCZ"] < ALPHA)).sum())
out["joint"] = dict(
    n_joint_pass=n_joint,
    expected_joint_if_ASD_all_null=float(ALPHA * (both["p_SCZ"] < ALPHA).sum()),
    note=("If every ASD locus were null, ASD would still pass at 5% of the loci "
          "where SCZ passes, producing this many joint passes by chance alone."),
)
print(f"  joint: observed {n_joint}; expected {out['joint']['expected_joint_if_ASD_all_null']:.1f} "
      f"if ASD carried no local signal at all")
json.dump(out, open("out/signal_bound.json", "w"), indent=2)
