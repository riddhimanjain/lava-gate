"""15_build_table8.py -- assemble Table 8 from the panel outputs."""
import argparse
import glob
import json
import os

import pandas as pd

PAIRS = [("P4", "ASD", "SCZ"), ("P1", "SCZ", "BIP"), ("P2a", "MDD", "BIP"),
         ("P2b", "MDD_noUKB", "BIP_noUKB"), ("P3", "ASD", "ADHD"), ("P5", "AN", "ASD")]
PENDING = "`PENDING`"


def fmt(v, spec="{:.4g}", pct=None):
    if v is None:
        return PENDING
    try:
        if pct is not None:
            return f"{v:,.0f} ({100*v/pct:.1f}%)"
        return spec.format(v)
    except (TypeError, ValueError):
        return PENDING


def load(panel_dir, group_dir):
    """Per-pair summaries and the R1 group rows, keyed by pair label."""
    summaries, rows = {}, {}
    for f in glob.glob(os.path.join(group_dir, "PANEL_SUMMARY_*.csv")):
        for r in pd.read_csv(f).to_dict("records"):
            rows[r["pair"]] = r
    for pair, _, _ in PAIRS:
        p = os.path.join(panel_dir, pair, "pair_summary.json")
        if os.path.exists(p):
            summaries[pair] = json.load(open(p))
    return summaries, rows


def build(panel_dir, group_dir):
    summ, rows = load(panel_dir, group_dir)
    labels = [p for p, _, _ in PAIRS]

    def cell(pair, fn):
        try:
            v = fn(summ.get(pair), rows.get(pair))
            return PENDING if v is None else v
        except (KeyError, TypeError, AttributeError, ZeroDivisionError):
            return PENDING

    def weaker_stronger(r):
        """Channel-2 loss for the weaker and the stronger trait, in that order."""
        if r is None:
            return None, None
        a, b = r["ch2_loss_1"], r["ch2_loss_2"]
        w = r["weaker_trait"]
        return (a, b) if w == r["traits"].split(" x ")[0] else (b, a)

    spec = [
        ("Processed", lambda s, r: f"{r['processed']:,} ({100*r['processed']/r['blocks']:.1f}%)"),
        ("Weaker trait estimable", lambda s, r: (
            f"{min(r['est_1'], r['est_2']):,} "
            f"({100*min(r['est_1'], r['est_2'])/r['processed']:.1f}%)")),
        ("Both estimable", lambda s, r: (
            f"{r['both_estimable']:,} ({100*r['both_estimable']/r['processed']:.1f}%)")),
        ("**Reported**", lambda s, r: f"**{r['reported']:,} ({r['pct_reported']:.2f}%)**"),
        ("Bonferroni survivors", lambda s, r: f"{r['bonferroni_survivors']:,}"),
        ("Ch. 2 loss, weaker", lambda s, r: f"{weaker_stronger(r)[0]:.1f}%"),
        ("Ch. 2 loss, stronger", lambda s, r: f"{weaker_stronger(r)[1]:.1f}%"),
        ("SNPs in pass", lambda s, r: f"{r['n_snps_common']:,}"),
        ("Median *K*", lambda s, r: f"{r['median_K']:.0f}"),
        ("Measured *r*~e~", lambda s, r: f"{s['r_e']:.4f}"),
        ("Predicted intercept", lambda s, r: f"{s['predicted_gate_intercept']:.4f}"),
        ("Fitted intercept", lambda s, r: f"{s['fitted_aggregate_intercept']:.4f}"),
        ("Aggregate slope", lambda s, r: f"{s['aggregate_slope']:.3f}"),
        ("C1 slope fold-change", lambda s, r: (
            f"{s['diagnostic']['C1_stability']['slope_fold_change']:.2f}"
            f"{'' if s['diagnostic']['C1_stability']['passed'] else ' ✗'}")),
        ("C2 pooled *z*, weaker", lambda s, r: (
            f"{min(s['diagnostic']['C2_floor']['pooled_z'].values()):+.1f}")),
        ("C3 null fraction, weaker", lambda s, r: (
            f"{max(s['diagnostic']['C3_signal']['frac_null'].values()):.2f}")),
        ("**Usable?**", lambda s, r: f"**{'Yes' if s['diagnostic']['usable'] else 'No'}**"),
    ]

    md = ["| Quantity | " + " | ".join(f"**{l}**" for l in labels) + " |",
          "|---|" + "---|" * len(labels)]
    table = {}
    for name, fn in spec:
        cells = [cell(p, fn) for p in labels]
        table[name] = dict(zip(labels, cells))
        md.append(f"| {name} | " + " | ".join(cells) + " |")

    n_pending = sum(v == PENDING for row in table.values() for v in row.values())
    return "\n".join(md), table, n_pending


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--group-summaries", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    md, table, n_pending = build(a.panel, a.group_summaries)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(md + "\n")
    with open(os.path.splitext(a.out)[0] + ".json", "w") as f:
        json.dump(table, f, indent=2)
    print(md)
    print(f"\n{n_pending} cell(s) still PENDING.")
    if n_pending:
        print("Those pairs have not been run. Do not fill them in by hand.")
