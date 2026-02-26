"""
Command-line interface for htr_hit_caller.

Usage:
    htr-hit-caller --case ERKi_combo.csv --control Veh_combo.csv \
        --null-col is_null --output-dir results/ [options]
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from .hit_calling import call_hits
from .plots import make_all_plots


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="htr-hit-caller",
        description="Bootstrap hit-calling for HT-recruit kinase inhibitor screens",
    )
    parser.add_argument("--case",    required=True, help="CSV file for case condition (treatment)")
    parser.add_argument("--control", required=True, help="CSV file for control condition (DMSO)")
    parser.add_argument("--null-col", required=True,
                        help="Boolean column name marking null/non-activating domains")
    parser.add_argument("--output-dir", default="htr_results",
                        help="Directory to write results CSV and plots (default: htr_results)")
    parser.add_argument("--case-name",    default="case")
    parser.add_argument("--control-name", default="control")
    parser.add_argument("--n-boot",         type=int,   default=1000)
    parser.add_argument("--thresh",         type=float, default=5.0)
    parser.add_argument("--add-counts",     type=int,   default=100)
    parser.add_argument("--fpr-strong",     type=float, default=0.005)
    parser.add_argument("--fpr-weak",       type=float, default=0.02)
    parser.add_argument("--min-per-bucket", type=int,   default=250)
    parser.add_argument("--min-per-replicate-total", type=int, default=700)
    parser.add_argument("--random-state",   type=int,   default=None)
    parser.add_argument("--no-plots", action="store_true", help="Skip saving plots")

    args = parser.parse_args(argv)

    case_df    = pd.read_csv(args.case,    index_col=0)
    control_df = pd.read_csv(args.control, index_col=0)

    results, delta_df, thresholds = call_hits(
        case=case_df,
        control=control_df,
        null_col=args.null_col,
        n_boot=args.n_boot,
        thresh=args.thresh,
        add_counts=args.add_counts,
        fpr_strong=args.fpr_strong,
        fpr_weak=args.fpr_weak,
        min_per_bucket=args.min_per_bucket,
        min_per_replicate_total=args.min_per_replicate_total,
        random_state=args.random_state,
        case_name=args.case_name,
        control_name=args.control_name,
    )

    import os
    os.makedirs(args.output_dir, exist_ok=True)

    results_path = os.path.join(args.output_dir, "results.csv")
    results.to_csv(results_path)
    print(f"Results saved to {results_path}")
    print(f"  Strong up hits:   {results['sig_up'].sum()}")
    print(f"  Strong down hits: {results['sig_down'].sum()}")
    print(f"  Weak up hits:     {results['sig_up_weak'].sum()}")
    print(f"  Weak down hits:   {results['sig_down_weak'].sum()}")

    if not args.no_plots:
        make_all_plots(
            results, delta_df, thresholds,
            output_dir=args.output_dir,
            case_name=args.case_name,
            control_name=args.control_name,
            null_col=args.null_col,
        )
        print(f"Plots saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
