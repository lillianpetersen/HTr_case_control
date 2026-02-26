"""
Example usage of htr_hit_caller with real LP012 kinase inhibitor screen data.

This script demonstrates the full pipeline:
1. Load the pre-processed count data (from kinase_inhibition_screen_results.csv)
2. Define the null set (non-activating domains)
3. Run call_hits()
4. Print summary statistics and save plots

Run from the htr_hit_caller directory after installing:
    pip install -e ".[dev]"
    python examples/example_usage.py
"""

import os
import pandas as pd
import numpy as np

import htr_hit_caller as htr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR = os.path.expanduser(
    "~/Library/CloudStorage/GoogleDrive-lilliank@stanford.edu/My Drive/project"
)
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "LP012_kinase_inh_HTr_analysis")
OUTPUT_DIR  = os.path.join(PROJECT_DIR, "figures",  "LP012_kinase_inh_HTr_analysis", "htr_pkg_output")

HITS_CSV = os.path.join(RESULTS_DIR, "kinase_inhibition_screen_results.csv")

# ---------------------------------------------------------------------------
# Load pre-processed data (already filtered, renamed, labeled)
# ---------------------------------------------------------------------------
print("Loading data...")
df = pd.read_csv(HITS_CSV, index_col=0)

# Remove the 8 known duplicate index entries
df = df[~df.index.duplicated(keep="first")]
print(f"  Total domains: {len(df)}")

# Extract per-condition count DataFrames
veh = df[[
    "counts OFF R1 Veh", "counts ON R1 Veh",
    "counts OFF R2 Veh", "counts ON R2 Veh",
]].rename(columns={
    "counts OFF R1 Veh": "counts OFF R1",
    "counts ON R1 Veh":  "counts ON R1",
    "counts OFF R2 Veh": "counts OFF R2",
    "counts ON R2 Veh":  "counts ON R2",
}).copy()

erki = df[[
    "counts OFF R1 ERKi", "counts ON R1 ERKi",
    "counts OFF R2 ERKi", "counts ON R2 ERKi",
]].rename(columns={
    "counts OFF R1 ERKi": "counts OFF R1",
    "counts ON R1 ERKi":  "counts ON R1",
    "counts OFF R2 ERKi": "counts OFF R2",
    "counts ON R2 ERKi":  "counts ON R2",
}).copy()

# ---------------------------------------------------------------------------
# Define null set: non-activating domains
# ---------------------------------------------------------------------------
activating = df["activating"].fillna(False).astype(bool)
veh["is_null"]  = ~activating
erki["is_null"] = ~activating

n_null = veh["is_null"].sum()
n_act  = (~veh["is_null"]).sum()
print(f"  Null (non-activating): {n_null} domains")
print(f"  Activating:            {n_act}  domains")

# ---------------------------------------------------------------------------
# Run hit calling
# ---------------------------------------------------------------------------
print("\nRunning bootstrapping hit caller (n_boot=1000)...")
results, delta_df, thresholds = htr.call_hits(
    case=erki,
    control=veh,
    null_col="is_null",
    n_boot=1000,
    fpr_strong=0.005,
    fpr_weak=0.02,
    min_per_bucket=250,
    min_per_replicate_total=700,
    random_state=42,
    case_name="ERKi",
    control_name="Veh",
)

print(f"  Domains after filtering: {len(results)}")
print(f"  Strong up hits:   {results['sig_up'].sum()}")
print(f"  Weak up hits:     {results['sig_up_weak'].sum()}")
print(f"  Strong down hits: {results['sig_down'].sum()}")
print(f"  Weak down hits:   {results['sig_down_weak'].sum()}")
print(f"\n  Thresholds: {thresholds}")

# ---------------------------------------------------------------------------
# Check known hits
# ---------------------------------------------------------------------------
print("\nChecking known hits...")
# ELK1;32: ERK-dependent transcription factor → less active under ERKi (sig_down)
known_hits = {
    "ELK1;32": "sig_down",
    "TLE3;60": "sig_down",
}
for domain, expected_hit_type in known_hits.items():
    if domain in results.index:
        pp = results.loc[domain, "post_prob"]
        is_hit = results.loc[domain, expected_hit_type]
        print(f"  {domain}: post_prob={pp:.4f}, {expected_hit_type}={is_hit}")
    else:
        print(f"  {domain}: filtered out")

# ---------------------------------------------------------------------------
# Save plots
# ---------------------------------------------------------------------------
print(f"\nSaving plots to {OUTPUT_DIR} ...")
htr.make_all_plots(
    results=results,
    delta_df=delta_df,
    thresholds=thresholds,
    output_dir=OUTPUT_DIR,
    case_name="ERKi",
    control_name="Veh",
    null_col="is_null",
)
print("Done.")

# ---------------------------------------------------------------------------
# Save results CSV
# ---------------------------------------------------------------------------
out_csv = os.path.join(OUTPUT_DIR, "ERKi_hit_calling_results.csv")
results.to_csv(out_csv)
print(f"Results saved to {out_csv}")
