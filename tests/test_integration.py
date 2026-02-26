"""
Integration test against real LP012 data.

Skipped automatically if data files are not present.

Uses the pre-processed kinase_inhibition_screen_results.csv which already
contains count columns in the correct format (counts OFF R1 Veh, etc.)
and the `activating` annotation column.
"""

import os
import pytest
import numpy as np
import pandas as pd

import htr_hit_caller as htr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_BASE = os.path.expanduser(
    "~/Library/CloudStorage/GoogleDrive-lilliank@stanford.edu"
    "/My Drive/project"
)
_RESULTS_DIR = os.path.join(_BASE, "results", "LP012_kinase_inh_HTr_analysis")
_HITS_CSV = os.path.join(_RESULTS_DIR, "kinase_inhibition_screen_results.csv")

_DATA_AVAILABLE = os.path.exists(_HITS_CSV)

_VEH_COUNT_RENAME = {
    "counts OFF R1 Veh": "counts OFF R1",
    "counts ON R1 Veh":  "counts ON R1",
    "counts OFF R2 Veh": "counts OFF R2",
    "counts ON R2 Veh":  "counts ON R2",
}
_ERKI_COUNT_RENAME = {
    "counts OFF R1 ERKi": "counts OFF R1",
    "counts ON R1 ERKi":  "counts ON R1",
    "counts OFF R2 ERKi": "counts OFF R2",
    "counts ON R2 ERKi":  "counts ON R2",
}


@pytest.fixture(scope="module")
def lp012_data():
    """Load and prepare LP012 DataFrames from the pre-processed hits CSV."""
    df = pd.read_csv(_HITS_CSV, index_col=0)

    # Drop duplicate index entries (8 known duplicates in this dataset)
    df = df[~df.index.duplicated(keep="first")]

    # Extract per-condition count DataFrames
    veh = df[list(_VEH_COUNT_RENAME.keys())].rename(columns=_VEH_COUNT_RENAME).copy()
    erki = df[list(_ERKI_COUNT_RENAME.keys())].rename(columns=_ERKI_COUNT_RENAME).copy()

    # Mark null domains: non-activating per the reference results
    activating = df["activating"].fillna(False).astype(bool)
    veh["is_null"]  = ~activating
    erki["is_null"] = ~activating

    return erki, veh, df


@pytest.mark.skipif(not _DATA_AVAILABLE, reason="LP012 data files not found")
class TestIntegrationLP012:

    def test_runs_without_error(self, lp012_data):
        erki, veh, _ = lp012_data
        results, delta_df, thresholds = htr.call_hits(
            case=erki, control=veh, null_col="is_null",
            n_boot=2000, random_state=42,
            min_per_bucket=250, min_per_replicate_total=700,
            case_name="ERKi", control_name="Veh",
        )
        assert results is not None
        assert len(results) > 0

    def test_reasonable_hit_counts(self, lp012_data):
        erki, veh, _ = lp012_data
        results, _, _ = htr.call_hits(
            case=erki, control=veh, null_col="is_null",
            n_boot=2000, random_state=42,
            min_per_bucket=250, min_per_replicate_total=700,
            case_name="ERKi", control_name="Veh",
        )
        n_up   = results["sig_up"].sum()
        n_down = results["sig_down"].sum()
        assert 10 <= n_up   <= 60,  f"Expected 10-60 strong up hits, got {n_up}"
        assert 10 <= n_down <= 80,  f"Expected 10-80 strong down hits, got {n_down}"

    def test_elk1_32_is_hit(self, lp012_data):
        """ELK1;32 is a known ERKi-sensitive down hit (less active under ERK inhibition)."""
        erki, veh, _ = lp012_data
        results, _, _ = htr.call_hits(
            case=erki, control=veh, null_col="is_null",
            n_boot=2000, random_state=42,
            min_per_bucket=250, min_per_replicate_total=700,
            case_name="ERKi", control_name="Veh",
        )
        assert "ELK1;32" in results.index, "ELK1;32 was filtered out"
        assert results.loc["ELK1;32", "sig_down"], \
            f"ELK1;32 post_prob = {results.loc['ELK1;32', 'post_prob']:.4f}"

    def test_null_post_prob_approximately_uniform(self, lp012_data):
        """Null domain post_probs should be spread across [0,1] with mean ~0.5."""
        erki, veh, _ = lp012_data
        results, _, _ = htr.call_hits(
            case=erki, control=veh, null_col="is_null",
            n_boot=2000, random_state=42,
            min_per_bucket=250, min_per_replicate_total=700,
            case_name="ERKi", control_name="Veh",
        )
        null_pp = results.loc[results["is_null"], "post_prob"]
        assert abs(null_pp.mean() - 0.5) < 0.05, \
            f"Null mean post_prob = {null_pp.mean():.3f}"
        # Std of Uniform[0,1] ≈ 0.289; allow some slack
        assert null_pp.std() > 0.15, \
            f"Null post_prob std = {null_pp.std():.3f} (expected spread across [0,1])"

    def test_post_prob_agrees_with_original_within_tolerance(self, lp012_data):
        """
        post_prob values should correlate well with the reference script output.
        Both use bootstrapping from the same data; results will differ due to
        different random seeds but should be highly correlated.
        """
        erki, veh, ref_df = lp012_data
        results, _, _ = htr.call_hits(
            case=erki, control=veh, null_col="is_null",
            n_boot=1000, random_state=42,
            min_per_bucket=250, min_per_replicate_total=700,
            case_name="ERKi", control_name="Veh",
        )
        ref_col = "post prob ERKi"
        if ref_col not in ref_df.columns:
            pytest.skip(f"Reference column '{ref_col}' not in hits CSV")

        shared = results.index.intersection(ref_df.index)
        ref  = ref_df.loc[shared, ref_col].dropna()
        ours = results.loc[ref.index, "post_prob"]

        # Should be highly correlated (different seeds, same underlying data)
        corr = np.corrcoef(ours.values, ref.values)[0, 1]
        assert corr > 0.95, f"Correlation with reference post_prob = {corr:.3f} (expected > 0.95)"
