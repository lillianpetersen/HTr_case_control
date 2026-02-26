"""
Unit tests for filtering.py: filter_counts.
"""

import pandas as pd
import pytest

from htr_hit_caller.filtering import filter_counts


def _make_df(**kwargs):
    """Make a single-row DataFrame with provided counts."""
    defaults = {
        "counts OFF R1": 1000,
        "counts ON R1":  1000,
        "counts OFF R2": 1000,
        "counts ON R2":  1000,
        "extra_col": "keep_me",
    }
    defaults.update(kwargs)
    return pd.DataFrame([defaults], index=["dom0"])


class TestFilterCounts:

    def test_keeps_sufficient_domain(self):
        df = _make_df()
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert "dom0" in out.index

    def test_removes_low_bucket(self):
        """A domain with one bucket below min_per_bucket is dropped."""
        df = _make_df(**{"counts ON R1": 100})  # below 250
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert len(out) == 0

    def test_removes_low_total_r1(self):
        """R1 total (ON+OFF) below min_per_replicate_total → drop."""
        # ON=300, OFF=300 → total=600 < 700
        df = _make_df(**{"counts OFF R1": 300, "counts ON R1": 300})
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert len(out) == 0

    def test_removes_low_total_r2(self):
        """R2 total below threshold → drop, even if R1 is fine."""
        df = _make_df(**{"counts OFF R2": 300, "counts ON R2": 300})
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert len(out) == 0

    def test_keeps_domain_exactly_at_threshold(self):
        """Domain exactly at min values should be retained (not strictly less)."""
        df = _make_df(
            **{
                "counts OFF R1": 250, "counts ON R1": 450,
                "counts OFF R2": 250, "counts ON R2": 450,
            }
        )
        # total R1 = 700 (= min), each bucket = 250 (= min)
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert "dom0" in out.index

    def test_preserves_extra_columns(self):
        """Non-count columns should be preserved in output."""
        df = _make_df()
        out = filter_counts(df)
        assert "extra_col" in out.columns

    def test_multiple_domains(self):
        """Only failing domains are dropped; passing ones are kept."""
        df = pd.DataFrame({
            "counts OFF R1": [1000,  100, 1000],
            "counts ON R1":  [1000, 1000, 1000],
            "counts OFF R2": [1000, 1000, 1000],
            "counts ON R2":  [1000, 1000, 1000],
        }, index=["good", "bad_bucket", "good2"])
        out = filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
        assert "good" in out.index
        assert "good2" in out.index
        assert "bad_bucket" not in out.index

    def test_does_not_modify_input(self):
        df = _make_df()
        original_len = len(df)
        _ = filter_counts(df)
        assert len(df) == original_len
