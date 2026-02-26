"""
Unit tests for enrichment.py: compute_rhos and bootstrap_rhos.
"""

import numpy as np
import pandas as pd
import pytest

from htr_hit_caller.enrichment import compute_rhos, bootstrap_rhos


def _make_counts_df(off_r1, on_r1, off_r2, on_r2, index=None):
    """Helper: build a small counts DataFrame."""
    data = {
        "counts OFF R1": off_r1,
        "counts ON R1":  on_r1,
        "counts OFF R2": off_r2,
        "counts ON R2":  on_r2,
    }
    if index is None:
        index = [f"dom{i}" for i in range(len(off_r1))]
    return pd.DataFrame(data, index=index)


# ---------------------------------------------------------------------------
# compute_rhos
# ---------------------------------------------------------------------------

class TestComputeRhos:

    def test_basic_known_values(self):
        """
        With two domains and uniform counts, rho should be predictable.

        Domain A: ON=200, OFF=100 in both reps → ON/OFF enriched
        Domain B: ON=100, OFF=200 → depleted

        Library totals: ON_total = 300, OFF_total = 300 per rep
        rho_A = log2((200/300) / (100/300)) = log2(2) = 1.0
        rho_B = log2((100/300) / (200/300)) = log2(0.5) = -1.0
        """
        df = _make_counts_df(
            off_r1=[100, 200], on_r1=[200, 100],
            off_r2=[100, 200], on_r2=[200, 100],
        )
        result = compute_rhos(df, thresh=0)
        assert pytest.approx(result.loc["dom0", "rho_R1"], abs=1e-6) == 1.0
        assert pytest.approx(result.loc["dom1", "rho_R1"], abs=1e-6) == -1.0
        assert pytest.approx(result.loc["dom0", "rho_mean"], abs=1e-6) == 1.0

    def test_returns_expected_columns(self):
        df = _make_counts_df([100], [200], [100], [200])
        result = compute_rhos(df)
        assert set(result.columns) == {"rho_R1", "rho_R2", "rho_mean", "rho_sd"}

    def test_clips_thresh(self):
        """Counts below thresh should be clipped to thresh before normalization."""
        # Domain 0: ON=1 (below thresh=5), should be clipped to 5
        # Domain 1: ON=1000
        df = _make_counts_df(
            off_r1=[5, 5],   on_r1=[1, 1000],
            off_r2=[5, 5],   on_r2=[1, 1000],
        )
        result_clipped = compute_rhos(df, thresh=5)
        result_no_clip = compute_rhos(df, thresh=0)
        # With thresh=5, dom0 ON is clipped from 1→5; without thresh dom0 ON=1
        # The rho for dom0 should differ
        assert result_clipped.loc["dom0", "rho_R1"] != result_no_clip.loc["dom0", "rho_R1"]

    def test_sd_zero_when_replicates_equal(self):
        """When R1 and R2 are identical, SD should be 0."""
        df = _make_counts_df([100], [200], [100], [200])
        result = compute_rhos(df)
        assert pytest.approx(result.loc["dom0", "rho_sd"], abs=1e-9) == 0.0

    def test_index_preserved(self):
        idx = ["geneA;1", "geneB;2", "geneC;3"]
        df = _make_counts_df([100]*3, [200]*3, [100]*3, [200]*3, index=idx)
        result = compute_rhos(df)
        assert list(result.index) == idx


# ---------------------------------------------------------------------------
# bootstrap_rhos
# ---------------------------------------------------------------------------

class TestBootstrapRhos:

    def _make_df(self, n=20, seed=0):
        rng = np.random.default_rng(seed)
        counts = rng.integers(500, 2000, size=(n, 4))
        return _make_counts_df(
            off_r1=counts[:, 0], on_r1=counts[:, 1],
            off_r2=counts[:, 2], on_r2=counts[:, 3],
        )

    def test_output_shape(self):
        df = self._make_df(n=15)
        n_boot = 50
        result = bootstrap_rhos(df, n_boot=n_boot, random_state=0)
        # 2 replicates × n_boot
        assert result.shape == (15, 2 * n_boot)

    def test_bootstrap_mean_close_to_observed(self):
        """Bootstrap mean rho should be close to the observed rho."""
        df = self._make_df(n=30)
        n_boot = 2000
        boot = bootstrap_rhos(df, n_boot=n_boot, thresh=5, add_counts=0, random_state=7)
        obs = compute_rhos(df, thresh=5)

        # Compare bootstrap column mean (mean over all boot samples per domain)
        boot_mean = boot.mean(axis=1)
        obs_mean = obs["rho_mean"].values

        # Should be within 0.3 log2 units for well-sampled domains
        assert np.abs(boot_mean - obs_mean).mean() < 0.3

    def test_reproducible_with_seed(self):
        df = self._make_df(n=10)
        b1 = bootstrap_rhos(df, n_boot=20, random_state=42)
        b2 = bootstrap_rhos(df, n_boot=20, random_state=42)
        np.testing.assert_array_equal(b1, b2)

    def test_different_seeds_give_different_results(self):
        df = self._make_df(n=10)
        b1 = bootstrap_rhos(df, n_boot=20, random_state=1)
        b2 = bootstrap_rhos(df, n_boot=20, random_state=2)
        assert not np.allclose(b1, b2)

    def test_no_nan_or_inf(self):
        df = self._make_df(n=20)
        boot = bootstrap_rhos(df, n_boot=50, random_state=0)
        assert np.all(np.isfinite(boot))
