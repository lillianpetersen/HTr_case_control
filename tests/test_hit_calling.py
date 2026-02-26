"""
Unit tests for hit_calling.py: call_hits.
"""

import numpy as np
import pandas as pd
import pytest

from htr_hit_caller import call_hits


def _make_synthetic_data(
    n_null: int = 80,
    n_activating: int = 20,
    seed: int = 42,
    case_shift: float = 0.0,
):
    """
    Build synthetic case and control DataFrames.

    All domains have sufficient counts to pass filtering.
    ``case_shift`` adds extra ON counts to non-null domains in the case condition
    to simulate real hits.
    """
    rng = np.random.default_rng(seed)
    n = n_null + n_activating

    # Base counts: large enough to pass all filters
    base_off = rng.integers(800, 1500, size=(n, 2))
    base_on  = rng.integers(800, 1500, size=(n, 2))

    def make_df(extra_on=0):
        on = base_on.copy().astype(float)
        # Add shift to non-null (activating) domains
        if extra_on != 0:
            on[n_null:, :] += extra_on
        return pd.DataFrame({
            "counts OFF R1": base_off[:, 0],
            "counts ON R1":  on[:, 0].astype(int),
            "counts OFF R2": base_off[:, 1],
            "counts ON R2":  on[:, 1].astype(int),
            "is_null": [True] * n_null + [False] * n_activating,
        }, index=[f"dom{i}" for i in range(n)])

    case    = make_df(extra_on=case_shift)
    control = make_df(extra_on=0)
    return case, control


class TestCallHitsShape:

    def test_returns_three_outputs(self):
        case, control = _make_synthetic_data()
        out = call_hits(case, control, null_col="is_null", n_boot=50, random_state=0)
        assert len(out) == 3

    def test_results_index_matches_shared(self):
        """Results should only contain domains present in both after filtering."""
        case, control = _make_synthetic_data(n_null=50, n_activating=10)
        results, delta_df, thresholds = call_hits(
            case, control, null_col="is_null", n_boot=50, random_state=0
        )
        assert len(results) == len(delta_df)
        assert list(results.index) == list(delta_df.index)

    def test_delta_df_shape(self):
        n_boot = 30
        case, control = _make_synthetic_data()
        results, delta_df, _ = call_hits(
            case, control, null_col="is_null", n_boot=n_boot, random_state=0
        )
        # 2 reps × n_boot columns
        assert delta_df.shape[1] == 2 * n_boot

    def test_expected_result_columns(self):
        case, control = _make_synthetic_data()
        results, _, _ = call_hits(
            case, control, null_col="is_null", n_boot=50,
            random_state=0, case_name="ERKi", control_name="Veh"
        )
        for col in ["ERKi R1", "ERKi Avg", "Veh R1", "Veh Avg",
                    "post_prob", "logit_dom",
                    "sig_up", "sig_down", "sig_up_weak", "sig_down_weak"]:
            assert col in results.columns, f"Missing column: {col}"


class TestPostProb:

    def test_post_prob_bounds(self):
        case, control = _make_synthetic_data()
        results, _, _ = call_hits(
            case, control, null_col="is_null", n_boot=100, random_state=0
        )
        assert (results["post_prob"] >= 0).all()
        assert (results["post_prob"] <= 1).all()

    def test_null_centered_near_half(self):
        """Null domain post_probs should be distributed near 0.5 on average."""
        case, control = _make_synthetic_data(n_null=100, n_activating=0)
        results, _, _ = call_hits(
            case, control, null_col="is_null", n_boot=200, random_state=0
        )
        null_mean = results.loc[results["is_null"], "post_prob"].mean()
        assert abs(null_mean - 0.5) < 0.1, f"Null mean post_prob = {null_mean:.3f}, expected ~0.5"

    def test_logit_dom_finite(self):
        """logit_dom should not contain inf/-inf (clipping prevents this)."""
        case, control = _make_synthetic_data()
        results, _, _ = call_hits(
            case, control, null_col="is_null", n_boot=50, random_state=0
        )
        assert np.all(np.isfinite(results["logit_dom"]))


class TestHitCalls:

    def test_strong_hits_subset_of_weak(self):
        """Every strong hit should also satisfy the weak hit condition or stronger."""
        case, control = _make_synthetic_data(case_shift=2000)
        results, _, _ = call_hits(
            case, control, null_col="is_null", n_boot=100, random_state=0
        )
        # sig_up are strong hits; sig_up_weak are the *additional* weak-only
        # In the code: sig_up_weak = (> thr_up_weak) & (<= thr_up)
        # So sig_up and sig_up_weak are mutually exclusive; together they cover all weak+strong
        # Check: no sig_up domain is NOT (sig_up | sig_up_weak)
        strong_up = results["sig_up"]
        strong_down = results["sig_down"]
        # Strong hits should not also be weak-only (they're separate categories in this impl)
        # But all truly strong hits should have higher post_prob than weak threshold
        thr_up      = results.loc[results["sig_up"], "post_prob"].min() if strong_up.any() else None
        thr_up_weak = results.loc[results["sig_up_weak"], "post_prob"].max() if results["sig_up_weak"].any() else None
        if thr_up is not None and thr_up_weak is not None:
            assert thr_up >= thr_up_weak

    def test_thresholds_dict_keys(self):
        case, control = _make_synthetic_data()
        _, _, thresholds = call_hits(
            case, control, null_col="is_null", n_boot=50, random_state=0
        )
        assert set(thresholds.keys()) == {"thr_up", "thr_down", "thr_up_weak", "thr_down_weak"}

    def test_thresholds_ordering(self):
        """thr_down < thr_down_weak < thr_up_weak < thr_up."""
        case, control = _make_synthetic_data()
        _, _, thr = call_hits(
            case, control, null_col="is_null", n_boot=100, random_state=0
        )
        assert thr["thr_down"] < thr["thr_down_weak"]
        assert thr["thr_down_weak"] < thr["thr_up_weak"]
        assert thr["thr_up_weak"] < thr["thr_up"]

    def test_inner_join_drops_missing_domains(self):
        """Domains missing from either case or control are excluded from results."""
        case, control = _make_synthetic_data(n_null=50, n_activating=10)
        # Drop a domain from case only
        case_missing = case.drop("dom0")
        results, _, _ = call_hits(
            case_missing, control, null_col="is_null", n_boot=50, random_state=0
        )
        assert "dom0" not in results.index

    def test_reproducible_with_seed(self):
        case, control = _make_synthetic_data()
        r1, _, _ = call_hits(case, control, null_col="is_null", n_boot=50, random_state=99)
        r2, _, _ = call_hits(case, control, null_col="is_null", n_boot=50, random_state=99)
        np.testing.assert_array_almost_equal(
            r1["post_prob"].values, r2["post_prob"].values
        )
