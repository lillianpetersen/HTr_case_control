"""
Main hit-calling pipeline for HT-recruit kinase inhibitor screens.
"""

from __future__ import annotations

from typing import Optional, Tuple, Dict

import numpy as np
import pandas as pd
from scipy.special import logit

from .enrichment import compute_rhos, bootstrap_rhos
from .filtering import filter_counts


def call_hits(
    case: pd.DataFrame,
    control: pd.DataFrame,
    null_col: str,
    n_boot: int = 1000,
    thresh: float = 5,
    add_counts: int = 100,
    fpr_strong: float = 0.005,
    fpr_weak: float = 0.02,
    min_per_bucket: int = 250,
    min_per_replicate_total: int = 700,
    random_state: Optional[int] = None,
    case_name: str = "case",
    control_name: str = "control",
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    """
    Full bootstrapping hit-calling pipeline.

    Compares gene activation domain activity between a kinase inhibitor
    condition (``case``) and a DMSO vehicle control (``control``) using
    a bootstrap posterior probability approach.

    Parameters
    ----------
    case : DataFrame
        Counts for the treatment condition. Must contain columns:
        'counts OFF R1', 'counts ON R1', 'counts OFF R2', 'counts ON R2'
        and a boolean column named ``null_col``.
    control : DataFrame
        Counts for the vehicle/control condition. Same column requirements.
    null_col : str
        Name of a boolean column in both DataFrames marking non-activating
        (null) domains used to calibrate hit thresholds.
    n_boot : int
        Number of bootstrap iterations per replicate per condition.
    thresh : float
        Floor clip applied to counts before normalization.
    add_counts : int
        Pseudocounts added to each bootstrap draw.
    fpr_strong : float
        Tail quantile for strong hits (e.g. 0.005 → 0.5% FPR).
    fpr_weak : float
        Tail quantile for weak hits (e.g. 0.02 → 2% FPR).
    min_per_bucket : int
        Minimum counts in any single ON/OFF bucket (filtering threshold).
    min_per_replicate_total : int
        Minimum total ON+OFF counts per replicate (filtering threshold).
    random_state : int or None
        Random seed for reproducibility.
    case_name : str
        Label used in output column names for the case condition.
    control_name : str
        Label used in output column names for the control condition.

    Returns
    -------
    results : DataFrame
        One row per domain (shared filtered index). Contains:
        - '{case_name} R1/R2/Avg/SD': rho scores for case
        - '{control_name} R1/R2/Avg/SD': rho scores for control
        - 'post_prob': posterior probability of dominance (0–1)
        - 'logit_dom': logit-transformed post_prob
        - 'sig_up', 'sig_down': strong hit calls (bool)
        - 'sig_up_weak', 'sig_down_weak': weak hit calls (bool)
    delta_df : DataFrame
        Bootstrap delta distributions, shape (n_domains, n_reps * n_boot).
        Index matches results.
    thresholds : dict
        Keys: 'thr_up', 'thr_down', 'thr_up_weak', 'thr_down_weak'.
    """
    # --- 1. Filter each condition independently ---
    case_filt = filter_counts(
        case, min_per_bucket=min_per_bucket,
        min_per_replicate_total=min_per_replicate_total
    )
    ctrl_filt = filter_counts(
        control, min_per_bucket=min_per_bucket,
        min_per_replicate_total=min_per_replicate_total
    )

    # Inner join: keep only domains present in both
    shared_idx = case_filt.index.intersection(ctrl_filt.index)
    case_filt = case_filt.loc[shared_idx]
    ctrl_filt = ctrl_filt.loc[shared_idx]

    count_cols = ["counts OFF R1", "counts ON R1", "counts OFF R2", "counts ON R2"]

    # --- 2. Bootstrap rhos ---
    # Derive two independent seeds so case and control get different draws
    rng_master = np.random.default_rng(random_state)
    seed_case = int(rng_master.integers(0, 2**31))
    seed_ctrl = int(rng_master.integers(0, 2**31))

    boot_case = bootstrap_rhos(
        case_filt[count_cols], n_boot=n_boot, thresh=thresh,
        add_counts=add_counts, random_state=seed_case
    )
    boot_ctrl = bootstrap_rhos(
        ctrl_filt[count_cols], n_boot=n_boot, thresh=thresh,
        add_counts=add_counts, random_state=seed_ctrl
    )

    # --- 3. Delta: case - control, shape (n_domains, n_reps * n_boot) ---
    delta = boot_case - boot_ctrl
    delta = np.asarray(delta, dtype=float)

    # --- 4. Build null distribution ---
    null_mask = case_filt[null_col].values.astype(bool)
    null_dist = delta[null_mask].ravel()
    null_dist.sort()
    n_null = null_dist.size

    # --- 5. Posterior probability of dominance ---
    # For each domain: fraction of null values < delta_b, averaged over bootstrap
    idx = np.searchsorted(null_dist, delta, side="right")
    post_prob = idx.mean(axis=1) / n_null

    # --- 6. Thresholds from null posterior probs ---
    neg_post_prob = post_prob[null_mask]
    thr_down      = float(np.quantile(neg_post_prob, fpr_strong))
    thr_down_weak = float(np.quantile(neg_post_prob, fpr_weak))
    thr_up        = float(np.quantile(neg_post_prob, 1 - fpr_strong))
    thr_up_weak   = float(np.quantile(neg_post_prob, 1 - fpr_weak))

    thresholds = {
        "thr_up": thr_up,
        "thr_down": thr_down,
        "thr_up_weak": thr_up_weak,
        "thr_down_weak": thr_down_weak,
    }

    # --- 7. Build results DataFrame ---
    eps = 1e-6
    logit_dom = logit(np.clip(post_prob, eps, 1 - eps))

    case_rhos = compute_rhos(case_filt[count_cols], thresh=thresh)
    ctrl_rhos = compute_rhos(ctrl_filt[count_cols], thresh=thresh)

    results = pd.DataFrame(index=shared_idx)
    results[f"{case_name} R1"]  = case_rhos["rho_R1"].values
    results[f"{case_name} R2"]  = case_rhos["rho_R2"].values
    results[f"{case_name} Avg"] = case_rhos["rho_mean"].values
    results[f"{case_name} SD"]  = case_rhos["rho_sd"].values

    results[f"{control_name} R1"]  = ctrl_rhos["rho_R1"].values
    results[f"{control_name} R2"]  = ctrl_rhos["rho_R2"].values
    results[f"{control_name} Avg"] = ctrl_rhos["rho_mean"].values
    results[f"{control_name} SD"]  = ctrl_rhos["rho_sd"].values

    results["post_prob"]    = post_prob
    results["logit_dom"]    = logit_dom

    results["sig_up"]        = post_prob > thr_up
    results["sig_up_weak"]   = (post_prob > thr_up_weak) & (post_prob <= thr_up)
    results["sig_down"]      = post_prob < thr_down
    results["sig_down_weak"] = (post_prob < thr_down_weak) & (post_prob >= thr_down)

    # Carry over the null column for downstream use
    results[null_col] = case_filt[null_col].values

    # --- Build delta_df ---
    delta_df = pd.DataFrame(delta, index=shared_idx)

    return results, delta_df, thresholds
