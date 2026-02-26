"""
Enrichment calculations: compute_rhos and bootstrap_rhos.
"""

import numpy as np
import pandas as pd


def compute_rhos(df: pd.DataFrame, thresh: float = 5) -> pd.DataFrame:
    """
    Compute log2(ON/OFF) enrichment (rho) per replicate.

    Parameters
    ----------
    df : DataFrame with columns:
        'counts OFF R1', 'counts ON R1', 'counts OFF R2', 'counts ON R2'
    thresh : float
        Floor clip applied to all count columns before normalization.

    Returns
    -------
    DataFrame with columns: rho_R1, rho_R2, rho_mean, rho_sd
    (indexed same as input df)
    """
    df2 = df.copy()

    for col in ["counts OFF R1", "counts ON R1", "counts OFF R2", "counts ON R2"]:
        df2[col] = df2[col].clip(lower=thresh)

    total_ON_R1  = df2["counts ON R1"].sum()
    total_OFF_R1 = df2["counts OFF R1"].sum()
    total_ON_R2  = df2["counts ON R2"].sum()
    total_OFF_R2 = df2["counts OFF R2"].sum()

    df2["rho_R1"] = np.log2(
        (df2["counts ON R1"] / total_ON_R1) / (df2["counts OFF R1"] / total_OFF_R1)
    )
    df2["rho_R2"] = np.log2(
        (df2["counts ON R2"] / total_ON_R2) / (df2["counts OFF R2"] / total_OFF_R2)
    )

    df2["rho_mean"] = df2[["rho_R1", "rho_R2"]].mean(axis=1)
    df2["rho_sd"]   = df2[["rho_R1", "rho_R2"]].std(axis=1)

    return df2[["rho_R1", "rho_R2", "rho_mean", "rho_sd"]]


def bootstrap_rhos(
    counts_df: pd.DataFrame,
    n_boot: int = 1000,
    thresh: float = 5,
    add_counts: int = 100,
    random_state=None,
) -> np.ndarray:
    """
    Bootstrap log2(ON/OFF) rhos for each domain using a multinomial model.

    For each replicate, draws n_boot multinomial samples from the observed
    count distribution, adds pseudocounts, and computes rho for each draw.
    The two replicates are flattened so the output has n_reps * n_boot columns.

    Parameters
    ----------
    counts_df : DataFrame with columns:
        'counts OFF R1', 'counts ON R1', 'counts OFF R2', 'counts ON R2'
    n_boot : int
        Number of bootstrap iterations per replicate.
    thresh : float
        Floor clip applied before computing multinomial probabilities.
    add_counts : int
        Pseudocounts added to each bootstrap draw before computing rho.
        Stabilizes log-ratio estimates for low-count domains.
    random_state : int or None
        Seed for reproducibility.

    Returns
    -------
    boot_rho : ndarray of shape (n_domains, n_reps * n_boot)
        Bootstrapped rho values. Columns are [R1_boot0, R1_boot1, ..., R2_boot0, ...].
    """
    rng = np.random.default_rng(random_state)

    counts = counts_df.copy()
    for c in ["counts OFF R1", "counts ON R1", "counts OFF R2", "counts ON R2"]:
        counts[c] = counts[c].clip(lower=thresh)

    n_domains = len(counts)
    n_rep = 2
    boot_rho = np.zeros((n_domains, n_rep, n_boot))

    for r, (on_col, off_col) in enumerate(
        [("counts ON R1", "counts OFF R1"), ("counts ON R2", "counts OFF R2")]
    ):
        on_counts  = counts[on_col].values
        off_counts = counts[off_col].values

        total_on  = on_counts.sum()
        total_off = off_counts.sum()

        p_on  = on_counts  / total_on
        p_off = off_counts / total_off

        for b in range(n_boot):
            on_b  = rng.multinomial(total_on,  p_on)  + add_counts
            off_b = rng.multinomial(total_off, p_off) + add_counts

            rho_b = np.log2((on_b / total_on) / (off_b / total_off))
            boot_rho[:, r, b] = rho_b

    # Flatten replicates: shape (n_domains, n_rep * n_boot)
    return boot_rho.reshape(n_domains, -1)
