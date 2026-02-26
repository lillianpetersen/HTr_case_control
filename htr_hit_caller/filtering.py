"""
Count filtering for HT-recruit data.
"""

import pandas as pd


_COUNT_COLS = ["counts OFF R1", "counts ON R1", "counts OFF R2", "counts ON R2"]
_R1_COLS = ["counts OFF R1", "counts ON R1"]
_R2_COLS = ["counts OFF R2", "counts ON R2"]


def filter_counts(
    df: pd.DataFrame,
    min_per_bucket: int = 250,
    min_per_replicate_total: int = 700,
) -> pd.DataFrame:
    """
    Filter domains by minimum count thresholds.

    Two criteria are applied:
    1. Any individual bucket (ON or OFF, either replicate) below ``min_per_bucket``
       → domain is dropped.
    2. Sum of ON+OFF for either replicate below ``min_per_replicate_total``
       → domain is dropped.

    Parameters
    ----------
    df : DataFrame with columns:
        'counts OFF R1', 'counts ON R1', 'counts OFF R2', 'counts ON R2'
    min_per_bucket : int
        Minimum counts in any single ON or OFF channel.
    min_per_replicate_total : int
        Minimum total ON+OFF counts per replicate.

    Returns
    -------
    Filtered copy of df (rows that fail either filter are dropped).
    """
    df2 = df.copy()

    # Drop domains where total per replicate is too low
    r1_total = df2[_R1_COLS].sum(axis=1)
    r2_total = df2[_R2_COLS].sum(axis=1)
    low_total = (r1_total < min_per_replicate_total) | (r2_total < min_per_replicate_total)
    df2 = df2.loc[~low_total]

    # Drop domains where any single bucket is too low
    low_bucket = (df2[_COUNT_COLS] < min_per_bucket).any(axis=1)
    df2 = df2.loc[~low_bucket]

    return df2
