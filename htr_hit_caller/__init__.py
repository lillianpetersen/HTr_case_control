"""
htr_hit_caller: Bootstrapping-based hit calling for HT-recruit kinase inhibitor screens.

Public API
----------
call_hits(case, control, null_col, ...)
    Full pipeline: filter → bootstrap → delta → posterior → hit calls.

compute_rhos(df, thresh=5)
    Log2(ON/OFF) enrichment per replicate.

bootstrap_rhos(counts_df, n_boot=1000, ...)
    Multinomial bootstrap of rho values.

filter_counts(df, min_per_bucket=250, min_per_replicate_total=700)
    Drop domains with insufficient counts.

make_all_plots(results, delta_df, thresholds, output_dir, ...)
    Save standard diagnostic plots as PDFs.
"""

from .hit_calling import call_hits
from .enrichment import compute_rhos, bootstrap_rhos
from .filtering import filter_counts
from .plots import (
    plot_replicates,
    plot_bootstrap_histogram,
    plot_scatter,
    plot_bootstrap_examples,
    plot_postprob_calibration,
    plot_bootstrap_convergence,
    make_all_plots,
)

__all__ = [
    "call_hits",
    "compute_rhos",
    "bootstrap_rhos",
    "filter_counts",
    "plot_replicates",
    "plot_bootstrap_histogram",
    "plot_scatter",
    "plot_bootstrap_examples",
    "plot_postprob_calibration",
    "plot_bootstrap_convergence",
    "make_all_plots",
]

__version__ = "0.1.0"
