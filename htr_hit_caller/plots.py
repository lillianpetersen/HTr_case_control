"""
Visualization utilities for HT-recruit hit calling results.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.special import logit


def plot_replicates(
    df: pd.DataFrame,
    case_name: str = "case",
    ax=None,
) -> plt.Figure:
    """
    R1 vs R2 scatter plot with Pearson correlation annotation.

    Parameters
    ----------
    df : DataFrame with columns '{case_name} R1' and '{case_name} R2'.
    case_name : str
        Condition label used in column names.
    ax : matplotlib Axes or None
        Axes to plot on. If None, a new figure is created.
    """
    r1_col = f"{case_name} R1"
    r2_col = f"{case_name} R2"

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = ax.get_figure()

    x = df[r1_col].dropna()
    y = df[r2_col].loc[x.index]
    r, _ = stats.pearsonr(x, y)

    ax.scatter(x, y, alpha=0.5, s=10)
    lim = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lim, lim, "k--", linewidth=1)
    ax.set_xlabel("Replicate 1 (log2 ON:OFF)")
    ax.set_ylabel("Replicate 2 (log2 ON:OFF)")
    ax.set_title(f"{case_name}  Pearson r = {r:.2f}")
    plt.tight_layout()
    return fig


def plot_bootstrap_histogram(
    results: pd.DataFrame,
    delta_df: pd.DataFrame,
    thresholds: Dict[str, float],
    case_name: str = "case",
    null_col: str = "is_null",
    ax=None,
) -> plt.Figure:
    """
    Histogram of logit posterior probability with threshold lines and hit counts.

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    delta_df : DataFrame returned by call_hits (unused here, kept for API consistency).
    thresholds : dict returned by call_hits.
    case_name : str
        Label for the case condition.
    null_col : str
        Column in results indicating null domains (bool).
    ax : matplotlib Axes or None
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
    else:
        fig = ax.get_figure()

    thr_down      = thresholds["thr_down"]
    thr_down_weak = thresholds["thr_down_weak"]
    thr_up        = thresholds["thr_up"]
    thr_up_weak   = thresholds["thr_up_weak"]

    ax.hist(results["logit_dom"], bins=80, color="steelblue", alpha=0.7)
    ax.axvline(0, color="k", linewidth=1)
    ax.axvline(logit(thr_down),      linestyle="--", color="k", linewidth=1.2)
    ax.axvline(logit(thr_up),        linestyle="--", color="k", linewidth=1.2)
    ax.axvline(logit(thr_down_weak), linestyle="--", color="grey", linewidth=1)
    ax.axvline(logit(thr_up_weak),   linestyle="--", color="grey", linewidth=1)

    n_up   = results["sig_up"].sum()
    n_down = results["sig_down"].sum()
    n_up_w = results["sig_up_weak"].sum()
    n_dn_w = results["sig_down_weak"].sum()

    ymax = ax.get_ylim()[1]
    ax.text(logit(thr_down), ymax * 0.98,
            f"99.5% CI\n{n_down} down", ha="right", va="top", fontsize=9)
    ax.text(logit(thr_up), ymax * 0.98,
            f"99.5% CI\n{n_up} up", ha="left", va="top", fontsize=9)
    ax.text(logit(thr_down_weak), ymax * 0.88,
            f"98% CI\n{n_dn_w} weak", ha="right", va="top", fontsize=9, color="grey")
    ax.text(logit(thr_up_weak), ymax * 0.88,
            f"98% CI\n{n_up_w} weak", ha="left", va="top", fontsize=9, color="grey")

    ax.set_xlabel(f"({case_name} − control) logit posterior probability")
    ax.set_ylabel("Number of domains")
    ax.set_title(f"{case_name}: bootstrapping hit calls")
    plt.tight_layout()
    return fig


def plot_scatter(
    results: pd.DataFrame,
    case_name: str = "case",
    control_name: str = "control",
    ax=None,
) -> plt.Figure:
    """
    Case vs control rho scatter with error bars; highlights hits.

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    case_name : str
    control_name : str
    ax : matplotlib Axes or None
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    x     = results[f"{control_name} Avg"]
    y     = results[f"{case_name} Avg"]
    xerr  = results[f"{control_name} SD"]
    yerr  = results[f"{case_name} SD"]

    # All domains (grey)
    ax.errorbar(x, y, xerr=xerr, yerr=yerr,
                fmt="o", alpha=0.3, markersize=3, linewidth=0.5,
                color="grey", label="all")

    # Strong hits
    sig = results["sig_up"] | results["sig_down"]
    if sig.any():
        ax.errorbar(x[sig], y[sig], xerr=xerr[sig], yerr=yerr[sig],
                    fmt="o", alpha=0.9, markersize=4, linewidth=0.8,
                    color="#33a02c", label="strong hit")

    # Weak hits
    sig_w = results["sig_up_weak"] | results["sig_down_weak"]
    if sig_w.any():
        ax.errorbar(x[sig_w], y[sig_w], xerr=xerr[sig_w], yerr=yerr[sig_w],
                    fmt="o", alpha=0.9, markersize=4, linewidth=0.8,
                    color="#b2df8a", label="weak hit")

    lim = [min(x.min(), y.min()) - 0.2, max(x.max(), y.max()) + 0.2]
    ax.plot(lim, lim, "k--", linewidth=1)
    ax.set_xlabel(f"{control_name} Avg rho")
    ax.set_ylabel(f"{case_name} Avg rho")
    ax.set_title(f"{control_name} vs {case_name}")
    ax.legend(loc="best", fontsize=8)
    plt.tight_layout()
    return fig


def plot_bootstrap_examples(
    results: pd.DataFrame,
    delta_df: pd.DataFrame,
    thresholds: Dict[str, float],
    null_col: str = "is_null",
    case_name: str = "case",
    control_name: str = "control",
) -> plt.Figure:
    """
    Histogram of bootstrap delta distributions for one representative domain
    per hit category, with the pooled null distribution shown in the
    background.  Short dotted vertical ticks mark the pseudocount-adjusted
    observed delta for each replicate (mean of that replicate's bootstrap
    draws, which equals the expected value of the pseudocount-adjusted
    estimator).

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    delta_df : DataFrame returned by call_hits.
    thresholds : dict returned by call_hits.
    null_col : str
    case_name : str  (unused, kept for API consistency)
    control_name : str  (unused, kept for API consistency)
    """
    null_mask_s = results[null_col]

    # --- Pick one representative domain per category ---
    candidates = {
        "strong up":   results.loc[results["sig_up"]        & ~null_mask_s],
        "strong down": results.loc[results["sig_down"]      & ~null_mask_s],
        "weak up":     results.loc[results["sig_up_weak"]   & ~null_mask_s],
        "weak down":   results.loc[results["sig_down_weak"] & ~null_mask_s],
        "null":        results.loc[null_mask_s],
    }
    examples = {}
    for cat, sub in candidates.items():
        if len(sub) == 0:
            continue
        if cat == "strong up":
            examples[cat] = sub["post_prob"].idxmax()
        elif cat == "strong down":
            examples[cat] = sub["post_prob"].idxmin()
        elif cat == "weak up":
            examples[cat] = sub["post_prob"].idxmax()
        elif cat == "weak down":
            examples[cat] = sub["post_prob"].idxmin()
        else:  # null: most typical (closest to 0.5)
            examples[cat] = (sub["post_prob"] - 0.5).abs().idxmin()

    # --- Pooled null delta distribution (background) ---
    null_deltas = delta_df.values[null_mask_s.values].ravel()

    colors = {
        "strong up": "#1a5fa8", "weak up":     "#74c2e1",
        "strong down": "#c0392b", "weak down": "#f0a070",
        "null": "#444444",
    }

    # Shared bins across all distributions
    all_vals = np.concatenate([delta_df.loc[d].values for d in examples.values()])
    bins = np.linspace(all_vals.min(), all_vals.max(), 60)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax_null = ax.twinx()

    # Pooled null on secondary y-axis
    ax_null.hist(null_deltas, bins=bins, histtype="stepfilled", color="grey", alpha=0.2,
                 edgecolor="#333333", linewidth=0.8, label="null (pooled)")
    ax_null.set_ylabel("Null pooled count", color="grey")
    ax_null.tick_params(axis="y", labelcolor="grey")

    # delta_df columns are ordered [R1_boot0..R1_bootN, R2_boot0..R2_bootN]
    n_boot_per_rep = delta_df.shape[1] // 2

    label_info = []  # collect (center_x, peak_y, domain, pp, color) for labelling later
    for cat, domain in examples.items():
        deltas = delta_df.loc[domain].values
        counts_hist, _ = np.histogram(deltas, bins=bins)
        ax.hist(deltas, bins=bins, histtype="stepfilled",
                color=colors[cat], alpha=0.4, lw=1.5, label=f"{cat}: {domain}")
        # Mark pseudocount-adjusted observed value per replicate
        obs_r1 = deltas[:n_boot_per_rep].mean()
        obs_r2 = deltas[n_boot_per_rep:].mean()
        for obs in (obs_r1, obs_r2):
            ax.axvline(obs, color=colors[cat], lw=1.5, linestyle=":", ymax=0.12)
        # Store label info: x = bootstrap mean, y = histogram peak
        center_x = deltas.mean()
        peak_y   = counts_hist.max()
        pp       = results.loc[domain, "post_prob"]
        label_info.append((center_x, peak_y, domain, pp, colors[cat]))

    # Draw labels above each histogram peak
    for center_x, peak_y, domain, pp, color in label_info:
        ax.text(center_x, peak_y, f"{domain}\n{pp:.3f}",
                color=color, ha="center", va="bottom", fontsize=6.5,
                fontweight="bold")

    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Bootstrap delta (case − control rho)")
    ax.set_ylabel("Count (per domain)")
    ax.set_title("Bootstrap delta distributions: example domains\n"
                 "(dotted line = observed value)")

    # Combined legend from both axes
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax_null.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=7)
    plt.tight_layout()
    return fig


def plot_postprob_calibration(
    results: pd.DataFrame,
    null_col: str = "is_null",
) -> plt.Figure:
    """
    Two-panel calibration plot.

    Left: post_prob distribution for null domains — should be approximately
    uniform on [0, 1], confirming that the null-based thresholding is
    well-calibrated.

    Right: post_prob distribution for activating (non-null) domains — should
    be bimodal with peaks near 0 and 1 if hits are present.

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    null_col : str
    """
    null_pp = results.loc[results[null_col],  "post_prob"]
    act_pp  = results.loc[~results[null_col], "post_prob"]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    axes[0].hist(null_pp, bins=25, color="#444444", alpha=0.75, edgecolor="none")
    axes[0].set_xlabel("Posterior probability")
    axes[0].set_ylabel("Count")
    axes[0].set_title(f"Null domains (n={len(null_pp)})\nExpected: uniform on [0, 1]")

    axes[1].hist(act_pp, bins=25, color="#1a5fa8", alpha=0.75, edgecolor="none")
    axes[1].set_xlabel("Posterior probability")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Activating domains (n={len(act_pp)})\nHits → peaks near 0 and 1")

    plt.tight_layout()
    return fig


def plot_bootstrap_convergence(
    results: pd.DataFrame,
    delta_df: pd.DataFrame,
    null_col: str = "is_null",
    example_domains: Optional[list] = None,
    ns: Optional[list] = None,
    n_reps: int = 2,
) -> plt.Figure:
    """
    Show how posterior probability stabilises as n_boot increases, for a
    handful of example domains.

    Uses the first ``n_reps * n`` columns of ``delta_df`` to simulate
    n_boot = n bootstrap draws, then recomputes post_prob.

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    delta_df : DataFrame returned by call_hits.
    null_col : str
    example_domains : list of index labels, or None (auto-selected).
    ns : list of ints (n_boot values to evaluate), or None (sensible default).
    n_reps : int
        Number of replicates used when building delta_df (default 2).
    """
    max_boot = delta_df.shape[1] // n_reps
    if ns is None:
        ns = [n for n in [10, 25, 50, 100, 250, 500, 1000, 2000] if n <= max_boot]

    null_mask_s = results[null_col]
    null_mask   = null_mask_s.values.astype(bool)

    # Auto-select examples: strongest up, strongest down, most-typical null
    if example_domains is None:
        example_domains = []
        up_sub = results.loc[results["sig_up"] & ~null_mask_s]
        if len(up_sub):
            example_domains.append(up_sub["post_prob"].idxmax())
        dn_sub = results.loc[results["sig_down"] & ~null_mask_s]
        if len(dn_sub):
            example_domains.append(dn_sub["post_prob"].idxmin())
        null_sub = results.loc[null_mask_s]
        if len(null_sub):
            example_domains.append((null_sub["post_prob"] - 0.5).abs().idxmin())

    domain_idxs = [results.index.get_loc(d) for d in example_domains]

    _colors = ["#1a5fa8", "#c0392b", "#444444", "#74c2e1", "#f0a070"]

    fig, ax = plt.subplots(figsize=(7, 4))
    for d_idx, domain in zip(domain_idxs, example_domains):
        pp_vals = []
        for n in ns:
            n_cols    = n_reps * n
            sub       = delta_df.values[:, :n_cols]
            null_dist = np.sort(sub[null_mask].ravel())
            n_null    = len(null_dist)
            d         = sub[d_idx]
            pp        = float(np.searchsorted(null_dist, d, side="right").mean() / n_null)
            pp_vals.append(pp)
        color = _colors[example_domains.index(domain) % len(_colors)]
        ax.plot(ns, pp_vals, marker="o", markersize=4, lw=1.5, color=color, label=domain)

    ax.axhline(0.5, color="grey", lw=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("n_boot")
    ax.set_ylabel("Posterior probability")
    ax.set_title("Bootstrap convergence: post_prob vs n_boot")
    ax.legend(fontsize=7)
    plt.tight_layout()
    return fig


def make_all_plots(
    results: pd.DataFrame,
    delta_df: pd.DataFrame,
    thresholds: Dict[str, float],
    output_dir: str,
    case_name: str = "case",
    control_name: str = "control",
    null_col: str = "is_null",
) -> None:
    """
    Save all standard plots as PDFs into ``output_dir``.

    Files created:
    - replicates_{case_name}.pdf
    - bootstrap_histogram_{case_name}.pdf
    - scatter_{case_name}_vs_{control_name}.pdf
    - bootstrap_examples_{case_name}.pdf
    - postprob_calibration_{case_name}.pdf
    - bootstrap_convergence_{case_name}.pdf

    Parameters
    ----------
    results : DataFrame returned by call_hits.
    delta_df : DataFrame returned by call_hits.
    thresholds : dict returned by call_hits.
    output_dir : str
        Directory where PDFs are saved (created if it doesn't exist).
    case_name : str
    control_name : str
    null_col : str
    """
    os.makedirs(output_dir, exist_ok=True)

    fig = plot_replicates(results, case_name=case_name)
    fig.savefig(os.path.join(output_dir, f"replicates_{case_name}.pdf"), transparent=True)
    plt.close(fig)

    fig = plot_bootstrap_histogram(
        results, delta_df, thresholds, case_name=case_name, null_col=null_col
    )
    fig.savefig(
        os.path.join(output_dir, f"bootstrap_histogram_{case_name}.pdf"), transparent=True
    )
    plt.close(fig)

    fig = plot_scatter(results, case_name=case_name, control_name=control_name)
    fig.savefig(
        os.path.join(output_dir, f"scatter_{case_name}_vs_{control_name}.pdf"),
        transparent=True,
    )
    plt.close(fig)

    fig = plot_bootstrap_examples(results, delta_df, thresholds, null_col=null_col,
                                  case_name=case_name, control_name=control_name)
    fig.savefig(
        os.path.join(output_dir, f"bootstrap_examples_{case_name}.pdf"), transparent=True
    )
    plt.close(fig)

