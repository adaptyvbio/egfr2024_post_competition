import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


COLUMN_NAME_DICT = {
    # Standard Chai1
    'chai_plddt': 'pLDDT',
    'chai_ptm': 'pTM',
    'chai_aggregate_score': 'Chai1 Aggregate Score',
    # iPTM
    'chai_iptm': 'iPTM',
    'chai_iptm_d0chn_ab_asym': 'iPTM d0chn Asym AB',
    'chai_iptm_d0chn_ba_asym': 'iPTM d0chn Asym BA',
    'chai_iptm_d0dom_ab_max': 'iPTM d0dom Max AB',
    'chai_iptm_d0dom_ba_max': 'iPTM d0dom Max BA',
    'chai_iptm_af_ab_asym': 'iPTM AF Asym AB',
    'chai_iptm_af_ba_asym': 'iPTM AF Asym BA',
    'chai_iptm_af_ab_max': 'iPTM AF Max AB',
    'chai_iptm_af_ba_max': 'iPTM AF Max BA',
    # ipSAE
    'chai_ipsae_ab_max': 'ipSAE Max AB',
    'chai_ipsae_ba_max': 'ipSAE Max BA',
    'chai_ipsae_ab_asym': 'ipSAE Asym AB',
    'chai_ipsae_ba_asym': 'ipSAE Asym BA',
    'chai_ipsae_d0chn_ab_max': r'ipSAE$_{d0chn}$ AB Max',
    'chai_ipsae_d0dom_ab_max': r'ipSAE$_{d0dom}$ AB Max',
    'chai_ipsae_d0chn_ba_max': r'ipSAE$_{d0chn}$ BA Max',
    'chai_ipsae_d0dom_ba_max': r'ipSAE$_{d0dom}$ BA Max',
    'chai_ipsae_d0chn_ab_asym': r'ipSAE$_{d0chn}$ AB Asym',
    'chai_ipsae_d0chn_ba_asym': r'ipSAE$_{d0chn}$ BA Asym',
    'chai_ipsae_d0dom_ab_asym': r'ipSAE$_{d0dom}$ AB Asym',
    'chai_ipsae_d0dom_ba_asym': r'ipSAE$_{d0dom}$ BA Asym',
    # Dock
    'chai_pdockq_ab_max': 'pDockQ AB Max',
    'chai_pdockq2_ab_max': 'pDockQ2 AB Max',
    'chai_pdockq_ba_max': 'pDockQ BA Max',
    'chai_pdockq2_ba_max': 'pDockQ2 BA Max',
    'chai_pdockq_ab_asym': 'pDockQ AB Asym',
    'chai_pdockq2_ab_asym': 'pDockQ2 AB Asym',
    'chai_pdockq_ba_asym': 'pDockQ BA Asym',
    'chai_pdockq2_ba_asym': 'pDockQ2 BA Asym',
    # LIS
    'chai_lis_ab_max': 'LIS Max AB',
    'chai_lis_ba_max': 'LIS Max BA',
    'chai_lis_ab_asym': 'LIS Asym AB',
    'chai_lis_ba_asym': 'LIS Asym BA',
    # nRes
    'chai_n0res_ab_max': 'Number of residues for d0 in ipSAE',
    'chai_n0chn_ab_max': r'Number of residues for d0 in $ipSAE_{d0chn}$',
    'chai_n0dom_ab_max': r'Number of residues for d0 in $ipSAE_{d0dom}$',
    'chai_d0res_ab_max': r'd0 for $ipSAE_{d0res}$',
    'chai_d0chn_ab_max': r'd0 for $ipSAE_{d0chn}$',
    'chai_d0dom_ab_max': r'd0 for $ipSAE_{d0dom}$',
    'chai_nres1_ab_max': 'Residues Chain 1 PAE<10',
    'chai_nres2_ab_max': 'Residues Chain 2 PAE<10',
    'chai_nres1_ab_asym': 'Residues Chain 1 PAE<10 Asym AB',
    'chai_nres1_ba_asym': 'Residues Chain 1 PAE<10 Asym BA',
    'chai_nres2_ab_asym': 'Residues Chain 2 PAE<10 Asym AB',
    'chai_nres2_ba_asym': 'Residues Chain 2 PAE<10 Asym BA',
    # Dist
    'chai_dist1_ab_max': 'Residues Chain 1 PAE<10 and Dist<10 Max AB',
    'chai_dist2_ab_max': 'Residues Chain 2 PAE<10 and Dist<10 Max AB',
    'chai_dist1_ab_asym': 'Residues Chain 1 PAE<10 and Dist<10 Asym AB',
    'chai_dist2_ab_asym': 'Residues Chain 2 PAE<10 and Dist<10 Asym AB',
    'chai_dist1_ba_asym': 'Residues Chain 1 PAE<10 and Dist<10 Asym BA',
    'chai_dist2_ba_asym': 'Residues Chain 2 PAE<10 and Dist<10 Asym BA',
    # Prodigy
    'prodigy_ppkd_af2': 'Prodigy (AF2)',
    'prodigy_ppkd_chai': 'Prodigy (Chai1)',
    # Synthyra
    'synteract3_human': r'$Synteract3_{Human}$ PPI Probability',
    'synteract3_multi': r'$Synteract3_{Multi}$ PPI Probability',
    'ppi-probability': 'Synteract2 PPI Probability',
    'predicted-pKd': 'Synteract2 Predicted pKd',
    # Adaptyv
    'pae_interaction': 'PAE Interaction',
    'esm_pll': 'ESM-PLL',
    'esm_pll_avg': 'ESM-PLL Average',
    'similarity_check': 'Similarity Check',
    'plddt': 'pLDDT',
    'iptm': 'iPTM',
    'ptm': 'pTM',
}


COLUMN_NAME_FINAL = {
    'esm_pll_avg': 'ESM-PLL Average',
    'chai_ipsae_d0chn_ab_max': 'ipSAE',
    'chai_iptm_d0chn_ab_asym': 'iPTM',
    'chai_nres2_ba_asym': 'Num Residues Chain PAE<10',
    'chai_lis_ab_asym': 'Local interaction score',
    'chai_dist2_ba_asym': 'Residues PAE<10 and Dist<10',
    'chai_pdockq2_ba_asym': 'pDockQ2',
    'chai_pdockq_ab_asym': 'pDockQ',
    'chai_plddt': 'pLDDT',
    'chai_ptm': 'pTM',
    'similarity_check': 'Similarity Check',
    'esm_pll': 'ESM-PLL',
}


Number = Union[int, float, np.number]



def _as_array(x: Union[Sequence[Number], np.ndarray, pd.Series]) -> np.ndarray:
    if isinstance(x, np.ndarray):
        return x.astype(float)
    if isinstance(x, pd.Series):
        return x.to_numpy(dtype=float)
    return np.asarray(list(x), dtype=float)


def _nanfinite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones_like(arrays[0], dtype=bool)
    for a in arrays:
        mask &= np.isfinite(a)
    return mask


def _compute_correlations(y_true_nz: np.ndarray, metric_nz: np.ndarray) -> Dict[str, float]:
    r_pearson, p_pearson = pearsonr(y_true_nz, metric_nz)
    r_spearman, p_spearman = spearmanr(y_true_nz, metric_nz)
    return {
        "pearson_r": float(r_pearson),
        "pearson_p": float(p_pearson),
        "spearman_r": float(r_spearman),
        "spearman_p": float(p_spearman),
    }


def _mw_effect_size(u_stat: float, n1: int, n2: int) -> float:
    # Rank-biserial correlation: r_rb = 1 - 2U/(n1*n2)
    # where U is the Mann-Whitney U statistic for the first group
    # This gives: r > 0 when first group has lower ranks (lower values)
    #             r < 0 when first group has higher ranks (higher values)
    if n1 <= 0 or n2 <= 0:
        return float("nan")
    # The maximum possible U is n1*n2, so U/(n1*n2) ranges from 0 to 1
    # r_rb ranges from 1 (first group all lower) to -1 (first group all higher)
    return 1.0 - 2.0 * (u_stat / float(n1 * n2))


def _compute_roc_auc(y_binary: np.ndarray, metric_vals: np.ndarray) -> float:
    """Compute ROC AUC for binary classification (0=non-binder, 1=binder)."""
    if y_binary.size == 0:
        return float("nan")
    # Check if we have both classes
    unique_classes = np.unique(y_binary)
    if len(unique_classes) < 2:
        return float("nan")
    try:
        return float(roc_auc_score(y_binary, metric_vals))
    except ValueError:
        return float("nan")


def _compute_pr_auc(y_binary: np.ndarray, metric_vals: np.ndarray) -> float:
    """Compute PR AUC (area under precision-recall curve) for binary classification (0=non-binder, 1=binder)."""
    if y_binary.size == 0:
        return float("nan")
    # Check if we have both classes
    unique_classes = np.unique(y_binary)
    if len(unique_classes) < 2:
        return float("nan")
    try:
        precision, recall, _ = precision_recall_curve(y_binary, metric_vals)
        return float(auc(recall, precision))
    except ValueError:
        return float("nan")


def _compare_distributions(metric_zero: np.ndarray, metric_nonzero: np.ndarray) -> Dict[str, float]:
    # Two-sided Mann–Whitney U
    # Compares non-binders (metric_zero) vs binders (metric_nonzero)
    if metric_zero.size == 0 or metric_nonzero.size == 0:
        return {"u_stat": float("nan"), "p_value": float("nan"), "effect_r": float("nan")}
    u_stat, p_value = mannwhitneyu(metric_zero, metric_nonzero, alternative="two-sided")
    # Compute effect size from perspective of binders (positive class)
    # If binders have higher values, effect_r should be positive
    # The formula r_rb = 1 - 2U/(n1*n2) gives r from perspective of first group (non-binders)
    # So: r > 0 when non-binders have lower values → binders have higher values ✓
    # To get from binders' perspective, we could flip sign, but current interpretation is correct
    effect_r = _mw_effect_size(u_stat, metric_zero.size, metric_nonzero.size)
    return {"u_stat": float(u_stat), "p_value": float(p_value), "effect_r": float(effect_r)}


def _format_corr_text(corr: Mapping[str, float]) -> str:
    return (
        f"Spearman r={corr['spearman_r']:.3f} (p={corr['spearman_p']:.2e})\n"
        f"Pearson r={corr['pearson_r']:.3f} (p={corr['pearson_p']:.2e})"
    )


def _format_mw_text(stats: Mapping[str, float], roc_auc: float, pr_auc: float) -> str:
    roc_str = f"ROC AUC={roc_auc:.3f}" if np.isfinite(roc_auc) else "ROC AUC=NaN"
    pr_str = f"PR AUC={pr_auc:.3f}" if np.isfinite(pr_auc) else "PR AUC=NaN"
    return f"{roc_str}, {pr_str}\np={stats['p_value']:.2e}, effect r={stats['effect_r']:.3f}"


def analyze_and_plot(true_value: Union[Sequence[Number], np.ndarray, pd.Series],
                     metrics: Union[Mapping[str, Union[Sequence[Number], np.ndarray, pd.Series]],
                                    Sequence[Union[Sequence[Number], np.ndarray, pd.Series]]],
                     metric_names: Optional[Sequence[str]] = None,
                     output_path: Optional[str] = None,
                     figsize: Optional[Tuple[float, float]] = None,
                     scatter_alpha: float = 0.7,
                     violin_inner: Optional[str] = "box",
                     point_size: float = 15.0,
                     use_final_names: bool = False) -> plt.Figure:
    """
    Build a horizontal figure with, per metric (one row per metric):
      - Left: scatter of nonzero true_value vs metric with Spearman/Pearson
      - Right: violin plots of metric for true==0 vs true!=0 with MW test

    Args:
        true_value: Iterable of floats (ground truth). Zeros define the negative group.
        metrics: Either a dict name->values or a list of metric arrays. If a list is
                 provided, metric_names must be provided or names will be auto-assigned.
        metric_names: Optional names for the metrics if a list is supplied.
        output_path: If provided, save the figure to this path.
        figsize: Optional custom figure size; defaults to a horizontal layout.
        scatter_alpha: Alpha for scatter points.
        violin_inner: Inner representation for seaborn.violinplot (e.g., "box", "quartile", None).
        point_size: Scatter point size.
        use_final_names: If True, use COLUMN_NAME_FINAL for display names; otherwise use COLUMN_NAME_DICT.

    Returns:
        The matplotlib Figure.
    """
    y_true = _as_array(true_value)

    # Normalize metrics into a dict of name -> array
    if isinstance(metrics, Mapping):
        metric_dict = {str(k): _as_array(v) for k, v in metrics.items()}
    else:
        metric_list = [ _as_array(v) for v in metrics ]  # type: ignore[arg-type]
        if metric_names is None:
            metric_names = [f"metric_{i+1}" for i in range(len(metric_list))]
        metric_dict = {str(name): arr for name, arr in zip(metric_names, metric_list)}

    # Harmonize lengths via masking
    n = y_true.shape[0]
    for name, arr in metric_dict.items():
        if arr.shape[0] != n:
            raise ValueError(f"Metric '{name}' length {arr.shape[0]} != true_value length {n}")

    # Store per-metric finite masks (each metric processed independently)
    # This ensures metrics with missing values don't affect other metrics
    finite_masks_by_metric: Dict[str, np.ndarray] = {}
    for metric_name, metric_vals in metric_dict.items():
        finite_masks_by_metric[metric_name] = _nanfinite_mask(y_true, metric_vals)

    num_metrics = len(metric_dict)
    if figsize is None:
        # width scales with columns; height scales with number of metrics
        figsize = (12.0, max(2.0, 2.0 * num_metrics))

    fig, axes = plt.subplots(num_metrics, 2, figsize=figsize, 
                             gridspec_kw={'hspace': 0.05, 'wspace': 0.15})
    if num_metrics == 1:
        axes = np.asarray([axes])  # shape (1, 2)

    # Set one title per column (no per-row titles)
    axes[0, 0].set_title("Correlation vs. binding affinity", fontsize=8, pad=3)
    axes[0, 1].set_title("Binders vs. Non-binders", fontsize=8, pad=3)

    # Precompute effect sizes and ROC AUC per metric (using per-metric masks), sort by descending |r|
    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    cmap = mpl.cm.Reds
    stats_by_metric: Dict[str, Dict[str, float]] = {}
    roc_auc_by_metric: Dict[str, float] = {}
    pr_auc_by_metric: Dict[str, float] = {}
    abs_r_by_metric: Dict[str, float] = {}
    for metric_name, metric_vals in metric_dict.items():
        # Apply per-metric finite mask
        finite_mask = finite_masks_by_metric[metric_name]
        y_true_masked = y_true[finite_mask]
        metric_vals_masked = metric_vals[finite_mask]
        
        # Zero vs nonzero split for this metric
        zero_mask = (y_true_masked == 0.0)
        nonzero_mask = ~zero_mask
        
        # Compute statistics using masked data
        stats = _compare_distributions(metric_vals_masked[zero_mask], metric_vals_masked[nonzero_mask])
        stats_by_metric[metric_name] = stats
        
        # Create binary labels for ROC AUC and PR AUC
        y_binary = np.where(zero_mask, 0, 1).astype(int)
        roc_auc = _compute_roc_auc(y_binary, metric_vals_masked)
        roc_auc_by_metric[metric_name] = roc_auc
        pr_auc = _compute_pr_auc(y_binary, metric_vals_masked)
        pr_auc_by_metric[metric_name] = pr_auc
        
        effect_r = stats.get("effect_r", float("nan"))
        abs_r_by_metric[metric_name] = float(abs(effect_r)) if np.isfinite(effect_r) else 0.0
    sorted_metric_names = sorted(metric_dict.keys(), key=lambda n: abs_r_by_metric[n], reverse=True)

    # Select the appropriate column name dictionary
    name_dict = COLUMN_NAME_FINAL if use_final_names else COLUMN_NAME_DICT

    for row_idx, metric_name in enumerate(sorted_metric_names):
        display_name = name_dict.get(metric_name, metric_name)
        metric_vals = metric_dict[metric_name]
        
        # Apply per-metric finite mask
        finite_mask = finite_masks_by_metric[metric_name]
        y_true_masked = y_true[finite_mask]
        metric_vals_masked = metric_vals[finite_mask]
        
        # Zero vs nonzero split for this metric
        zero_mask = (y_true_masked == 0.0)
        nonzero_mask = ~zero_mask
        
        # Nonzero subset for correlations
        y_true_nz = y_true_masked[nonzero_mask]
        metric_nz = metric_vals_masked[nonzero_mask]

        # Left subplot: scatter + lowess-like trend via regplot (no assumption of equality)
        ax_scatter = axes[row_idx, 0]
        sns.scatterplot(x=y_true_nz, y=metric_nz, ax=ax_scatter, s=point_size, alpha=scatter_alpha)
        sns.regplot(x=y_true_nz, y=metric_nz, scatter=False, ax=ax_scatter, line_kws={"color": "grey"})
        if len(display_name) > 25: # for long names, use smaller font size
            display_name_font_size = 4
        else:
            display_name_font_size = 7
        ax_scatter.set_ylabel(display_name, fontsize=display_name_font_size)
        ax_scatter.tick_params(labelsize=7)
        ax_scatter.margins(x=0.05, y=0.05)

        corr = _compute_correlations(y_true_nz, metric_nz)
        ax_scatter.text(
            0.02, 0.98, _format_corr_text(corr), transform=ax_scatter.transAxes,
            va="top", ha="left", bbox=dict(facecolor="white", alpha=0.8), fontsize=7
        )

        # Right subplot: violin plots comparing distributions where true==0 vs true!=0
        ax_violin = axes[row_idx, 1]
        group_labels = np.where(zero_mask, "Non-binders", "Binders")
        df_plot = pd.DataFrame({
            "group": group_labels,
            metric_name: metric_vals_masked,
        })

        # Retrieve precomputed stats and color from |r|
        stats = stats_by_metric[metric_name]
        abs_r = abs_r_by_metric[metric_name]
        color = cmap(norm(abs_r))

        # Color both violins by the effect-size category
        sns.violinplot(
            data=df_plot,
            x="group",
            y=metric_name,
            hue="group",
            ax=ax_violin,
            inner=violin_inner,
            palette={"Non-binders": color, "Binders": color},
            dodge=False,
        )
        # Remove legend if added due to hue
        leg = ax_violin.get_legend()
        if leg is not None:
            leg.remove()
        sns.stripplot(data=df_plot, x="group", y=metric_name, ax=ax_violin, color="k", alpha=0.25, size=1.5)
        ax_violin.set_xlabel("")
        ax_violin.set_ylabel("")
        ax_violin.tick_params(labelsize=7)
        ax_violin.margins(x=0.05, y=0.05)

        roc_auc = roc_auc_by_metric[metric_name]
        pr_auc = pr_auc_by_metric[metric_name]
        ax_violin.text(
            0.02, 0.98, _format_mw_text(stats, roc_auc, pr_auc), transform=ax_violin.transAxes,
            va="top", ha="left", bbox=dict(facecolor="white", alpha=0.8), fontsize=7
        )

        # Only show x-axis labels on the bottom row for each column
        if row_idx < num_metrics - 1:
            ax_scatter.set_xlabel("")
            ax_scatter.tick_params(labelbottom=False)
            ax_violin.set_xlabel("")
            ax_violin.tick_params(labelbottom=False)
        else:
            ax_scatter.set_xlabel("Binding affinity (pKd)", fontsize=8)
            ax_violin.set_xlabel("Binders, Non-binders", fontsize=8)

    for r in range(num_metrics):
        axes[r, 0].grid(True, linestyle=":", alpha=0.4)
        axes[r, 1].grid(True, linestyle=":", alpha=0.4)

    # Figure-level colorbar for |r| -> Reds colormap
    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    cmap = mpl.cm.Reds
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, location="right", fraction=0.015, pad=0.005)
    cbar.set_label("Effect size |r|", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # Tighten layout to minimize whitespace
    plt.tight_layout(pad=0.2, h_pad=0.1, w_pad=0.1)
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    return fig


def _parse_metrics_from_df(df: pd.DataFrame, metric_cols: List[str]) -> Dict[str, pd.Series]:
    metrics = {}
    for col in metric_cols:
        if col not in df.columns:
            raise KeyError(f"Metric column '{col}' not found in CSV.")
        metrics[col] = df[col]
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Analyze metrics vs. true baseline: correlation (nonzero) + distribution split.")
    parser.add_argument("--csv", type=str, help="Path to CSV file")
    parser.add_argument("--true_col", type=str, help="Column name for true baseline")
    parser.add_argument("--metric_cols", type=str, nargs="+", help="One or more metric column names")
    parser.add_argument("--output", "-o", type=str, default=None, help="Path to save the figure (PNG/PDF)")
    parser.add_argument("--figwidth", type=float, default=None, help="Optional figure width in inches")
    parser.add_argument("--rowheight", type=float, default=2.0, help="Height per metric row in inches")
    parser.add_argument("--final", action="store_true", help="Use final column names")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if args.true_col not in df.columns:
        raise KeyError(f"True column '{args.true_col}' not found in CSV.")

    metrics = _parse_metrics_from_df(df, args.metric_cols)

    if args.figwidth is None:
        figwidth = 12.0
    else:
        figwidth = args.figwidth
    figsize = (figwidth, max(2.0, args.rowheight * max(1, len(metrics))))

    analyze_and_plot(
        true_value=df[args.true_col],
        metrics=metrics,
        output_path=args.output,
        figsize=figsize,
        use_final_names=args.final,
    )

    if args.output is None:
        plt.show()


if __name__ == "__main__":
    main()


