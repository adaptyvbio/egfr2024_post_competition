import subprocess
import pandas as pd
import numpy as np
from plot import _compute_correlations, _compare_distributions, _compute_roc_auc, _compute_pr_auc, _as_array, _nanfinite_mask


base_chai_cols = [
    'plddt',
    'ptm',
    'aggregate_score',
]

chai_path = 'data/chai_out.csv'
chai_df = pd.read_csv(chai_path)

all_columns = chai_df.columns.tolist()
ipsae_columns = [f'chai_{col}' for col in all_columns if 'ipsae' in col.lower()]
dock_columns = [f'chai_{col}' for col in all_columns if 'dock' in col.lower()]
lis_columns = [f'chai_{col}' for col in all_columns if 'lis' in col.lower()]
iptm_columns = ['chai_iptm', 'chai_iptm_d0chn_ab_asym', 'chai_iptm_d0chn_ba_asym', 'chai_iptm_d0chn_ab_max']
dist_columns = [f'chai_{col}' for col in all_columns if 'dist' in col.lower()]
nres_columns = [f'chai_{col}' for col in all_columns if 'nres' in col.lower()]
base_chai_cols = [f'chai_{col}' for col in base_chai_cols]

final_cols = [
    'esm_pll_avg',
    'chai_ipsae_d0chn_ab_max',
    'chai_iptm_d0chn_ab_asym',
    'chai_nres2_ba_asym',
    'chai_lis_ab_asym',
    'chai_dist2_ba_asym',
    'chai_pdockq2_ba_asym',
    'chai_pdockq_ab_asym',
    'chai_plddt',
    'chai_ptm',
    'similarity_check',
    'esm_pll',
]

column_dict = {
    'base_chai_cols': base_chai_cols,
    'ipsae_cols': ipsae_columns,
    'dock_cols': dock_columns,
    'lis_cols': lis_columns,
    'iptm_cols': iptm_columns,
    'dist_cols': dist_columns,
    'nres_cols': nres_columns,
    'prodigy_cols': [
        'prodigy_ppkd_af2',
        'prodigy_ppkd_chai',
    ],
    'synthyra_cols': [
        'synteract3_human',
        'synteract3_multi',
        'ppi-probability',
        'predicted-pKd'
    ],
    'adaptyv_cols': [
        'pae_interaction',
        'esm_pll',
        'esm_pll_avg',
        'iptm',
        'plddt',
        'similarity_check',
    ]
}


DATA_DIR = 'data'
PLOT_DIR = 'plots'


def compute_metrics_for_column(df, true_col, metric_col):
    """Compute all metrics for a single metric column."""
    y_true = _as_array(df[true_col])
    metric_vals = _as_array(df[metric_col])
    
    # Apply finite mask
    finite_mask = _nanfinite_mask(y_true, metric_vals)
    y_true = y_true[finite_mask]
    metric_vals = metric_vals[finite_mask]
    
    # Split into zero vs nonzero
    zero_mask = (y_true == 0.0)
    nonzero_mask = ~zero_mask
    
    # Compute correlations (only for nonzero values)
    y_true_nz = y_true[nonzero_mask]
    metric_nz = metric_vals[nonzero_mask]
    
    if len(y_true_nz) > 0:
        corr = _compute_correlations(y_true_nz, metric_nz)
    else:
        corr = {
            "pearson_r": float("nan"),
            "pearson_p": float("nan"),
            "spearman_r": float("nan"),
            "spearman_p": float("nan"),
        }
    
    # Compute distribution comparison (Mann-Whitney U)
    metric_zero = metric_vals[zero_mask]
    metric_nonzero = metric_vals[nonzero_mask]
    dist_stats = _compare_distributions(metric_zero, metric_nonzero)
    
    # Compute ROC AUC and PR AUC
    y_binary = np.where(zero_mask, 0, 1).astype(int)
    roc_auc = _compute_roc_auc(y_binary, metric_vals)
    pr_auc = _compute_pr_auc(y_binary, metric_vals)
    
    return {
        'metric': metric_col,
        'pearson_r': corr['pearson_r'],
        'pearson_p': corr['pearson_p'],
        'spearman_r': corr['spearman_r'],
        'spearman_p': corr['spearman_p'],
        'effect_r': dist_stats['effect_r'],
        'effect_p': dist_stats['p_value'],
        'roc_auc': roc_auc,
        'pr_auc': pr_auc,
    }


if __name__ == "__main__":
    # Load the data
    df = pd.read_csv(f'{DATA_DIR}/final.csv')
    true_col = 'pkd'
    
    # Collect all metrics
    all_metrics = []
    
    for group_name, cols in column_dict.items():
        for col in cols:
            if col in df.columns:
                metrics = compute_metrics_for_column(df, true_col, col)
                metrics['group'] = group_name
                all_metrics.append(metrics)
            else:
                print(f"Warning: Column '{col}' not found in dataframe, skipping...")
    
    # Create DataFrame and save
    metrics_df = pd.DataFrame(all_metrics)
    # Reorder columns
    metrics_df = metrics_df[['group', 'metric', 'pearson_r', 'pearson_p', 'spearman_r', 'spearman_p', 
                             'effect_r', 'effect_p', 'roc_auc', 'pr_auc']]
    metrics_df.to_csv(f'{DATA_DIR}/metrics_summary.csv', index=False)
    print(f"Metrics saved to {DATA_DIR}/metrics_summary.csv")
    
    # Generate plots
    for key, cols in column_dict.items():
        subprocess.run([
            'py',
            'plot.py',
            '--csv',
            f'{DATA_DIR}/final.csv',
            '--true_col',
            'pkd',
            '--metric_cols',
            *cols,
            '--output',
            f'{PLOT_DIR}/{key}.png'
        ])
    
    subprocess.run([
        'py',
        'plot.py',
        '--csv',
        f'{DATA_DIR}/final.csv',
        '--true_col',
        'pkd',
        '--metric_cols',
        *final_cols,
        '--output',
        f'{PLOT_DIR}/final.png',
        '--final',
    ])