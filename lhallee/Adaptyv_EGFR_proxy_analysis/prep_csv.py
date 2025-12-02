import pandas as pd
from datasets import load_dataset


DATA_DIR = 'data'


def map_metric_to_dataset(base_df: pd.DataFrame, new_df: pd.DataFrame, output_csv: str):
    # Match SeqA and SeqB columns and create a copy to avoid SettingWithCopyWarning
    base_df = base_df[base_df['SeqA'].isin(new_df['SeqA']) & base_df['SeqB'].isin(new_df['SeqB'])].copy()
    # Merge on SeqA and SeqB to ensure proper alignment
    # Only merge columns that don't already exist in base_df
    cols_to_merge = [col for col in new_df.columns if col not in base_df.columns]
    if cols_to_merge:
        merge_cols = ['SeqA', 'SeqB']
        base_df = base_df.merge(new_df[merge_cols + cols_to_merge], on=merge_cols, how='inner')
    # Sort by pkd
    base_df = base_df.sort_values('pkd', ascending=False)
    # Save to new csv
    base_df.to_csv(output_csv, index=False)


base_cols = [
    'SeqA',
    'SeqB',
    'plddt',
    'ptm',
    'plddt',
    'aggregate_score',
]

chai_path = f'{DATA_DIR}/chai_out.csv'
chai_df = pd.read_csv(chai_path)

all_columns = chai_df.columns.tolist()

ipsae_columns = [col for col in all_columns if 'ipsae' in col.lower()]
dock_columns = [col for col in all_columns if 'dock' in col.lower()]
lis_columns = [col for col in all_columns if 'lis' in col.lower()]
iptm_columns = ['iptm', 'iptm_d0chn_ab_asym', 'iptm_d0chn_ba_asym', 'iptm_d0chn_ab_max']
dist_columns = [col for col in all_columns if 'dist' in col.lower()]
nres_columns = [col for col in all_columns if 'nres' in col.lower()]

print(all_columns)
print("-" * 100)
print(ipsae_columns)
print("-" * 100)
print(dock_columns)
print("-" * 100)
print(lis_columns)
print("-" * 100)
print(iptm_columns)
print("-" * 100)
print(dist_columns)
print("-" * 100)
print(nres_columns)

cols_to_keep = base_cols + ipsae_columns + dock_columns + lis_columns + iptm_columns + dist_columns + nres_columns

data = load_dataset("Synthyra/AdaptyvBioRound2EGFR", split='train')
data.to_pandas().to_csv(f"{DATA_DIR}/data.csv", index=False)

base_csv = f"{DATA_DIR}/data.csv"
new_csv = f"{DATA_DIR}/synteract2_out.csv"
output_csv = f"{DATA_DIR}/final.csv"

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)
map_metric_to_dataset(base_df, new_df, output_csv)

base_csv = f"{DATA_DIR}/final.csv"
new_csv = f"{DATA_DIR}/synteract3_out.csv"

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)
map_metric_to_dataset(base_df, new_df, output_csv)

new_csv = f"{DATA_DIR}/prodigy_ppkd_af2.csv"

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)

map_metric_to_dataset(base_df, new_df, output_csv)

new_csv = f"{DATA_DIR}/prodigy_ppkd_chai.csv"

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)

map_metric_to_dataset(base_df, new_df, output_csv)

new_csv = f'{DATA_DIR}/chai_out.csv'

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)
new_df = new_df[cols_to_keep]
# rename all cols except SeqA and SeqB to say chai_{col_name}
new_df = new_df.rename(columns={col: f'chai_{col}' for col in new_df.columns if col not in ['SeqA', 'SeqB']})

map_metric_to_dataset(base_df, new_df, output_csv)

new_csv = f'{DATA_DIR}/adaptyv_results.csv'

base_df = pd.read_csv(base_csv)
new_df = pd.read_csv(new_csv)
new_df['esm_pll_avg'] = new_df['esm_pll'] / new_df['sequence'].apply(len)

base_df = base_df[base_df['name'].isin(new_df['name'])].copy()
base_df = base_df.merge(new_df, on='name', how='inner')
base_df = base_df.sort_values('pkd', ascending=False)
base_df.to_csv(output_csv, index=False)