import pandas as pd
import ast
import glob

from loguru import logger
import os
import shutil


# Round 1: submissions file and affinity file
round_1_path = "./data/round_1/submissions/submissions_annotated.csv"
round_1_affinity = './data/round_1/data_package/result_summary.csv'
submissions_round_1 = pd.read_csv(round_1_path).drop_duplicates(subset='sequence')
submissions_round_1 = submissions_round_1.rename(columns={'model_names': 'model_category', 'category': 'design_category'})
affinity_round_1 = pd.read_csv(round_1_affinity)
affinity_round_1_missing_columns = pd.read_csv('./data/round_1/data_package/result_summary_missing_columns.csv')
affinity_round_1 = pd.merge(affinity_round_1, affinity_round_1_missing_columns, on = 'name', how = 'left')


remaining_metrics = pd.read_csv('./data/round_1/data_package/remaining_metrics.csv')
remaining_metrics = remaining_metrics.rename(columns={'name': 'id'})

remaining_pll = pd.read_csv('./data/round_1/data_package/remaining_pll.csv')[['id', 'esm_pll']]
remaining_metrics = pd.merge(remaining_metrics, remaining_pll, on='id', how='left')

# Merge round 1 submissions with affinity data
submissions_round_1 = pd.merge(
    submissions_round_1,
    affinity_round_1,
    on='sequence',
    how='left'
)


# Round 2: submissions file and affinity file
round_2_path = "./data/round_2/submissions/submissions_annotated.csv"
round_2_affinity = './data/round_2/data_package/result_summary.csv'
submissions_round_2 = pd.read_csv(round_2_path).drop_duplicates(subset='sequence')
submissions_round_2 = submissions_round_2.rename(columns={'main_model_category': 'model_category', 'category': 'design_category'})
affinity_round_2 = pd.read_csv(round_2_affinity)
submissions_round_2 = pd.merge(
    submissions_round_2,
    affinity_round_2,
    on='sequence',
    how='left'
)
submissions_round_2 = submissions_round_2.rename(columns = {'username_x': 'username'})

# Pooled submissions with ids
submissions_pooled = pd.read_csv('./data/raw/submissions/submissions_annotated.csv')
# Create a mapping of sequence to id from submissions_pooled
sequence_to_id = submissions_pooled.set_index('sequence')['id'].to_dict()

# Update ids in submissions_round_1 based on sequence matches
submissions_round_1['id'] = submissions_round_1['sequence'].map(sequence_to_id)
submissions_round_2['id'] = submissions_round_2['sequence'].map(sequence_to_id)

# Convert string representation of dict to actual dict
submissions_pooled['scores'] = submissions_pooled['scores'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

# Get all unique score keys across all dictionaries
score_keys = set()
for scores in submissions_pooled['scores'].dropna():
    score_keys.update(scores.keys())

for key in ['esm_pll', 'pae_interaction', 'plddt', 'iptm', 'ptm']:
    submissions_pooled[key] = submissions_pooled['scores'].apply(lambda x: x.get(key) if isinstance(x, dict) else None)
    submissions_round_1[key] = submissions_round_1['id'].map(submissions_pooled.set_index('id')[key])
    submissions_round_2[key] = submissions_round_2['id'].map(submissions_pooled.set_index('id')[key])

# Rank the round 1 submissions to match ranking of the second round
submissions_round_1 = submissions_round_1.sort_values('pae_interaction')
submissions_round_1['top_100'] = submissions_round_1['pae_interaction'].rank() <= 100

submissions_round_1 = submissions_round_1.drop('pae_interaction', axis=1)
submissions_round_1 = submissions_round_1.drop('plddt', axis=1)
submissions_round_1 = submissions_round_1.drop('iptm', axis=1)
submissions_round_1 = submissions_round_1.drop('ptm', axis=1)
submissions_round_1 = submissions_round_1.drop('esm_pll', axis=1)


submissions_round_1 = pd.merge(
    submissions_round_1,
    remaining_metrics,
    on='id',
    how='left'
)



# Fix boolean indexing syntax
submissions_round_1.loc[(submissions_round_1['selected'] == 'Yes') & (submissions_round_1['top_100'] == True), 'selected'] = 'Top 100'
submissions_round_1.loc[(submissions_round_1['selected'] == 'Yes') & (submissions_round_1['top_100'] == False), 'selected'] = 'Adaptyv selection'

# Fix round 2 selection update
submissions_round_2.loc[submissions_round_2['selected'] == 'Adaptyv Selection', 'selected'] = 'Adaptyv selection'


submissions_round_1['round'] = 1
submissions_round_2['round'] = 2
# Get shared columns between round 1 and 2 submissions
shared_cols = list(set(submissions_round_1.columns) & set(submissions_round_2.columns))


# Select only shared columns and concatenate
submissions = pd.concat([
    submissions_round_1[shared_cols],
    submissions_round_2[shared_cols]
], ignore_index=True)


# Standardize design categories by removing redundant labels and + filter suffix
submissions['design_category'] = submissions['design_category'].replace({
    'diversified binder + filter': 'diversified binder',
    'de novo + filter': 'de novo',
    'de novo + filter ': 'de novo',  # Extra space variant
    'hallucination + filter': 'hallucination'
})

submissions['de_novo'] = submissions['id'].map(
     submissions_pooled.set_index('id')['is_de_novo']
).map({True: 'De novo', False: 'Existing binder'})

submissions = submissions.drop('similarity_check', axis=1)
homology_data = pd.read_csv('./data/raw/homology/homology_data.csv')
submissions = pd.merge(submissions, homology_data, on='id', how='left')
submissions['sequence_length'] = submissions['sequence'].apply(len)
# Remove socials column if it exists
# Define desired column order - ensure 'round' is included
ordered_columns = [
    'id', 'round', 'sequence', 'selected', 'kd', 'similarity_check',
    'design_process', 'design_category', 'model_category', 'sequence_length'
]

# Remove socials and name if they exist in submissions columns
if 'socials' in submissions.columns:
    submissions = submissions.drop('socials', axis=1)

if 'dna' in submissions.columns:
    submissions = submissions.drop('dna', axis=1)

#if 'name' in submissions.columns:
#    submissions = submissions.drop('name', axis=1)

#if 'username' in submissions.columns:
#    submissions = submissions.drop('username', axis=1)

# Add remaining columns that aren't in the ordered list
remaining_cols = [col for col in submissions.columns if col not in ordered_columns]
ordered_columns.extend(remaining_cols)

submissions = submissions[ordered_columns]

# Some more renaming for consistency
submissions['binding'] = submissions['binding'].replace({
    'FALSE': 'No',
    'false': 'No',
    'False': 'No',
    'TRUE': 'Yes',
    'true': 'Yes',
    'True': 'Yes',
    'unknown': 'Unknown'
})

for col in ['expression', 'binding_strength', 'model_category', 'design_category']:
    if col in submissions.columns:
        submissions[col] = submissions[col].str.capitalize()

submissions['normalized_esm_pll'] = submissions['esm_pll'] / submissions['sequence_length']


# Read and concatenate all competition metrics CSVs
metrics_files = glob.glob('./data/processed/competition_metrics_csvs/*.csv')

metrics_dfs = []
for file in metrics_files:
    try:
        df = pd.read_csv(file)
        metrics_dfs.append(df)
    except Exception as e:
        continue

if not metrics_dfs:
    raise RuntimeError("No metrics files were successfully loaded")

metrics_combined = pd.concat(metrics_dfs, ignore_index=True)

# Merge with submissions and drop specified columns
submissions = pd.merge(submissions, metrics_combined, on='id', how='left')
columns_to_drop = ['aa_counts', 'interface_residues', 'status']
submissions = submissions.drop(columns=columns_to_drop, errors='ignore')

# Replace None values with NaN in expression and binding_strength columns
for col in ['expression', 'binding_strength']:
    if col in submissions.columns:
        submissions[col] = submissions[col].replace({'None': pd.NA, None: pd.NA})


# Read and merge Wells Wood TSV data
wells_wood_df = pd.read_csv('./data/wells_wood/foldseek_adaptyv_destress_binder_merged.tsv', sep='\t')

# Prefix all columns except sequence with wells_wood_
rename_cols = {col: f'wells_wood_{col}' for col in wells_wood_df.columns if col != 'sequence'}
wells_wood_df = wells_wood_df.rename(columns=rename_cols)

# Merge with submissions on sequence
submissions = pd.merge(submissions, wells_wood_df, on='sequence', how='left')


submissions.to_csv('./data/processed/all_submissions.csv', index=False)
submissions.to_parquet('./data/processed/all_submissions.parquet', compression='zstd', compression_level=3)

# Print submissions columns for debugging
print(submissions.columns)

# ####
# # Processing structures
# data_path = './data/raw/structures'
# structure_dirs = ['001_2024', '002_2024']
# structure_paths = [f"./data/raw/structures/{d}" for d in structure_dirs]
#
# from loguru import logger
# import os
# import shutil
#
#
#
# def save_pdbs_per_round(structure_dirs, submissions):
#     """Save PDB files per round, only including structures with two chains"""
#     # Track counts per round
#     round_counts = {1: 0, 2: 0}
#
#     # Track which files we've already copied to avoid duplicates
#     copied_files = set()
#
#     # Track missing files
#     missing_files = []
#
#     from Bio.PDB import PDBParser, Selection
#
#     for directory in structure_dirs:
#         for filename in os.listdir(directory):
#             file_id = filename.split('_')[0]
#             if file_id in submissions['id'].values and "rank_001" in filename:
#                 # Check if structure has two chains
#                 parser = PDBParser(QUIET=True)
#                 try:
#                     structure = parser.get_structure("complex", os.path.join(directory, filename))
#                     chains = Selection.unfold_entities(structure[0], 'C')
#
#                     if len(chains) != 2:
#                         logger.warning(f"Skipping {filename} - found {len(chains)} chains instead of 2")
#                         continue
#
#                     round = submissions[submissions['id'] == file_id]['round'].values[0]
#                     output_dir = f"./data/round_{round}/structures/alphafold"
#                     output_path = f"{output_dir}/{file_id}.pdb"
#
#                     # Only copy if we haven't already copied this file_id
#                     if file_id not in copied_files:
#                         os.makedirs(output_dir, exist_ok=True)
#                         shutil.copy(os.path.join(directory, filename), output_path)
#                         round_counts[round] += 1
#                         copied_files.add(file_id)
#                 except Exception as e:
#                     logger.error(f"Error processing {filename}: {str(e)}")
#                     continue
#
#     # Log counts and check for missing files
#     logger.info(f"Saved {round_counts[1]} two-chain PDB files for round 1")
#     logger.info(f"Saved {round_counts[2]} two-chain PDB files for round 2")
#
#     # Check for missing files in each round
#     for round_num in [1, 2]:
#         expected_count = submissions[submissions['round'] == round_num].shape[0]
#         if round_counts[round_num] != expected_count:
#             logger.warning(f"Expected {expected_count} PDB files for round {round_num}, but found {round_counts[round_num]} with two chains")
#             # Find which files are missing
#             round_submissions = submissions[submissions['round'] == round_num]
#             for _, row in round_submissions.iterrows():
#                 submission_id = row['id']
#                 if submission_id not in copied_files:
#                     missing_files.append({
#                         'id': submission_id,
#                         'round': round_num
#                     })
#                     logger.warning(f"Missing or invalid PDB file for submission {submission_id} in round {round_num}")
#
#     # Count actual files in output directories
#     for round_num in [1, 2]:
#         output_dir = f"./data/round_{round_num}/structures/alphafold"
#         if os.path.exists(output_dir):
#             actual_count = len([f for f in os.listdir(output_dir) if f.endswith('.pdb')])
#             logger.info(f"Found {actual_count} two-chain PDB files in {output_dir}")
#             if actual_count != round_counts[round_num]:
#                 logger.warning(f"Mismatch in round {round_num}: copied {round_counts[round_num]} but found {actual_count}")
#                 # Check which files are present in directory but not in copied_files
#                 for filename in os.listdir(output_dir):
#                     if filename.endswith('.pdb'):
#                         file_id = filename.split('.')[0]
#                         if file_id not in copied_files:
#                             logger.warning(f"Found unexpected file {filename} in {output_dir}")
#
#
# structure_dirs = ['./data/raw/structures/001_2024_remaining']
# save_pdbs_per_round(structure_dirs, submissions)
