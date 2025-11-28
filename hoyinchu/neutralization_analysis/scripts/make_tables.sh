## Path to your Foldseek executable
FOLDSEEK_PATH="/data1/lareauc/users/chuh/softwares/Foldseek/foldseek/bin/foldseek"
## Absolute path to this directory
PROJECT_DIRECTORY="/data1/lareauc/users/chuh/ProteinDesign/AdaptyvBio/egfr_neutralisation_analysis"

## Compute the binder tm-score to EGF
python ${PROJECT_DIRECTORY}/scripts/run_foldseek.py \
  --pdb_dir ${PROJECT_DIRECTORY}/data/structure_predictions \
  --reference_pdb ${PROJECT_DIRECTORY}/data/extra_predictions/Human_EGF_chain_A_only.pdb \
  --out ${PROJECT_DIRECTORY}/data/tables/tmscores_to_egf.csv \
  --foldseek_path ${FOLDSEEK_PATH}

## Compute the binder tm-score to round 1 winner
python ${PROJECT_DIRECTORY}/scripts/run_foldseek.py \
  --pdb_dir ${PROJECT_DIRECTORY}/data/structure_predictions \
  --reference_pdb ${PROJECT_DIRECTORY}/data/extra_predictions/round1_winner.pdb \
  --out ${PROJECT_DIRECTORY}/data/tables/tmscores_to_round1_winner.csv \
  --foldseek_path ${FOLDSEEK_PATH}

## Compute the interface residues
python ${PROJECT_DIRECTORY}/scripts/identify_interface_residues.py \
  --pdb_dir ${PROJECT_DIRECTORY}/data/structure_predictions \
  --distance_cutoff 4 \
  --output_csv ${PROJECT_DIRECTORY}/data/tables/interchain_residue_indices.csv



