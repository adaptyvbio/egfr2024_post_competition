# BioLM EGFR Competition Analysis

Analysis notebooks for protein language model log-probabilities, DNA feature importance, and molecule type comparisons.

## Setup

```bash
cd biolm/
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Generate Figures

```bash
cd biolm/

# ESM log-probability analysis → plots/figure13_log_probabilities.png
jupyter nbconvert --to notebook --execute protein_llm_log_probability_plot.ipynb

# DNA feature importance → plots/figure18_feature_importance.png
jupyter nbconvert --to notebook --execute dna_feature_importance_plot.ipynb

# Molecule type ANOVA → plots/anova_molecule_types.csv
jupyter nbconvert --to notebook --execute molecule_type_anova_table.ipynb

# Molecule type distribution → plots/molecule_type_distribution.png
jupyter nbconvert --to notebook --execute molecule_type_plot.ipynb
```

Data files in `data/` are tracked with Git LFS. If not installed: `brew install git-lfs && git lfs install` (macOS) or `sudo apt install git-lfs && git lfs install` (Ubuntu).
