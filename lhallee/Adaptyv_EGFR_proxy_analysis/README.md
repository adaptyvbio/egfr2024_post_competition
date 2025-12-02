# Adaptyv Bio 2024 EGFR Competition Binding Proxy Analysis

Computational evaluation of protein binder design metrics against experimental binding affinity data from the [Adaptyv Bio EGFR Binder Design Competitions](https://proteinbase.com/competitions).

---

## Overview

This repo evaluates computational metrics for predicting protein–protein binding affinity. We analyze predictions from multiple structure prediction and sequence-based methods against experimentally measured binding affinities (pKd) from the Adaptyv Bio EGFR competition.

The analysis includes:
- **Correlation analysis** between computational metrics and experimentally measured pKd values (for binders)
- **Distribution comparison** between binders and non-binders using Mann-Whitney U tests
- **Classification performance** via ROC-AUC and PR-AUC for distinguishing binders from non-binders

---

## Dataset

Data is sourced from the [Adaptyv Bio EGFR Binder Design Competitions](https://proteinbase.com/competitions/adaptyv-egfr-binder):
- **Round 1**: 200 computationally designed protein binders ([Competition Link](https://proteinbase.com/competitions/adaptyv-egfr-binder))
- **Round 2**: 400 additional designs ([Competition Link](https://proteinbase.com/competitions/adaptyv-egfr-binder2))

The target is the extracellular domain of human EGFR (Epidermal Growth Factor Receptor), a critical cancer therapy target.

---

## Metrics Evaluated

### Structure Prediction (Chai-1)

We chose to leverage Chai-1 for structure prediction because of its excellent performance without MSAs, which is valuable for _de novo_ sequences

Calculated from structure predictions with standard hyperparameters (3 recycle steps, 200 diffusion steps, 5 diffusion samples)
- The sample for each complex with the highest Chai-1 aggregate score is utilized for the rest of the analysis
- ipSAE variants, iPTM variants, pDockQ, pDockQ2, LIS, etc. are calculated via the script [here](https://github.com/DunbrackLab/IPSAE/blob/main/ipsae.py)

| Category | Metrics |
|----------|---------|
| **Confidence** | pLDDT, pTM, aggregate score |
| **Interface pTM (iPTM)** | Standard iPTM, d0chn/d0dom variants, asymmetric/max aggregations |
| **Interface pSAE (ipSAE)** | d0chn/d0dom variants with asymmetric/max aggregations |
| **Docking Quality** | pDockQ, pDockQ2 (asymmetric/max) |
| **Local Interaction Score** | LIS (asymmetric/max) |
| **Interface Residues** | nRes (residues with PAE<10), Dist (PAE<10 and distance<10Å) |

### Structure-Based Binding Affinity
- **PRODIGY**: Predicted pKd from AlphaFold2 and Chai-1 structures

### Competition Baseline Metrics
- PAE interaction score
- ESM pseudo-log-likelihood (ESM-PLL)
- Sequence similarity check

### Sequence-Based PPI Prediction
- **Synteract2**: PPI probability and predicted pKd
- **Synteract3**: Human and multi-species PPI probability models

---

## Repository Structure

```
BinderAnalysis/
├── data/                      # Processed datasets and metrics
│   ├── final.csv              # Merged dataset with all metrics
│   ├── metrics_summary.csv    # Summary statistics for all metrics
│   └── ...
├── plots/                     # Generated visualization figures
├── af2_preds/                 # AlphaFold2 structure predictions (PDB)
├── chai1_preds/               # Chai-1 structure predictions (CIF)
├── prep_csv.py                # Data preparation and merging
├── plot.py                    # Core plotting and statistical analysis
├── plot_plots.py              # Batch plot generation
├── run_prodigy.py             # PRODIGY binding affinity predictions
├── synteract3.py              # Synteract3 PPI predictions
├── synthyra_api/              # Synthyra API utilities
├── prodigy_local/             # Local PRODIGY implementation
├── metrics/                   # Regression and classification utilities
└── requirements.txt           # Python dependencies
```

---

## Reproducing Plots

The original competition metrics, PRODIGY, Synteract2, and Synteract3 metrics have been precomputed. You can easily redo the calculation from the PRODIGY model via `python run_prodigy.py`

To regenerate all analysis plots:

```bash
# 1. Prepare and merge all metrics into final.csv
python prep_csv.py

# 2. Generate all comparison plots
python plot_plots.py
```

This will produce:
- Individual metric group plots in `plots/` (e.g., `ipsae_cols.png`, `dock_cols.png`)
- Combined summary plot (`plots/final.png`)
- Metrics summary table (`data/metrics_summary.csv`)

---

## Statistical Methods

Each metric is evaluated via:

1. **Correlation (binders only)**: Pearson and Spearman correlation with experimental pKd
2. **Distribution separation**: Mann-Whitney U test comparing metric distributions for binders vs. non-binders
3. **Classification**: ROC-AUC and PR-AUC for binary binder/non-binder classification
4. **Effect size**: Rank-biserial correlation (r) indicating practical significance

Plots are sorted by effect size magnitude for easy comparison across metrics.

---

## Citations

Work utilized in the this repo

### Competition preprint
```bibtex
```

### AlphaFold2
```bibtex
@article{jumper2021alphafold,
  title={Highly accurate protein structure prediction with AlphaFold},
  author={Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek, Augustin and Potapenko, Anna and others},
  journal={Nature},
  volume={596},
  number={7873},
  pages={583--589},
  year={2021},
  doi={10.1038/s41586-021-03819-2}
}
```

### Chai-1
```bibtex
@article{chai2024,
  title={Chai-1: Decoding molecular interactions in proteins, nucleic acids, and molecules},
  author={{Chai Discovery team}},
  journal={bioRxiv},
  year={2024},
  doi={10.1101/2024.10.10.615955}
}
```

### ipSAE
```bibtex
@article{chakravarty2025ipsae,
  title={ipSAE: predicted inter-chain solvent-accessible surface area energy for evaluation of AlphaFold-Multimer and Chai-1 models of protein complexes},
  author={Chakravarty, Devlina and Porter, Lauren L and Dunbrack, Roland L},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.02.10.637595}
}
```

### pDockQ
```bibtex
@article{bryant2022ppi,
  title={Improved prediction of protein-protein interactions using AlphaFold2},
  author={Bryant, Patrick and Pozzati, Gabriele and Elofsson, Arne},
  journal={Nature Communications},
  volume={13},
  number={1},
  pages={1265},
  year={2022},
  doi={10.1038/s41467-022-28865-w}
}
```

### pDockQ2
```bibtex
@article{zhu2023pdockq2,
  title={Evaluation of AlphaFold-Multimer prediction on multi-chain protein complexes},
  author={Zhu, Wensi and Shenoy, Aditi and Kundrotas, Petras and Elofsson, Arne},
  journal={Bioinformatics},
  volume={39},
  number={7},
  pages={btad424},
  year={2023},
  doi={10.1093/bioinformatics/btad424}
}
```

### LIS (Local Interaction Score)
```bibtex
@article{kim2024lis,
  title={Benchmarking AlphaFold2 on peptide binding structure prediction},
  author={Kim, Junsu and Chung, Woojin and Jang, Junho and Kim, Jinsol and Kim, Kwangho and Seok, Chaok},
  journal={bioRxiv},
  year={2024},
  doi={10.1101/2024.02.19.580970}
}
```

### Synteract2
```bibtex
@artciple{hallee2025diffusionsequencemodelsenhanced,
    title={Diffusion Sequence Models for Enhanced Protein Representation and Generation}, 
    author={Logan Hallee and Nikolaos Rafailidis and David B. Bichara and Jason P. Gleghorn},
    journal={arXiv},
    year={2025},
    doi={10.48550/arXiv.2506.08293}, 
}
```

---

## Links

- **Competition Round 1**: [Adaptyv EGFR Binder Design Competition 1](https://proteinbase.com/competitions/adaptyv-egfr-binder)
- **Competition Round 2**: [Adaptyv EGFR Binder Design Competition 2](https://proteinbase.com/competitions/adaptyv-egfr-binder2)
- **Dataset**: [Synthyra/AdaptyvBioRound2EGFR](https://huggingface.co/datasets/Synthyra/AdaptyvBioRound2EGFR)

---

## License

This project is provided for research purposes. Please refer to the individual tool licenses (PRODIGY, Chai-1, etc.) for their respective terms.
