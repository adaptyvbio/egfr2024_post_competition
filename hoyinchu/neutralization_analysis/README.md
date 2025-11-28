# AdaptyvBio EGFR Binder Competition Round 2 Neutralisation Analysis

This repo contains the follow-up analysis of the neutralisation data from the Round 2 AdaptyvBio EGFR Binder Design Competition (descirbed in the preprint [here](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2)). Most of the analysis is contained in the `EGFR_neutralisation_analysis.ipynb` notebook

![Neutralisation Coefficient Ranking](./figures/neutralisation_coefficient_ranking.pdf)

## Data Sources & Description:
`data/extract_predictions`
1. `Human_EGF.pdb`: This prediction pdb is obtained by running colabfold notebook with all default settings using the identical EGF and EGFR sequences reported in the competition. The prediction meta data can be found in the zip file `EGF_EGFR_61de7.result.zip` and the rank_001 model is used.
2. `Human_EGF_chain_A/B_only.pdb`: Same as `Human_EGF.pdb` but only chain A (EGF) or B (EGFR) is preserved and reindexed.

`data/tables`
1. `foldseek_adaptyv_destress_binder_merged.tsv`: Destress & Other Biophysical Metrics of the binders obtained from https://github.com/wells-wood-research/adaptyv-bio-analysis/blob/main/data/foldseek_adaptyv_destress_binder_merged.tsv
2. `interchain_residue_indices.csv`: Set of residues within 4 angstrom of chain A and chain B
3. `neutralisation_coefficient.csv`: Manually extracted neutralisation coefficient
4. `tmscores_to_egf.csv`: TM-score of the binders to EGF
5. `neutralisation_df_filtered`: A merged dataframe of other tables in this directory filtered for binders tested with neutralisation assay. Constructed by the notebook

`data/structure_predictions`: This should be downloaded from [here](https://api.adaptyvbio.com/storage/v1/object/public/egfr_design_competition_2/structure_predictions.zip)

`data/pymol_bfactor_pdbs`: PDBs with bfactors modified for visualisation purposes

## Figures

All figures are produced from the notebook `EGFR_neutralisation_analysis.ipynb` with the exception of the EGFR structure which was produced using ChimeraX. 


## Reproducing the analysis
Pre-requisite: [Foldseek](https://github.com/steineggerlab/foldseek) (to compute binder TM-score against Human EGF)

1. Download and unzip the structure prediction zip files from https://github.com/adaptyvbio/egfr_competition_2 and place the output directory under `data/`
2. Copy the `Human_EGF.pdb` file to the unzipped `structure_predictions` directory
3. Run `bash make_tables.sh` under the `scripts` directory and provide a path to the FoldSeek executable
4. Now you have all the data required for running the Jupyter notebook `EGFR_neutralisation_analysis.ipynb`!


Related Repos:
- https://github.com/adaptyvbio/egfr_competition_1
- https://github.com/adaptyvbio/egfr_competition_2

## References

- Tudor-Stefan Cotet, Igor Krawczuk, Filippo Stocco, Noelia Ferruz, Anthony Gitter, Yoichi Kurumida, Lucas de Almeida Machado, Francesco Paesani, Cianna N. Calia, Chance A. Challacombe, Nikhil Haas, Ahmad Qamar, Bruno E. Correia, Martin Pacesa, Lennart Nickel, Kartic Subr, Leonardo V. Castorina, Maxwell J. Campbell, Constance Ferragu, Patrick Kidger, Logan Hallee, Christopher W. Wood, Michael J. Stam, Tadas Kluonis, Suleyman Mert Unal, Elian Belot, Alexander Naka, Adaptyv Competition Organizers. Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition. bioRxiv (2025) doi: https://doi.org/10.1101/2025.04.17.648362

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all.
Nature Methods (2022) doi: 10.1038/s41592-022-01488-1

- Jumper et al. "Highly accurate protein structure prediction with AlphaFold."
Nature (2021) doi: 10.1038/s41586-021-03819-2

