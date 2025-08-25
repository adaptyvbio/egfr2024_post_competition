import os
from loguru import logger
import pyrosetta as pr
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, LayerSelector
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.simple_metrics.metrics import TotalEnergyMetric, SasaMetric

from src.utils.relax import pr_relax
from src.utils.interface import hotspot_residues
from src.setup.pyrosetta import initialize_pyrosetta

import math


def score_interface(pdb_file, binder_chain="A", target_chain="B", save_path="/root/competition_vol"):
    
    try:
        pose = pr.pose_from_pdb(pdb_file)
        logger.info(f"Loaded pose with {pose.total_residue()} residues")
    except Exception as e:
        logger.error(f"Failed to load pose: {str(e)}")
        return None, None, None

    # analyze interface statistics
    logger.info("Setting up InterfaceAnalyzerMover")
    iam = InterfaceAnalyzerMover()
    iam.set_interface("A_B")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    
    logger.info("Applying InterfaceAnalyzerMover")
    iam.apply(pose)
    logger.info("Successfully applied InterfaceAnalyzerMover")

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    interface_residues_set = hotspot_residues(pdb_file, binder_chain = binder_chain, target_chain = target_chain)
    logger.info(f"Found {len(interface_residues_set)} interface residues")
    
    interface_residues_pdb_ids = []
    
    # Iterate over the interface residues
    logger.info("Processing interface residues")
    for pdb_res_num, aa_type in interface_residues_set.items():
        #logger.debug(f"Processing residue {pdb_res_num} of type {aa_type}")
        interface_AA[aa_type] += 1
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    interface_nres = len(interface_residues_pdb_ids)
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)
    logger.info(f"Processed {interface_nres} interface residues")

    # Calculate hydrophobicity
    logger.info("Calculating interface hydrophobicity")
    hydrophobic_aa = set('ACFILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100 if interface_nres != 0 else 0
    logger.info(f"Interface hydrophobicity: {interface_hydrophobicity}%")

    # retrieve statistics
    logger.info("Getting interface statistics")
    interfacescore = iam.get_all_data()
    #logger.debug(f"Raw interface score data: {interfacescore}")
    
    try:
        interface_sc = interfacescore.sc_value
        logger.info(f"Shape complementarity: {interface_sc}")
    except Exception as e:
        logger.error(f"Error getting sc_value: {str(e)}")
        interface_sc = 0
        
    try:
        interface_interface_hbonds = interfacescore.interface_hbonds
        logger.info(f"Interface H-bonds: {interface_interface_hbonds}")
    except Exception as e:
        logger.error(f"Error getting interface_hbonds: {str(e)}")
        interface_interface_hbonds = 0
        
    try:
        interface_dG = iam.get_interface_dG()
        logger.info(f"Interface dG: {interface_dG}")
    except Exception as e:
        logger.error(f"Error getting interface_dG: {str(e)}")
        interface_dG = 0
        
    try:
        interface_dSASA = iam.get_interface_delta_sasa()
        logger.info(f"Interface dSASA: {interface_dSASA}")
    except Exception as e:
        logger.error(f"Error getting interface_delta_sasa: {str(e)}")
        interface_dSASA = 0
        
    try:
        interface_packstat = iam.get_interface_packstat()
        logger.info(f"Interface packstat: {interface_packstat}")
    except Exception as e:
        logger.error(f"Error getting interface_packstat: {str(e)}")
        interface_packstat = 0
        
    try:
        interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100
        logger.info(f"Interface dG/SASA ratio: {interface_dG_SASA_ratio}")
    except Exception as e:
        logger.error(f"Error getting dG_dSASA_ratio: {str(e)}")
        interface_dG_SASA_ratio = 0

    # Calculate buried unsatisfied H-bonds
    try:
        buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
        interface_delta_unsat_hbonds = buns_filter.report_sm(pose)
        logger.info(f"Found {interface_delta_unsat_hbonds} buried unsatisfied H-bonds")
    except Exception as e:
        logger.error(f"Error calculating buried unsatisfied H-bonds: {str(e)}")
        interface_delta_unsat_hbonds = 0

    # Calculate H-bond percentages
    logger.info("Calculating H-bond percentages")
    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None
    logger.info(f"H-bond percentage: {interface_hbond_percentage}")

    # calculate binder energy score
    logger.info("Calculating binder energy score")
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)
    logger.info(f"Binder score: {binder_score}")

    # calculate binder SASA fraction
    logger.info("Calculating binder SASA")
    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)
    interface_binder_fraction = (interface_dSASA / binder_sasa) * 100 if binder_sasa > 0 else 0
    logger.info(f"Binder SASA: {binder_sasa}, Interface fraction: {interface_binder_fraction}")

    # calculate surface hydrophobicity
    logger.info("Calculating surface hydrophobicity")
    try:
        binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p 
                        for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]
        logger.info("Successfully split chains")
    except Exception as e:
        logger.error(f"Error splitting chains: {str(e)}")
        raise

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(binder_pose)
    logger.info(f"Selected {sum(1 for x in surface_res if x)} surface residues")

    exp_apol_count = 0
    total_count = 0 

    logger.info("Counting surface hydrophobic residues")
    for i in range(1, len(surface_res) + 1):
        if surface_res[i]:
            res = binder_pose.residue(i)
            if res.is_apolar() or res.name() in ['PHE', 'TRP', 'TYR']:
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count/total_count if total_count > 0 else 0
    logger.info(f"Surface hydrophobicity: {surface_hydrophobicity} ({exp_apol_count}/{total_count})")

    # Initialize results dictionary with default values
    interface_scores = {
        'binder_score': binder_score,
        'surface_hydrophobicity': surface_hydrophobicity,
        'interface_sc': interface_sc,
        'interface_packstat': interface_packstat,
        'interface_dG': interface_dG,
        'interface_dSASA': interface_dSASA,
        'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
        'interface_fraction': interface_binder_fraction,
        'interface_hydrophobicity': interface_hydrophobicity,
        'interface_nres': interface_nres,
        'interface_interface_hbonds': interface_interface_hbonds,
        'interface_hbond_percentage': interface_hbond_percentage,
        'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
        'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }


    return interface_scores, interface_AA, interface_residues_pdb_ids_str
