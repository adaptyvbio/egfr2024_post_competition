from Bio.PDB import PDBParser, DSSP
from collections import defaultdict
from src.utils.interface import hotspot_residues
from loguru import logger
import os


def calc_ss_percentage(pdb_file, dssp_path, binder_chain="A", target_chain="B", atom_distance_cutoff=4.0):
    try:
        # Check if files exist
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        if not os.path.exists(dssp_path):
            raise FileNotFoundError(f"DSSP executable not found: {dssp_path}")
            
        logger.info(f"Processing PDB file: {pdb_file}")
        logger.info(f"Using DSSP executable at: {dssp_path}")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        model = structure[0]

        # Log available chains
        available_chains = [c.id for c in model]
        logger.info(f"Available chains in structure: {available_chains}")

        if binder_chain not in model:
            raise ValueError(f"Chain {binder_chain} not found in structure. Available chains: {available_chains}")

        try:
            dssp = DSSP(model, pdb_file, dssp=dssp_path)
        except Exception as e:
            logger.error(f"DSSP calculation failed: {str(e)}")
            raise

        ss_counts = defaultdict(int)
        ss_interface_counts = defaultdict(int)
        plddts_interface = []
        plddts_ss = []

        chain = model[binder_chain]
        
        try:
           
            interacting_residues = set(hotspot_residues(
                pdb_file, 
                binder_chain=binder_chain,
                target_chain=target_chain,
                atom_distance_cutoff=atom_distance_cutoff
            ).keys())
            
            logger.info(f"Found {len(interacting_residues)} interacting residues")
            
            if not interacting_residues:
                logger.warning("No interacting residues found")
        except Exception as e:
            logger.error(f"Failed to calculate hotspot residues: {str(e)}")
            raise

        for residue in chain:
            try:
                residue_id = residue.id[1]
                if (binder_chain, residue_id) in dssp:
                    ss = dssp[(binder_chain, residue_id)][2]
                    ss_type = 'loop'
                    if ss in ['H', 'G', 'I']:
                        ss_type = 'helix'
                    elif ss == 'E':
                        ss_type = 'sheet'

                    ss_counts[ss_type] += 1

                    if ss_type != 'loop':
                        avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                        plddts_ss.append(avg_plddt_ss)

                    if residue_id in interacting_residues:
                        ss_interface_counts[ss_type] += 1
                        avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                        plddts_interface.append(avg_plddt_residue)
            except Exception as e:
                logger.error(f"Error processing residue {residue_id}: {str(e)}")
                continue

        total_residues = sum(ss_counts.values())
        total_interface_residues = sum(ss_interface_counts.values())

        if total_residues == 0:
            raise ValueError("No residues were processed successfully")

        helix_pct = round((ss_counts['helix'] / total_residues) * 100, 2) if total_residues > 0 else 0
        sheet_pct = round((ss_counts['sheet'] / total_residues) * 100, 2) if total_residues > 0 else 0
        loop_pct = round(((total_residues - ss_counts['helix'] - ss_counts['sheet']) / total_residues) * 100, 2) if total_residues > 0 else 0

        i_helix_pct = round((ss_interface_counts['helix'] / total_interface_residues) * 100, 2) if total_interface_residues > 0 else 0
        i_sheet_pct = round((ss_interface_counts['sheet'] / total_interface_residues) * 100, 2) if total_interface_residues > 0 else 0
        i_loop_pct = round(((total_interface_residues - ss_interface_counts['helix'] - ss_interface_counts['sheet']) / total_interface_residues) * 100, 2) if total_interface_residues > 0 else 0

        i_plddt = round(sum(plddts_interface) / len(plddts_interface) / 100, 2) if plddts_interface else 0
        ss_plddt = round(sum(plddts_ss) / len(plddts_ss) / 100, 2) if plddts_ss else 0

        return (helix_pct, sheet_pct, loop_pct, 
                i_helix_pct, i_sheet_pct, i_loop_pct, 
                i_plddt, ss_plddt)

    except Exception as e:
        logger.error(f"Error in calc_ss_percentage: {str(e)}")
        raise