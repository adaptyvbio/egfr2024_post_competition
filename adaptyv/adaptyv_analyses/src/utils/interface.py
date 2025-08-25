import numpy as np
from Bio.PDB import PDBParser, Selection
from scipy.spatial import cKDTree
from loguru import logger
from src.utils.constants import THREE_TO_ONE_MAP


def hotspot_residues(trajectory_pdb, binder_chain="A", target_chain="B", atom_distance_cutoff=4.0):
    """Analyze interface residues with detailed logging"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", trajectory_pdb)
        
        # Check if chains exist
        if binder_chain not in structure[0]:
            raise ValueError(f"Binder chain '{binder_chain}' not found in structure")
        if target_chain not in structure[0]:
            raise ValueError(f"Target chain '{target_chain}' not found in structure")

        binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
        if not binder_atoms:
            raise ValueError(f"No atoms found in binder chain {binder_chain}")
            
        target_atoms = Selection.unfold_entities(structure[0][target_chain], 'A')
        if not target_atoms:
            raise ValueError(f"No atoms found in target chain {target_chain}")

        binder_coords = np.array([atom.coord for atom in binder_atoms])
        target_coords = np.array([atom.coord for atom in target_atoms])

        binder_tree = cKDTree(binder_coords)
        target_tree = cKDTree(target_coords)

        interacting_residues = {}
        pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

        logger.debug(f"Processing {len(pairs)} potential interactions")
        
        for binder_idx, close_indices in enumerate(pairs):
            binder_residue = binder_atoms[binder_idx].get_parent()
            binder_resname = binder_residue.get_resname()
            if binder_resname in THREE_TO_ONE_MAP:
                aa_single_letter = THREE_TO_ONE_MAP[binder_resname]
                if close_indices:  
                    logger.debug(f"Binder residue {binder_residue.id[1]} ({binder_resname}) has {len(close_indices)} contacts")
                    interacting_residues[binder_residue.id[1]] = aa_single_letter

        logger.debug(f"Found {len(interacting_residues)} interface residues: {sorted(interacting_residues.keys())}")
        return interacting_residues

    except Exception as e:
        logger.error(f"Error in hotspot_residues: {str(e)}")
        raise  # Re-raise the exception with proper context

def analyze_interface_contacts(pdb_file, binder_chain="A", target_chain="B", atom_distance_cutoff=4.0):
    """Analyze interface contacts using KDTree"""
    try:
        logger.info(f"Analyzing interface contacts for {pdb_file}")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", pdb_file)

        binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
        binder_coords = np.array([atom.coord for atom in binder_atoms])

        target_atoms = Selection.unfold_entities(structure[0][target_chain], 'A')
        target_coords = np.array([atom.coord for atom in target_atoms])

        binder_tree = cKDTree(binder_coords)
        target_tree = cKDTree(target_coords)

        pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)
        
        binder_residues = {}
        target_residues = {}
        target_residue_indices = set()
        
        atom_contacts = sum(len(close_indices) for close_indices in pairs)
        
        for binder_idx, close_indices in enumerate(pairs):
            if close_indices:
                binder_residue = binder_atoms[binder_idx].get_parent()
                binder_resname = binder_residue.get_resname()
                
                if binder_resname in THREE_TO_ONE_MAP:
                    aa_single_letter = THREE_TO_ONE_MAP[binder_resname]
                    binder_residues[binder_residue.id[1]] = aa_single_letter
                    
                    for close_idx in close_indices:
                        target_residue = target_atoms[close_idx].get_parent()
                        target_resname = target_residue.get_resname()
                        if target_resname in THREE_TO_ONE_MAP:
                            target_aa = THREE_TO_ONE_MAP[target_resname]
                            target_residues[target_residue.id[1]] = target_aa
                            target_residue_indices.add(target_residue.id[1])

        target_residue_list = sorted(list(target_residue_indices))
        binder_residue_list = sorted(list(binder_residues.keys()))
        total_interface_residues = len(binder_residues) + len(target_residues)
        
        return {
            'binder_residues': binder_residues,
            'target_residues': target_residues,
            'target_residue_list': target_residue_list,
            'binder_residue_list': binder_residue_list,
            'total_interface_residues': total_interface_residues,
            'atom_contacts': atom_contacts
        }
        
    except Exception as e:
        logger.error(f"Error in analyze_interface_contacts: {str(e)}", exc_info=True)
        return None 