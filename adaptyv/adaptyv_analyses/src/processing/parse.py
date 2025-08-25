import os
import glob
from loguru import logger

def get_pdb_files(structure_dirs):
    """Get dictionary of submission ID to PDB file path"""
    logger.info(f"Searching for PDB files in: {structure_dirs}")
    pdb_files = {}

    for directory in structure_dirs:
        if not os.path.exists(directory):
            logger.warning(f"Directory does not exist: {directory}")
            continue
            
        for filename in os.listdir(directory):
        
            # Extract submission ID from filename (everything before first underscore)
            submission_id = filename.split('.')[0]
            
            pdb_files[submission_id] = os.path.join(directory, filename)
            
    logger.info(f"Found {len(pdb_files)} PDB files")
    return pdb_files


def check_missing_relaxed_structures(structure_dirs, save_path, data_path):
    """Check for missing relaxed structures and available alternatives"""
    logger.info("Checking for missing relaxed structures...")

    missing_relaxed = []
    alternative_ranks = {}
    
    # Get all unrelaxed structures
    unrelaxed_pdbs = []
    for dir in structure_dirs:
        unrelaxed_pdbs.extend(glob.glob(f"{data_path}/{dir}/*unrelaxed_rank_001_*.pdb"))
    
    # Check each unrelaxed structure
    for pdb in unrelaxed_pdbs:
        base_name = os.path.basename(pdb)
        submission_id = base_name.split('_unrelaxed')[0]
        
        # Expected relaxed path
        relaxed_dir = f"{save_path}/relaxed"
        relaxed_pdb = os.path.join(relaxed_dir, base_name.replace('unrelaxed', 'relaxed'))
        
        if not os.path.exists(relaxed_pdb):
            missing_relaxed.append(submission_id)
            
            # Check for alternative ranks
            base_dir = os.path.dirname(pdb)
            alt_ranks = glob.glob(f"{base_dir}/{submission_id}_unrelaxed_rank_*_*.pdb")
            if len(alt_ranks) > 1:  # More than just rank_001
                alternative_ranks[submission_id] = [os.path.basename(p) for p in alt_ranks]
        
    # Log results
    logger.info(f"Found {len(missing_relaxed)} missing relaxed structures")
    for submission_id in missing_relaxed:
        if submission_id in alternative_ranks:
            logger.info(f"Submission {submission_id} missing relaxed structure but has alternatives: {alternative_ranks[submission_id]}")
        else:
            logger.info(f"Submission {submission_id} missing relaxed structure with no alternatives")
            
    return missing_relaxed, alternative_ranks


def check_missing_relaxed_structures(structure_dirs, save_path, data_path):
    """Check for missing relaxed structures and available alternatives"""
    logger.info("Checking for missing relaxed structures...")

    missing_relaxed = []
    alternative_ranks = {}
    
    # Get all unrelaxed structures
    unrelaxed_pdbs = []
    for dir in structure_dirs:
        unrelaxed_pdbs.extend(glob.glob(f"{data_path}/{dir}/*unrelaxed_rank_001_*.pdb"))
    
    # Check each unrelaxed structure
    for pdb in unrelaxed_pdbs:
        base_name = os.path.basename(pdb)
        submission_id = base_name.split('_unrelaxed')[0]
        
        # Expected relaxed path
        relaxed_dir = f"{save_path}/relaxed"
        relaxed_pdb = os.path.join(relaxed_dir, f"{submission_id}.pdb") # different methods: either submission_id or submission_id_unrelaxed_rank_001_001.pdb
        
        if not os.path.exists(relaxed_pdb):
            missing_relaxed.append(submission_id)
            
            # Check for alternative ranks
            base_dir = os.path.dirname(pdb)
            alt_ranks = glob.glob(f"{base_dir}/{submission_id}_unrelaxed_rank_*_*.pdb")
            if len(alt_ranks) > 1:  # More than just rank_001
                alternative_ranks[submission_id] = [os.path.basename(p) for p in alt_ranks]
        
    # Log results
    logger.info(f"Found {len(missing_relaxed)} missing relaxed structures")
    for submission_id in missing_relaxed:
        if submission_id in alternative_ranks:
            logger.info(f"Submission {submission_id} missing relaxed structure but has alternatives: {alternative_ranks[submission_id]}")
        else:
            logger.info(f"Submission {submission_id} missing relaxed structure with no alternatives")
            
    return missing_relaxed, alternative_ranks
