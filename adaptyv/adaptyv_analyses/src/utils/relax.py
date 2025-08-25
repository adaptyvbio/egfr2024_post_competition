import os
from loguru import logger
import pyrosetta as pr
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from src.utils.pdb import clean_pdb


def pr_relax(pdb_file, relaxed_pdb_path):
    """Perform FastRelax on a PDB structure"""
    try:
            
        if os.path.exists(relaxed_pdb_path):
            logger.info(f"Relaxed structure {relaxed_pdb_path} already exists")
            return True
            
        logger.info(f"Starting relaxation of {pdb_file}")
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        # Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True)
        mmf.set_bb(True)
        mmf.set_jump(False)

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf)
        fastrelax.max_iter(200)
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        
        fastrelax.apply(pose)
        
        # Align and copy B factors
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)
        
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)
        
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)
        logger.info("Relaxation completed successfully")
        
        if os.path.exists(relaxed_pdb_path):
            file_size = os.path.getsize(relaxed_pdb_path)
            logger.info(f"Successfully saved relaxed structure ({file_size} bytes)")
            return True
        else:
            logger.error("Failed to save relaxed structure - file does not exist")
            return False
        
    except Exception as e:
        logger.error(f"Error in pr_relax: {str(e)}", exc_info=True)
        return False