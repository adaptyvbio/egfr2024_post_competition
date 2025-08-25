import os
import pandas as pd
from loguru import logger
from pathlib import Path
import multiprocess as mp
import functools
from tqdm import tqdm
from src.setup.pyrosetta import initialize_pyrosetta
from src.utils.pdb import clean_pdb, add_cryst1_record
from src.utils.interface import analyze_interface_contacts
from src.utils.scores import score_interface
from src.utils.dssp import calc_ss_percentage
from src.processing.parse import get_pdb_files
from src.utils.relax import pr_relax

def save_batch_results(all_results, submissions, save_path):
    if len(all_results) > 0:
        logger.info(f"Processed {len(all_results)} submissions successfully")
        
        results_dir = f"{save_path}/results"
        os.makedirs(results_dir, exist_ok=True)
        
        output_path = os.path.join(results_dir, "final_results.csv")
        
        results_df = pd.DataFrame(all_results)
        
        merged_df = submissions.merge(results_df, on='id', how='left')
        
        merged_df.to_csv(output_path, index=False)
        logger.info(f"Saved results to {output_path}")
        return True
    else:
        logger.error("No submissions were processed successfully in batch")
        return False

def process_submission(submission_data, pdb_files, dssp_path, save_path = '/root/competition_vol', binder_chain="A", target_chain="B"):
    """Process a single submission with structure relaxation, DSSP and interface scoring"""
    try:
        submission_id = submission_data['id']
        submission_sequence = submission_data['sequence']
        
        if submission_id not in pdb_files:
            error_msg = f"No PDB file found for submission {submission_id}"
            logger.warning(error_msg)
            return submission_id, submission_sequence, None, None, None, error_msg
        
        pdb_path = pdb_files[submission_id]
        
        # First relax the structure using pr_relax
        try:
            relaxed_dir = os.path.join(save_path, "relaxed")
            os.makedirs(relaxed_dir, exist_ok=True)
            relaxed_path = os.path.join(relaxed_dir, f"{submission_id}.pdb")
            
            success = pr_relax(pdb_path, relaxed_path)
            if not success:
                raise ValueError("Structure relaxation failed")
                
            # Use relaxed structure for subsequent analysis
            pdb_path = relaxed_path
            logger.info(f"Successfully relaxed structure for {submission_id}")
            
        except Exception as e:
            error_msg = f"Structure relaxation failed for {submission_id}: {str(e)}"
            logger.error(error_msg)
            return submission_id, submission_sequence, None, None, None, error_msg

        # Get DSSP results
        try:
            dssp_results = calc_ss_percentage(pdb_path, dssp_path, binder_chain=binder_chain, target_chain=target_chain)
            if not dssp_results:
                raise ValueError("DSSP calculation returned None")
            logger.info(f"DSSP results for {submission_id}: {dssp_results}")
        except Exception as e:
            error_msg = f"DSSP calculation failed for {submission_id}: {str(e)}"
            logger.error(error_msg)
            return submission_id, submission_sequence, None, None, None, error_msg

        # Get interface results
        try:
            interface_results = score_interface(pdb_path, binder_chain=binder_chain, target_chain=target_chain)
            if not interface_results:
                raise ValueError("Interface scoring returned None")
            logger.info(f"Interface results for {submission_id}: {interface_results}")
        except Exception as e:
            error_msg = f"Interface scoring failed for {submission_id}: {str(e)}"
            logger.error(error_msg)
            return submission_id, submission_sequence, None, None, None, error_msg

        # Get contact analysis
        try:
            contact_results = analyze_interface_contacts(pdb_path, binder_chain=binder_chain, target_chain=target_chain)
            if not contact_results:
                raise ValueError("Contact analysis returned None")
            logger.info(f"Contact analysis for {submission_id}: {contact_results}")
        except Exception as e:
            error_msg = f"Contact analysis failed for {submission_id}: {str(e)}"
            logger.error(error_msg)
            return submission_id, submission_sequence, None, None, None, error_msg

        return submission_id, submission_sequence, dssp_results, interface_results, contact_results, None
        
    except Exception as e:
        error_msg = f"Error processing submission {submission_id}: {str(e)}"
        logger.error(error_msg, exc_info=True)
        return submission_id, submission_sequence, None, None, None, error_msg

def process_submissions_batch(input_csv: str, structure_dirs: list, start_idx: int, end_idx: int, save_path: str, data_path: str):
    """Process a batch of submissions"""
    logger.info(f"Starting DSSP annotation process for batch {start_idx}-{end_idx}")
    
   
    # Adjust paths for volume
    volume_input_csv = f"{data_path}/{Path(input_csv).name}"
    volume_structure_dirs = [f"{data_path}/{d}" for d in structure_dirs]
    
    # Read the submissions data
    logger.info(f"Reading submissions from {volume_input_csv}")
    try:
        submissions = pd.read_csv(volume_input_csv)
        # Filter to just the batch we want
        submissions = submissions.iloc[start_idx:end_idx].copy()
        logger.info(f"Loaded {len(submissions)} submissions for batch")
        
        if len(submissions) == 0:
            logger.error(f"No submissions found in range {start_idx}-{end_idx}")
            return False
            
        # Get PDB files
        pdb_files = get_pdb_files(volume_structure_dirs)
        logger.info(f"Found {len(pdb_files)} PDB files")
        
        # Process submissions
        dssp_path = "/root/bindcraft/functions/dssp"
        all_results = []
        
        # Create submission data list - pass Series instead of tuples
        submission_data = [row for _, row in submissions.iterrows()]
        
        # Process submissions in parallel using multiprocess Pool
        num_cores = mp.cpu_count()
        logger.info(f"Using {num_cores} CPU cores for parallel processing")
        
        # Initialize PyRosetta in each worker process
        def init_worker():
            initialize_pyrosetta()
        
        with mp.Pool(processes=num_cores, initializer=init_worker) as pool:
            # Create partial function with fixed arguments
            process_func = functools.partial(
                process_submission, 
                pdb_files=pdb_files, 
                save_path=save_path,
                dssp_path=dssp_path,
                binder_chain="A",
                target_chain="B"
            )
            
            # Process submissions in parallel with progress bar
            results = []
            with tqdm(total=len(submission_data), desc=f"Processing batch {start_idx}-{end_idx}") as pbar:
                for result in pool.imap_unordered(process_func, submission_data):
                    results.append(result)
                    pbar.update()
        
        # Process results
        for result in results:
            try:
                submission_id, sequence, dssp_results, interface_results, contact_results, error = result
                
                if error:
                    logger.error(f"Error processing submission {submission_id}: {error}")
                    logger.error(f"Sequence: {sequence}")
                    # Add failed submission with error message
                    all_results.append({
                        'id': submission_id,
                        'error': error,
                        'status': 'failed'
                    })
                    continue
                    
                if not all([dssp_results, interface_results, contact_results]):
                    error_msg = f"Missing results for submission {submission_id}"
                    logger.error(error_msg)
                    all_results.append({
                        'id': submission_id,
                        'error': error_msg,
                        'status': 'failed'
                    })
                    continue
                    
                # Unpack the interface results tuple
                interface_metrics, aa_counts, interface_residues = interface_results
                
                result_dict = {
                    'id': submission_id,
                    'helix_pct': dssp_results[0],
                    'sheet_pct': dssp_results[1], 
                    'loop_pct': dssp_results[2],
                    'i_helix_pct': dssp_results[3],
                    'i_sheet_pct': dssp_results[4],
                    'i_loop_pct': dssp_results[5],
                    'i_plddt': dssp_results[6],
                    'ss_plddt': dssp_results[7],
                    **interface_metrics,  # Unpack the metrics dictionary
                    'aa_counts': str(aa_counts),  # Convert to string to ensure serialization
                    'interface_residues': str(interface_residues),  # Convert to string to ensure serialization
                    **contact_results,
                    'status': 'success'
                }
                all_results.append(result_dict)
                logger.info(f"Successfully processed submission {submission_id}")
            except Exception as e:
                logger.error(f"Error processing result: {str(e)}", exc_info=True)
                if 'submission_id' in locals():
                    all_results.append({
                        'id': submission_id,
                        'error': str(e),
                        'status': 'failed'
                    })
                continue

        # Convert results to DataFrame and ensure consistent columns
        if len(all_results) > 0:
            logger.info(f"Processed {len(all_results)} submissions")
            try:
                results_df = pd.DataFrame(all_results)
                # Ensure all failed entries have NaN for numeric columns
                numeric_columns = results_df.select_dtypes(include=['float64', 'int64']).columns
                results_df.loc[results_df['status'] == 'failed', numeric_columns] = pd.NA
                return results_df.to_dict('records')
            except Exception as e:
                logger.error(f"Error converting results to DataFrame: {str(e)}", exc_info=True)
                return all_results
        else:
            logger.error("No submissions were processed successfully in batch")
            return []
        
    except Exception as e:
        logger.error(f"Error in batch processing: {str(e)}", exc_info=True)
        return False