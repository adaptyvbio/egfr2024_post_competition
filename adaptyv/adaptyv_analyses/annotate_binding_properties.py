import os
import sys
import pandas as pd
from loguru import logger
from modal import App, Image, Secret, Volume, Mount
from tqdm import tqdm
import multiprocess as mp
import functools
from pathlib import Path

logger.add("annotate_competition.log", rotation="10 MB", level="INFO")
logger.add(sys.stderr, format="{time} {level} {message}", level="INFO")
logger.info("Starting annotation setup")

from src.setup.pyrosetta import set_up_pyrosetta, initialize_pyrosetta
from src.utils.scores import score_interface
from src.utils.dssp import calc_ss_percentage
from src.utils.pdb import clean_pdb, add_cryst1_record
from src.utils.relax import pr_relax
from src.utils.interface import analyze_interface_contacts
from src.processing.parse import get_pdb_files, check_missing_relaxed_structures
from src.processing.batch import process_submissions_batch

competition_vol = Volume.from_name("competition-output", create_if_missing=True)
structures_mount_1 = Mount.from_local_dir(
    "/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/round_1/structures/alphafold",
    remote_path="/root/competition_data/001_2024"
)
structures_mount_2 = Mount.from_local_dir(
    "/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/round_2/structures/alphafold",
    remote_path="/root/competition_data/002_2024"
)
sequences_mount = Mount.from_local_dir(
    "/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/processed",
    remote_path="/root/competition_data"
)
TIMEOUT = 3600 * 12  # 12 hours

image = (
    Image.debian_slim()
    .apt_install("git", "wget", "aria2", "ffmpeg", "dssp")
    .pip_install(
        "pdb-tools==2.4.8", 
        "ffmpeg-python==0.2.0",
        "plotly==5.18.0",
        "kaleido==0.2.1",
        "pyarrow",
        "fastparquet",
        "boto3",
        "python-dotenv",
        "loguru",
        "biopython",
        "numpy",
        "pandas",
        "tqdm",
        "scipy",
        "multiprocess",
    )
    .pip_install("git+https://github.com/sokrypton/ColabDesign.git")
    .pip_install("pyrosettacolabsetup")
    .run_function(set_up_pyrosetta)
    #.run_function(initialize_pyrosetta)
    .run_commands(
        "git clone https://github.com/martinpacesa/BindCraft /root/bindcraft"
        " && chmod +x /root/bindcraft/functions/dssp"
        " && chmod +x /root/bindcraft/functions/DAlphaBall.gcc"
    )
    .pip_install("jax[cuda]")
)

app = App("annotate-competition", image=image)

@app.function(
    image=image,
    cpu=32,
    timeout=TIMEOUT,
    volumes={
        "/root/competition_output": competition_vol
    },
    mounts=[structures_mount_1, structures_mount_2, sequences_mount]
)
def run_batch_processing(input_csv: str, structure_dirs: list, start_idx: int, end_idx: int, save_path: str, data_path: str):
    """Process a batch of submissions using Modal"""
    logger.info(f"Processing batch {start_idx}-{end_idx}")
       
    # Initialize PyRosetta at the start of processing
    initialize_pyrosetta()
    
    # Check if required directories exist
    logger.info("Checking required directories...")
    
    # Check data path and structure directories
    if not os.path.exists(data_path):
        raise RuntimeError(f"Data path does not exist: {data_path}")
        
    for dir in structure_dirs:
        full_path = os.path.join(data_path, dir)
        if not os.path.exists(full_path):
            raise RuntimeError(f"Structure directory does not exist: {full_path}")
    
    # Check save path and create if needed
    if not os.path.exists(save_path):
        logger.info(f"Creating save directory: {save_path}")
        os.makedirs(save_path, exist_ok=True)

    if not os.path.exists(os.path.join(save_path, "relaxed")):
        logger.info(f"Creating relaxed directory: {save_path}/relaxed")
        os.makedirs(os.path.join(save_path, "relaxed"), exist_ok=True)
        
    if not os.path.exists(os.path.join(save_path, "results")):
        logger.info(f"Creating results directory: {save_path}/results")
        os.makedirs(os.path.join(save_path, "results"), exist_ok=True)
        
    # Check input CSV exists
    input_csv_path = os.path.join(data_path, input_csv)
    if not os.path.exists(input_csv_path):
        raise RuntimeError(f"Input CSV file does not exist: {input_csv_path}")
        
    logger.info("All required directories and files exist")
    
    # Save batch results to Modal volume
    batch_results = process_submissions_batch(input_csv, structure_dirs, start_idx, end_idx, save_path, data_path)
    batch_output_path = os.path.join(save_path, "results", f"batch_{start_idx}_{end_idx}.csv")
    pd.DataFrame(batch_results).to_csv(batch_output_path, index=False)
    logger.info(f"Saved batch results to {batch_output_path}")
    
    return batch_output_path  # Return the path instead of the results

@app.local_entrypoint()
def main():
    """Main entry point for the Modal app"""
    logger.info("Starting annotation pipeline")
    
    try:
        # Read submissions file to get total count
        submissions = pd.read_csv('/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/processed/all_submissions.csv')
        total_submissions = len(submissions)
        batch_size = 64
        
        logger.info(f"Found {total_submissions} total submissions")
        
        # Process batches and collect paths
        batch_paths = []
        total_batches = (total_submissions + batch_size - 1) // batch_size
        
        with tqdm(total=total_batches, desc="Processing batches") as batch_pbar:
            for start_idx in range(0, total_submissions, batch_size):
                end_idx = min(start_idx + batch_size, total_submissions)
                logger.info(f"Processing batch {start_idx}-{end_idx}")
                
                # Run batch processing and collect the path
                batch_path = run_batch_processing.remote(
                    input_csv="all_submissions.csv",
                    structure_dirs=['001_2024', '002_2024'],
                    start_idx=start_idx,
                    end_idx=end_idx,
                    save_path="/root/competition_output",
                    data_path="/root/competition_data"
                )
                batch_paths.append(batch_path)
                batch_pbar.update(1)
        
        # Create a new Modal function to combine results
        @app.function(
            volumes={"/root/competition_output": competition_vol}
        )
        def combine_results(batch_paths):
            all_results = []
            for path in batch_paths:
                try:
                    batch_df = pd.read_csv(path)
                    all_results.append(batch_df)
                except Exception as e:
                    logger.error(f"Error reading batch results from {path}: {str(e)}")
            
            if not all_results:
                raise RuntimeError("No results were successfully processed")
                
            final_results = pd.concat(all_results, ignore_index=True)
            output_path = "/root/competition_output/results/final_results.csv"
            final_results.to_csv(output_path, index=False)
            return output_path
        
        # Combine results inside Modal
        final_path = combine_results.remote(batch_paths)
        logger.info(f"Successfully processed all submissions. Results saved to {final_path}")
            
    except Exception as e:
        logger.error(f"Error in main: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()