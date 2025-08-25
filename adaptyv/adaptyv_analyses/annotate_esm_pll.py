import os
import sys
import pandas as pd
import torch
import math
import esm
from loguru import logger
from modal import App, Image, Secret, Volume, Mount
from tqdm import tqdm

# Set up logging
logger.add("annotate_competition.log", rotation="10 MB", level="INFO")
logger.add(sys.stderr, format="{time} {level} {message}", level="INFO")
logger.info("Starting annotation setup")

def compute_pll(sequence, model, alphabet, device, batch_converter):
    """Compute PLL score for a sequence using the provided model and settings"""
    data = [("protein", sequence)]
    *_, batch_tokens = batch_converter(data)
    log_probs = []
    for i in range(len(sequence)):
        batch_tokens_masked = batch_tokens.clone()
        batch_tokens_masked[0, i + 1] = alphabet.mask_idx
        with torch.no_grad():
            token_probs = torch.log_softmax(
                model(batch_tokens_masked.to(device))["logits"], dim=-1
            )
        log_probs.append(token_probs[0, i + 1, alphabet.get_idx(sequence[i])].item())
    return math.fsum(log_probs)

# Image definition
image = (
    Image.debian_slim()
    .pip_install(
        "fair-esm",
        "torch", 
        "pandas",
        "tqdm",
        "loguru"
    )
)

app = App("compute-esm-pll", image=image)

# Define volume and mount
competition_vol = Volume.from_name("competition-output", create_if_missing=True)
sequences_mount = Mount.from_local_dir(
    "/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/processed/submissions",
    remote_path="/root/competition_data"
)

@app.function(
    image=image,
    gpu="A10G",
    timeout=3600,
    volumes={"/root/competition_output": competition_vol},
    mounts=[sequences_mount]
)
def process_submissions(input_csv: str):
    """Process round 1 submissions and compute ESM2 PLL scores"""
    try:
        # Set up device and model inside the function
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info(f"Using device: {device}")
        
        # Load model
        logger.info("Loading ESM2 model...")
        model, alphabet = esm.pretrained.load_model_and_alphabet("esm2_t33_650M_UR50D")
        model = model.to(device)
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        logger.info("Model loaded successfully")

        # Read submissions from mounted path
        df = pd.read_csv(os.path.join("/root/competition_data", input_csv))
        
        # Filter for round 1
        df_round1 = df[df['round'] == 1].copy()
        logger.info(f"Processing {len(df_round1)} submissions from round 1")

        # Compute PLL scores
        plls = []
        for seq in tqdm(df_round1['sequence'], desc="Computing PLLs"):
            pll = compute_pll(seq, model, alphabet, device, batch_converter)
            plls.append(pll)

        # Add scores to dataframe
        df_round1['esm_pll'] = plls
        
        # Save results to volume
        output_path = '/root/competition_output/round1_esm_pll_scores.csv'
        df_round1.to_csv(output_path, index=False)
        logger.info(f"Results saved to volume at {output_path}")
        
        return output_path
        
    except Exception as e:
        logger.error(f"Error in processing: {str(e)}", exc_info=True)
        raise e

@app.local_entrypoint()
def main():
    """Main entry point"""
    try:
        input_csv = 'all_submissions.csv'
        output_path = process_submissions.remote(input_csv)
        logger.info(f"Processing complete. Results saved to {output_path}")
    except Exception as e:
        logger.error(f"Error in main: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()