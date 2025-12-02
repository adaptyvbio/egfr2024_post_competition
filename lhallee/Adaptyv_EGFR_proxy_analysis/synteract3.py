import torch
import pandas as pd
from transformers import AutoModel
from tqdm.auto import tqdm


"""
Synteract3 is a proprietary model from Synthyra.
If you'd like access to it, please Synthyra at info@synthyra.com
"""

DATA_DIR = 'data'
MODEL_DICT = {
    'synteract3_human': "Synthyra/###",
    'synteract3_multi': "Synthyra/###",
}


if __name__ == "__main__":
    # py -m synteract3
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    df = pd.read_csv(f'{DATA_DIR}/raw_data.csv')
    
    for model_name, model_path in MODEL_DICT.items():
        model = AutoModel.from_pretrained(model_path, trust_remote_code=True)
        model = model.eval().to(device)

        synteract_probs = []
        for seqa, seqb in tqdm(zip(df['SeqA'].tolist(), df['SeqB'].tolist()), total=len(df), desc=model_name):
            prob_ab = model.get_ppi_logits(seqa, seqb).detach().cpu().item()
            prob_ba = model.get_ppi_logits(seqb, seqa).detach().cpu().item()
            prob = (prob_ab + prob_ba) / 2
            synteract_probs.append(prob)

        df[model_name] = synteract_probs

        model.cpu()
        del model
        torch.cuda.empty_cache()

    df.to_csv(f'{DATA_DIR}/synteract3_out.csv', index=False)