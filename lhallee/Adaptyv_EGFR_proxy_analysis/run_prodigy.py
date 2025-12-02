import subprocess
import os
import pandas as pd
import re
import numpy as np
from glob import glob
from tqdm.auto import tqdm
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.MMCIFParser import MMCIFParser
from concurrent.futures import ProcessPoolExecutor


def process_structure(file: str) -> tuple:
    """Run prodigy on a single PDB file and extract sequences.

    Returns a tuple: (SeqA, SeqB, kd)
    """
    # Always normalize to forward slashes for subprocess arg
    file = file.replace('\\', '/')

    # Run prodigy CLI as a subprocess
    env = os.environ.copy()
    # Force UTF-8 for child process so prints/logging are encodable
    env["PYTHONIOENCODING"] = "utf-8"
    result = subprocess.run(
        ["py", "-m", "prodigy_local.cli", file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # merge stderr into stdout to avoid missing lines
        text=True,
        encoding="utf-8",
        errors="replace",  # avoid decode failures on special characters
        env=env,
    )
    output_str = result.stdout

    # Extract KD robustly; ignore degree symbol variations and whitespace
    # Example line:
    # [++] predicted dissociation constant (M) at 25.0˚C:  1.23e-06
    m = re.search(r"dissociation constant \(M\) at\s*[0-9.]+[^\n]*C:\s*([0-9.eE+-]+)", output_str)
    if not m:
        raise ValueError(f"Could not parse KD from prodigy output for {file}.\nOutput was:\n{output_str}")
    kd = float(m.group(1))
    pkd = -np.log10(kd)
    # Parse sequences from structure file (PDB or mmCIF)
    parser = MMCIFParser() if file.lower().endswith('.cif') else PDBParser(QUIET=True)
    ppb = PPBuilder()
    structure = parser.get_structure('structure', file)
    chains = list(structure.get_chains())
    
    # Extract sequences from first two chains using PPBuilder
    seq_a = ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(chains[0])])
    seq_b = ''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(chains[1])])

    return (seq_a, seq_b, pkd)


if __name__ == "__main__":
    # py -m run_prodigy
    af2_pdb_files = list(glob("./af2_preds/*.pdb"))
    chai_cif_files = list(glob("./chai1_preds/*.cif"))

    file_dict = {
        'af2': af2_pdb_files,
        'chai': chai_cif_files,
    }

    # Choose a reasonable number of workers
    max_workers = max(1, os.cpu_count() - 2)
    for name, files in file_dict.items():
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # executor.map preserves input order
            ordered_results = list(
                tqdm(
                    executor.map(process_structure, files, chunksize=1),
                    total=len(files)
                )
            )

        # save to df
        seqs_a = [r[0] for r in ordered_results]
        seqs_b = [r[1] for r in ordered_results]
        kds = [r[2] for r in ordered_results]
        df = pd.DataFrame({'SeqA': seqs_a, 'SeqB': seqs_b, f'prodigy_ppkd_{name}': kds})
        df.to_csv(f'prodigy_ppkd_{name}.csv', index=False)