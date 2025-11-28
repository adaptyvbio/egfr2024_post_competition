import os
import subprocess
import tempfile
import pandas as pd

def compute_tm_scores(pdb_dir, reference_pdb, foldseek_path='foldseek', tmp_dir=None):
    """
    Compute TM-scores between a reference PDB and all PDBs in a given directory using Foldseek.

    Args:
        pdb_dir (str): Path to directory containing PDB files.
        reference_pdb (str): Path to the reference PDB file.
        foldseek_path (str): Path to Foldseek executable (default assumes it's in PATH).
        tmp_dir (str): Path to temporary working directory (optional).

    Returns:
        pd.DataFrame: DataFrame with columns ['query', 'target', 'alntmscore']
    """
    # Make a temp directory
    work_dir = tempfile.mkdtemp(dir=tmp_dir)

    reference_db = os.path.join(work_dir, "reference_db")
    target_db = os.path.join(work_dir, "target_db")
    output_alignment = os.path.join(work_dir, "alignment.m8")

    # Step 1: Create Foldseek databases
    subprocess.run([foldseek_path, "createdb", reference_pdb, reference_db], check=True)
    subprocess.run([foldseek_path, "createdb", pdb_dir, target_db], check=True)

    # Step 2: Run Foldseek align
    subprocess.run([
        foldseek_path, "easy-search", reference_db, target_db, output_alignment, work_dir,
        "--format-output", "query,target,fident,alnlen,alntmscore", "--alignment-type", "1",
        "--prefilter-mode", "0", "-e", "1e5","--exhaustive-search", # Goes through all pairs
    ], check=True)

    # Step 3: Read the output
    results = pd.read_csv(output_alignment, sep='\t', header=None)
    results.columns = ['query', 'target', 'fident', 'alnlen', 'tm_score']
    return results
    
def post_process_foldseek_df(df):
    """
    Post process the dataframe output from foldseek for downstream analysis
    """
    ## Keep only chain A (these are the binders)
    df["target_chain"] = df["target"].apply(lambda x: x[-1])
    df_filtered = df[df["target_chain"]=="A"].copy().reset_index(drop=True)
    ## Remove the chain identifiers
    df_filtered["target"] = df_filtered["target"].apply(lambda x: x[:-2])
    df_filtered_subset = df_filtered[["target","tm_score"]]
    df_filtered_subset.columns = ["name","tm_score"]
    return df_filtered_subset


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute TM-score of PDBs against a reference using Foldseek.")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing PDB files to compare.")
    parser.add_argument("--reference_pdb", required=True, help="Path to the reference PDB file.")
    parser.add_argument("--out", required=True, help="Output csv file")
    parser.add_argument("--foldseek_path", default="foldseek", help="Path to foldseek binary (default assumes it's in PATH).")
    parser.add_argument("--tmp_dir", default=None, help="Temporary directory to use.")
    
    args = parser.parse_args()

    df = compute_tm_scores(args.pdb_dir, args.reference_pdb, args.foldseek_path, args.tmp_dir)
    df_processed = post_process_foldseek_df(df)
    print(df_processed)
    df_processed.to_csv(args.out,index=False)
    print(f"TM-scores to EGF wrote to {args.out}")