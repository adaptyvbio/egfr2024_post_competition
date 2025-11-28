#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select

class ChainSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

def extract_and_reindex_chain(pdb_path, chain_id, output_path):
    """
    Extracts a specified chain from a PDB file, reindexes atom and residue numbers, and writes to a new PDB file.

    Args:
        pdb_path (str or Path): Path to the input PDB file.
        chain_id (str): Chain identifier to extract (e.g., "A").
        output_path (str or Path): Path to write the new reindexed PDB file.
    """
    pdb_path = Path(pdb_path)
    output_path = Path(output_path)
    
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)

    # Select only the desired chain
    chain = None
    for model in structure:
        for ch in model:
            if ch.id == chain_id:
                chain = ch
                break

    if chain is None:
        raise ValueError(f"Chain {chain_id} not found in the structure.")

    # Reindex atoms and residues
    for res_idx, residue in enumerate(chain.get_residues(), start=1):
        residue.id = (' ', res_idx, ' ')  # (' ', resseq, icode)

    # Save the chain
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), select=ChainSelector(chain_id))

    print(f"✅ Chain {chain_id} extracted, reindexed, and saved to {output_path}.")

def main():
    parser = argparse.ArgumentParser(description="Extract a chain from a PDB file and reindex atom and residue numbers.")
    parser.add_argument("--pdb", type=str, help="Input PDB file path")
    parser.add_argument("--chain", type=str, help="Chain ID to extract (e.g., 'A')")
    parser.add_argument("--output", type=str, help="Output PDB file path")
    args = parser.parse_args()

    extract_and_reindex_chain(args.pdb, args.chain, args.output)

if __name__ == "__main__":
    main()