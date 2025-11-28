import pandas as pd
import argparse
import os
from Bio.PDB import PDBParser, NeighborSearch

def find_contacting_residues(pdb_file, chain_A_id, chain_B_id, distance_cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    model = structure[0]  # assuming first model
    
    try:
        chain_A = model[chain_A_id]
        chain_B = model[chain_B_id]
    except KeyError:
        print(f"Warning: Chains {chain_A_id} or {chain_B_id} not found in {pdb_file}")
        return None

    atoms_A = list(chain_A.get_atoms())
    atoms_B = list(chain_B.get_atoms())
    
    neighbor_search = NeighborSearch(atoms_A + atoms_B)

    residues_in_A = set()
    residues_in_B = set()
    
    for atom_A in atoms_A:
        close_atoms = neighbor_search.search(atom_A.coord, distance_cutoff, level='A')
        for atom_B in close_atoms:
            if atom_B.get_parent().get_parent().id == chain_B_id:
                res_A = atom_A.get_parent()
                res_B = atom_B.get_parent()
                # Store both residues
                residues_in_A.add(res_A.get_id()[1])
                residues_in_B.add(res_B.get_id()[1])
    
    return residues_in_A, residues_in_B

def main():
    parser = argparse.ArgumentParser(description="Find contacting residues between two chains across multiple PDB files.")
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory containing PDB files")
    parser.add_argument("--chain_A_id", type=str,default="A", help="Chain ID for chain A")
    parser.add_argument("--chain_B_id", type=str,default="B", help="Chain ID for chain B")
    parser.add_argument("--distance_cutoff", type=float, default=4.0, help="Distance cutoff in angstroms (default: 4.0)")
    parser.add_argument("--output_csv", type=str, default="contact_residues_summary.csv", help="Output CSV filename (default: contact_residues_summary.csv)")
    
    args = parser.parse_args()
    
    pdb_files = [os.path.join(args.pdb_dir, f) for f in os.listdir(args.pdb_dir) if f.endswith('.pdb')]
    if not pdb_files:
        print(f"No PDB files found in directory {args.pdb_dir}")
        return

    summary_rows = []
    for idx,pdb_file in enumerate(pdb_files):
        print(f"Processing {pdb_file} ({idx+1}/{len(pdb_files)})")
        result = find_contacting_residues(pdb_file, args.chain_A_id, args.chain_B_id, args.distance_cutoff)
        if result is None:
            continue
        residues_in_A, residues_in_B = result
        
        # Format residue IDs with chain letter
        residues_in_A_formatted = [f"{args.chain_A_id}{res}" for res in sorted(residues_in_A)]
        residues_in_B_formatted = [f"{args.chain_B_id}{res}" for res in sorted(residues_in_B)]
        
        row = {
            'name': os.path.basename(pdb_file).replace(".pdb",""),
            'Chain_A_interface_residues': ",".join(residues_in_A_formatted),
            'Chain_B_interface_residues': ",".join(residues_in_B_formatted)
        }
        summary_rows.append(row)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(args.output_csv, index=False)
        print(f"Saved summarized contact residues to {args.output_csv}")
    else:
        print("No contacts found across all PDBs.")

if __name__ == "__main__":
    main()