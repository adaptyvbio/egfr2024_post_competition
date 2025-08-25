import os

def clean_pdb(pdb_file):
    """Clean PDB file to keep only relevant lines"""
    with open(pdb_file, 'r') as f_in:
        relevant_lines = [line for line in f_in if line.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END'))]
    
    with open(pdb_file, 'w') as f_out:
        f_out.writelines(relevant_lines)

def add_cryst1_record(pdb_file):
    """Add CRYST1 record to PDB file if missing"""
    with open(pdb_file, 'r') as f:
        content = f.read()
    
    if not content.startswith('CRYST1'):
        cryst1_record = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"
        with open(pdb_file, 'w') as f:
            f.write(cryst1_record + content)
