#!/usr/bin/env python3
"""
Extract pocket around reactive residue from PDB file
"""
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def extract_pocket(pdb_path: str, residue: str, cutoff: float = 12.0, output_path: str = None):
    """
    Extract pocket residues within cutoff distance of reactive residue.

    Args:
        pdb_path: Path to PDB file
        residue: Residue specifier (e.g., 'CYS145' or 'CYS145:A')
        cutoff: Distance cutoff in Angstroms (default: 12.0)
        output_path: Output PDB path (optional)
    """
    # Parse residue specifier
    if ':' in residue:
        res_name_num, chain = residue.split(':')
    else:
        res_name_num = residue
        chain = None

    # Extract residue number
    res_num = int(''.join(c for c in res_name_num if c.isdigit()))
    res_name = ''.join(c for c in res_name_num if not c.isdigit())

    print(f"Extracting pocket around {res_name}{res_num}" + (f" chain {chain}" if chain else ""))
    print(f"Cutoff: {cutoff}Å")

    # Read PDB file
    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    # Find reactive residue atoms
    reactive_coords = []
    for line in lines:
        if not line.startswith('ATOM'):
            continue

        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_num = int(line[22:26].strip())

        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        # Check if this is the reactive residue
        if residue_num == res_num and residue_name == res_name:
            if chain is None or chain_id == chain:
                reactive_coords.append([x, y, z])

    if not reactive_coords:
        raise ValueError(f"Reactive residue {residue} not found in {pdb_path}")

    reactive_coords = np.array(reactive_coords)
    reactive_center = reactive_coords.mean(axis=0)

    print(f"Found reactive residue at: {reactive_center}")
    print(f"Number of atoms in reactive residue: {len(reactive_coords)}")

    # Extract pocket atoms within cutoff
    pocket_lines = []
    pocket_residues = set()

    # Keep HEADER, TITLE, etc.
    for line in lines:
        if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK')):
            pocket_lines.append(line)

    # Process ATOM lines
    for line in lines:
        if not line.startswith('ATOM'):
            continue

        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        coord = np.array([x, y, z])
        dist = np.linalg.norm(coord - reactive_center)

        if dist <= cutoff:
            pocket_lines.append(line)

            # Track unique residues
            residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            residue_num = int(line[22:26].strip())
            pocket_residues.add((chain_id, residue_name, residue_num))

    # Add END
    pocket_lines.append("END\n")

    print(f"\nPocket extraction complete:")
    print(f"  Total atoms: {sum(1 for l in pocket_lines if l.startswith('ATOM'))}")
    print(f"  Total residues: {len(pocket_residues)}")

    # Save if output path provided
    if output_path:
        with open(output_path, 'w') as f:
            f.writelines(pocket_lines)
        print(f"  Saved to: {output_path}")

    return pocket_lines


def main():
    parser = argparse.ArgumentParser(description="Extract pocket from PDB")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--residue", required=True, help="Reactive residue (e.g., CYS145)")
    parser.add_argument("--cutoff", type=float, default=12.0, help="Distance cutoff (Å)")
    parser.add_argument("--output", required=True, help="Output PDB file")

    args = parser.parse_args()

    extract_pocket(args.pdb, args.residue, args.cutoff, args.output)


if __name__ == "__main__":
    main()
