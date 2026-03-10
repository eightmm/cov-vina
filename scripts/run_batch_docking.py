#!/usr/bin/env python3
"""
Batch covalent docking from .smi file

Reads molecules from a .smi file and runs covalent docking for each.
Output format: {output_dir}/{molecule_name}/final_poses.sdf
"""
import argparse
import os
from pathlib import Path

from cov_vina import run_covalent_pipeline


def parse_smi_file(smi_path: str):
    """
    Parse .smi file and return list of (SMILES, name) tuples.

    Format:
        SMILES_string    molecule_name
        # Comments starting with # are ignored
    """
    molecules = []
    with open(smi_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                smiles = parts[0]
                name = parts[1]
                molecules.append((smiles, name))
            elif len(parts) == 1:
                # No name provided, use index
                smiles = parts[0]
                name = f"molecule_{len(molecules)+1}"
                molecules.append((smiles, name))

    return molecules


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch covalent docking from .smi file"
    )

    # Required
    parser.add_argument("-p", "--protein", required=True,
                        help="Path to protein pocket PDB file")
    parser.add_argument("-s", "--smi_file", required=True,
                        help="Path to .smi file with SMILES and names")
    parser.add_argument("-r", "--reactive_residue", required=True,
                        help="Reactive residue (e.g., 'CYS145')")

    # Output
    parser.add_argument("-o", "--out_dir", default="batch_results",
                        help="Base output directory (default: batch_results)")

    # Docking parameters
    parser.add_argument("-n", "--num_confs", type=int, default=200,
                        help="Number of conformers per molecule (default: 200)")
    parser.add_argument("--opt_steps", type=int, default=100,
                        help="Optimization steps (default: 100)")
    parser.add_argument("--optimize", action="store_true", default=True,
                        help="Enable gradient optimization (default: True)")
    parser.add_argument("--save_all", action="store_true",
                        help="Save all poses instead of top poses")

    # Visualization
    parser.add_argument("--make_gif", action="store_true",
                        help="Generate trajectory.gif for each molecule")

    args = parser.parse_args()

    # Parse molecules
    molecules = parse_smi_file(args.smi_file)
    print(f"Found {len(molecules)} molecules in {args.smi_file}")

    # Run docking for each molecule
    for i, (smiles, name) in enumerate(molecules, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(molecules)}] Processing: {name}")
        print(f"SMILES: {smiles}")
        print('='*60)

        # Create output directory
        mol_out_dir = os.path.join(args.out_dir, name)
        os.makedirs(mol_out_dir, exist_ok=True)

        try:
            results = run_covalent_pipeline(
                protein_pdb=args.protein,
                query_ligand=smiles,
                reactive_residue=args.reactive_residue,
                output_dir=mol_out_dir,
                num_confs=args.num_confs,
                optimize=args.optimize,
                opt_steps=args.opt_steps,
                save_all_poses=args.save_all if args.save_all else None,
                verbose=False,  # Reduce clutter
            )

            # Rename output file to standard format
            old_path = os.path.join(mol_out_dir, "covalent_poses_all.sdf")
            new_path = os.path.join(mol_out_dir, "final_poses.sdf")
            if os.path.exists(old_path):
                os.rename(old_path, new_path)

            print(f"✓ Success: {name}")
            print(f"  Best score: {results['best_score']:.3f} kcal/mol")
            print(f"  Warhead: {results['warhead_type']}")
            print(f"  Output: {new_path}")

        except Exception as e:
            print(f"✗ Failed: {name}")
            print(f"  Error: {e}")
            continue

    print(f"\n{'='*60}")
    print("Batch docking complete!")
    print(f"Results saved to: {args.out_dir}/")
    print('='*60)
