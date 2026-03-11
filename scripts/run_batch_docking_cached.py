#!/usr/bin/env python3
"""
Batch covalent docking with pocket caching for efficiency.

Key optimization:
- Load pocket once and reuse for all ligands
- Avoids redundant pocket feature computation
- ~2-3x speedup for multiple ligands on same target
"""
import argparse
import os
import time
from pathlib import Path

from cov_vina import run_covalent_pipeline, load_pocket_for_caching


def parse_smi_file(smi_path: str):
    """Parse .smi file and return list of (SMILES, name) tuples."""
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
                smiles = parts[0]
                name = f"molecule_{len(molecules)+1}"
                molecules.append((smiles, name))

    return molecules


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch covalent docking with pocket caching"
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

    # Performance
    parser.add_argument("--pocket_cutoff", type=float, default=12.0,
                        help="Pocket extraction cutoff (Å, default: 12.0)")

    args = parser.parse_args()

    # Parse molecules
    molecules = parse_smi_file(args.smi_file)
    print(f"Found {len(molecules)} molecules in {args.smi_file}")

    # Pre-load pocket once (OPTIMIZATION: 2-3x speedup)
    print(f"\nPre-loading protein pocket: {args.protein}")
    print(f"Reactive residue: {args.reactive_residue}")

    cache_start = time.time()
    cached_pocket = load_pocket_for_caching(
        protein_pdb=args.protein,
        reactive_residue=args.reactive_residue,
        pocket_cutoff=args.pocket_cutoff,
        verbose=True,
    )
    cache_time = time.time() - cache_start
    print(f"✓ Pocket loaded and cached in {cache_time:.2f}s")
    print(f"  Device: {cached_pocket['device']}")
    print(f"  Pocket atoms: {cached_pocket['pocket_bundle'].mol.GetNumAtoms()}")

    # Run docking for each molecule
    total_start = time.time()
    successful = 0
    failed = 0

    for i, (smiles, name) in enumerate(molecules, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(molecules)}] Processing: {name}")
        print(f"SMILES: {smiles}")
        print('='*60)

        # Create output directory
        mol_out_dir = os.path.join(args.out_dir, name)
        os.makedirs(mol_out_dir, exist_ok=True)

        mol_start = time.time()

        try:
            results = run_covalent_pipeline(
                protein_pdb=args.protein,  # Still needed for pocket.pdb output
                query_ligand=smiles,
                reactive_residue=args.reactive_residue,
                output_dir=mol_out_dir,
                pocket_cutoff=args.pocket_cutoff,
                num_confs=args.num_confs,
                optimize=args.optimize,
                opt_steps=args.opt_steps,
                save_all_poses=args.save_all if args.save_all else None,
                _cached_pocket=cached_pocket,  # Use cached pocket!
                verbose=True,  # Show progress for each ligand
            )

            # Rename output file to standard format
            old_path = os.path.join(mol_out_dir, "covalent_poses_all.sdf")
            new_path = os.path.join(mol_out_dir, "final_poses.sdf")
            if os.path.exists(old_path):
                os.rename(old_path, new_path)

            mol_time = time.time() - mol_start
            successful += 1

            print(f"✓ Success: {name}")
            print(f"  Best score: {results['best_score']:.3f} kcal/mol")
            print(f"  Warhead: {results['warhead_type']}")
            print(f"  Runtime: {mol_time:.2f}s")
            print(f"  Output: {new_path}")

        except Exception as e:
            mol_time = time.time() - mol_start
            failed += 1

            print(f"✗ Failed: {name}")
            print(f"  Error: {e}")
            print(f"  Runtime: {mol_time:.2f}s")

            # Save error log
            error_log = os.path.join(mol_out_dir, "error.log")
            with open(error_log, 'w') as f:
                f.write(f"SMILES: {smiles}\n")
                f.write(f"Error: {str(e)}\n")
            continue

    total_time = time.time() - total_start

    print(f"\n{'='*60}")
    print("Batch docking complete!")
    print(f"  Successful: {successful}/{len(molecules)}")
    print(f"  Failed: {failed}/{len(molecules)}")
    print(f"  Total time: {total_time:.2f}s")
    print(f"  Avg time/ligand: {total_time/len(molecules):.2f}s")
    print(f"Results saved to: {args.out_dir}/")
    print('='*60)
