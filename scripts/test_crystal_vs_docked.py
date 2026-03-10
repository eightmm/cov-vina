"""
Compare docked poses with crystal structure.

This script:
1. Takes a query ligand SMILES
2. Generates random conformers
3. Runs covalent docking pipeline
4. Compares best pose with initial conformer (if reference provided)
"""

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import numpy as np

from cov_vina import run_covalent_pipeline


def calculate_rmsd(coords1, coords2):
    """Calculate RMSD between two coordinate sets."""
    return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))


def main():
    parser = argparse.ArgumentParser(description="Test covalent docking with random conformers")
    parser.add_argument("-q", "--query", required=True, help="Query ligand SMILES")
    parser.add_argument("-p", "--pocket", required=True, help="Protein PDB file")
    parser.add_argument("-r", "--residue", required=True, help="Reactive residue (e.g. CYS145)")
    parser.add_argument("-o", "--output", default="test_output", help="Output directory")
    parser.add_argument("--steps", type=int, default=100, help="Optimization steps")
    parser.add_argument("--num_confs", type=int, default=50, help="Number of conformers to generate")

    args = parser.parse_args()

    print(f"Query SMILES: {args.query}")
    print(f"Reactive residue: {args.residue}")
    print(f"Optimization steps: {args.steps}")
    print(f"Number of conformers: {args.num_confs}")

    # Run pipeline
    print("\n=== Running Covalent Docking Pipeline ===")
    results = run_covalent_pipeline(
        protein_pdb=args.pocket,
        query_ligand=args.query,
        reactive_residue=args.residue,
        output_dir=args.output,
        num_confs=args.num_confs,
        optimize=True,
        opt_steps=args.steps,
        verbose=True,
        pocket_cutoff=12.0,
    )

    print(f"\n=== Results ===")
    print(f"Best score: {results['best_score']:.3f} kcal/mol")
    print(f"Number of poses saved: {results['num_poses']}")

    # Load output SDF and show score distribution
    output_sdf = f"{args.output}/covalent_poses_all.sdf"
    suppl = Chem.SDMolSupplier(output_sdf, removeHs=False)

    scores = []
    initial_scores = []

    for mol in suppl:
        if mol is not None:
            if mol.HasProp("Vina_Score"):
                scores.append(float(mol.GetProp("Vina_Score")))
            if mol.HasProp("Vina_Score_Initial"):
                initial_scores.append(float(mol.GetProp("Vina_Score_Initial")))

    if scores:
        print(f"\nScore statistics:")
        print(f"  Best (lowest):  {min(scores):.3f} kcal/mol")
        print(f"  Worst (highest): {max(scores):.3f} kcal/mol")
        print(f"  Mean: {np.mean(scores):.3f} kcal/mol")
        print(f"  Std:  {np.std(scores):.3f} kcal/mol")

    if initial_scores:
        improvements = [init - final for init, final in zip(initial_scores, scores)]
        print(f"\nOptimization improvements:")
        print(f"  Best improvement: {max(improvements):.3f} kcal/mol")
        print(f"  Mean improvement: {np.mean(improvements):.3f} kcal/mol")
        print(f"  Poses improved: {sum(1 for x in improvements if x > 0)}/{len(improvements)}")

    print(f"\n✓ Results saved to {output_sdf}")


if __name__ == "__main__":
    main()
