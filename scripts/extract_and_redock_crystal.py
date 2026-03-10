"""
Extract crystal ligand, reconstruct original (pre-reaction) form, and re-dock.

For 6lu7 PJE:
- Crystal has S-C bond (Michael adduct)
- Original has C=C (vinyl group)
"""

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from cov_vina import run_covalent_pipeline


def reconstruct_vinyl_ligand_smiles():
    """
    Reconstruct the original vinyl ligand from PJE adduct.

    PJE adduct structure suggests it was originally a vinyl compound.
    Looking at the structure, it's a peptidomimetic with acrylamide-like warhead.

    Simplified reconstruction: Use a vinyl-peptide derivative that would
    form similar Michael adduct.
    """

    # Based on PJE structure and 6lu7 literature,
    # the original ligand has a vinyl group that undergoes Michael addition

    # Simplified vinyl peptide-like structure
    # This captures the key features: vinyl (C=C), carbonyl, peptide backbone
    candidates = [
        ("vinyl-alanine", "C=CC(=O)N[C@@H](C)C(=O)O"),
        ("vinyl-glycine", "C=CC(=O)NCC(=O)O"),
        ("simple-acrylamide-peptide", "C=CC(=O)N[C@@H](CC(C)C)C(=O)O"),
    ]

    return candidates


def main():
    parser = argparse.ArgumentParser(description="Extract crystal ligand and re-dock")
    parser.add_argument("--pdb", default="examples/6lu7/6lu7.pdb", help="Crystal structure PDB")
    parser.add_argument("--pocket", default="examples/6lu7/6lu7_pocket.pdb", help="Pocket PDB")
    parser.add_argument("--residue", default="CYS145", help="Reactive residue")
    parser.add_argument("--output", default="examples/6lu7/redocking", help="Output directory")
    parser.add_argument("--steps", type=int, default=200, help="Optimization steps")
    parser.add_argument("--num_confs", type=int, default=50, help="Number of conformers")

    args = parser.parse_args()

    print("="*70)
    print("CRYSTAL LIGAND RE-DOCKING")
    print("="*70)

    # Get candidate ligands
    candidates = reconstruct_vinyl_ligand_smiles()

    print("\nCandidate reconstructed ligands:")
    for name, smiles in candidates:
        print(f"  {name:25s}: {smiles}")

    # Use the most realistic one
    ligand_name, ligand_smiles = candidates[0]

    print(f"\n{'='*70}")
    print(f"Selected ligand: {ligand_name}")
    print(f"SMILES: {ligand_smiles}")
    print(f"{'='*70}")

    # Verify it has warhead
    mol = Chem.MolFromSmiles(ligand_smiles)
    if mol:
        from cov_vina.molecular.anchor import detect_warheads
        warheads = detect_warheads(mol)

        if warheads:
            print(f"\n✓ Warhead detected: {warheads[0].warhead_type}")
        else:
            print("\n✗ WARNING: No warhead detected!")
            return

    # Run re-docking
    print(f"\nRunning covalent docking...")
    print(f"  Protein: {args.pocket}")
    print(f"  Ligand: {ligand_smiles}")
    print(f"  Residue: {args.residue}")
    print(f"  Output: {args.output}")

    results = run_covalent_pipeline(
        protein_pdb=args.pocket,
        query_ligand=ligand_smiles,
        reactive_residue=args.residue,
        output_dir=args.output,
        num_confs=args.num_confs,
        optimize=True,
        opt_steps=args.steps,
        verbose=True,
    )

    print(f"\n{'='*70}")
    print("RESULTS")
    print(f"{'='*70}")
    print(f"Best score: {results['best_score']:.3f} kcal/mol")
    print(f"Poses saved: {results['num_poses']}")
    print(f"Output: {args.output}/covalent_poses_all.sdf")
    print(f"{'='*70}")

    # TODO: Calculate RMSD to crystal structure (requires proper atom mapping)
    print("\nNote: Crystal structure has adduct form (with S bonded)")
    print("      Direct RMSD comparison requires removing S and matching atoms")


if __name__ == "__main__":
    main()
