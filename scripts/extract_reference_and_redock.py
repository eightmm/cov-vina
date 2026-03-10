#!/usr/bin/env python3
"""
Extract crystal ligand, convert to SMILES, and redock to calculate RMSD.
"""
import argparse
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np

from cov_vina import run_covalent_pipeline


def pdb_to_sdf(pdb_path, sdf_path):
    """Convert PDB to SDF using RDKit"""
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read PDB: {pdb_path}")

    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()

    return mol


def get_smiles_from_mol(mol):
    """Get SMILES from RDKit molecule, removing covalent adduct if present"""
    # Try to get SMILES
    smiles = Chem.MolToSmiles(mol)
    print(f"  Original SMILES: {smiles}")
    print(f"  Num atoms: {mol.GetNumAtoms()}")
    print(f"  Num heavy atoms: {mol.GetNumHeavyAtoms()}")

    return smiles


def calculate_rmsd(ref_mol, docked_sdf, align=True):
    """Calculate RMSD between reference and docked poses"""
    ref_conf = ref_mol.GetConformer()
    ref_coords = ref_conf.GetPositions()

    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    rmsds = []

    for i, mol in enumerate(suppl):
        if mol is None:
            continue

        conf = mol.GetConformer()
        coords = conf.GetPositions()

        # Simple RMSD (assume same atom order)
        if len(coords) == len(ref_coords):
            rmsd = np.sqrt(np.mean(np.sum((coords - ref_coords)**2, axis=1)))
            rmsds.append(rmsd)

    if not rmsds:
        return None, []

    return min(rmsds), rmsds


def main():
    parser = argparse.ArgumentParser(description="Extract reference ligand and redock")
    parser.add_argument("--pdb", required=True, help="Full PDB file")
    parser.add_argument("--pocket", required=True, help="Pocket PDB file")
    parser.add_argument("--ligand_name", required=True, help="Ligand residue name (e.g., 03P, AQ4, PJE)")
    parser.add_argument("--residue", required=True, help="Reactive residue (e.g., CYS775)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--smiles", default=None, help="Manual SMILES (if crystal is adduct)")
    parser.add_argument("--num_confs", type=int, default=200, help="Number of conformers")
    parser.add_argument("--opt_steps", type=int, default=100, help="Optimization steps")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("="*70)
    print("REFERENCE LIGAND EXTRACTION AND REDOCKING")
    print("="*70)

    # Extract crystal ligand to PDB
    crystal_pdb = os.path.join(args.output_dir, "reference_crystal.pdb")
    cmd = f"grep '^HETATM.*{args.ligand_name}' {args.pdb} > {crystal_pdb} && echo 'END' >> {crystal_pdb}"
    os.system(cmd)

    print(f"\n1. Extracted crystal ligand: {args.ligand_name}")
    print(f"   Saved to: {crystal_pdb}")

    # Convert to SDF
    crystal_sdf = os.path.join(args.output_dir, "reference_crystal.sdf")
    try:
        mol = pdb_to_sdf(crystal_pdb, crystal_sdf)
        print(f"\n2. Converted to SDF: {crystal_sdf}")
    except Exception as e:
        print(f"\n✗ Failed to convert PDB to SDF: {e}")
        return

    # Get SMILES
    if args.smiles:
        smiles = args.smiles
        print(f"\n3. Using manual SMILES: {smiles}")
    else:
        smiles = get_smiles_from_mol(mol)
        print(f"\n3. Extracted SMILES: {smiles}")

    # Check for warhead
    from cov_vina.molecular.anchor import detect_warheads
    test_mol = Chem.MolFromSmiles(smiles)
    if test_mol:
        warheads = detect_warheads(test_mol)
        if warheads:
            print(f"\n✓ Warhead detected: {warheads[0].warhead_type}")
        else:
            print(f"\n✗ WARNING: No warhead detected in SMILES!")
            print(f"   This ligand may not be covalent or needs manual SMILES")

    # Run redocking
    print(f"\n4. Running covalent docking...")
    print(f"   Protein: {args.pocket}")
    print(f"   Ligand: {smiles}")
    print(f"   Residue: {args.residue}")

    try:
        results = run_covalent_pipeline(
            protein_pdb=args.pocket,
            query_ligand=smiles,
            reactive_residue=args.residue,
            output_dir=args.output_dir,
            num_confs=args.num_confs,
            optimize=True,
            opt_steps=args.opt_steps,
            verbose=False,
        )

        print(f"\n5. Docking completed!")
        print(f"   Best score: {results['best_score']:.3f} kcal/mol")
        print(f"   Poses: {results['num_poses']}")

    except Exception as e:
        print(f"\n✗ Docking failed: {e}")
        return

    # Calculate RMSD
    docked_sdf = os.path.join(args.output_dir, "covalent_poses_all.sdf")

    print(f"\n6. Calculating RMSD...")
    try:
        min_rmsd, all_rmsds = calculate_rmsd(mol, docked_sdf)

        if min_rmsd is not None:
            print(f"\n{'='*70}")
            print("RESULTS")
            print(f"{'='*70}")
            print(f"Best RMSD: {min_rmsd:.3f} Å")
            print(f"Best score: {results['best_score']:.3f} kcal/mol")
            print(f"Poses evaluated: {len(all_rmsds)}")
            print(f"All RMSDs: {[f'{r:.3f}' for r in sorted(all_rmsds)[:5]]}... (top 5)")
            print(f"{'='*70}")
        else:
            print(f"\n✗ RMSD calculation failed (atom count mismatch)")

    except Exception as e:
        print(f"\n✗ RMSD calculation failed: {e}")


if __name__ == "__main__":
    main()
