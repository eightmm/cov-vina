#!/usr/bin/env python
"""
Calculate RMSD between crystal ligand and docked poses using MCS alignment.
"""
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign


def calculate_mcs_rmsd(ref_mol, query_mol, verbose=False):
    """
    Calculate RMSD between reference and query molecules using MCS alignment.

    Args:
        ref_mol: Reference molecule (e.g., crystal structure)
        query_mol: Query molecule (e.g., docked pose)
        verbose: Print detailed information

    Returns:
        rmsd: RMSD in Angstroms, or None if MCS fails
    """
    # Find maximum common substructure
    mcs_result = rdFMCS.FindMCS(
        [ref_mol, query_mol],
        timeout=2,
        completeRingsOnly=False,
        ringMatchesRingOnly=False,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
    )

    if mcs_result.numAtoms == 0:
        if verbose:
            print("No MCS found")
        return None

    # Get MCS molecule
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)

    if mcs_mol is None:
        if verbose:
            print(f"Invalid MCS SMARTS: {mcs_smarts}")
        return None

    # Get atom mappings
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)
    query_match = query_mol.GetSubstructMatch(mcs_mol)

    if len(ref_match) == 0 or len(query_match) == 0:
        if verbose:
            print("MCS match failed")
        return None

    # Extract coordinates
    ref_conf = ref_mol.GetConformer()
    query_conf = query_mol.GetConformer()

    ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in ref_match])
    query_coords = np.array([query_conf.GetAtomPosition(i) for i in query_match])

    # Calculate RMSD after optimal alignment
    # Use Kabsch algorithm to find optimal rotation

    # Center both structures
    ref_center = ref_coords.mean(axis=0)
    query_center = query_coords.mean(axis=0)
    ref_coords_centered = ref_coords - ref_center
    query_coords_centered = query_coords - query_center

    # Kabsch algorithm: find optimal rotation matrix
    # R = argmin ||P - QR||
    # Solution: R = V @ U.T where U, S, V = SVD(P.T @ Q)
    H = query_coords_centered.T @ ref_coords_centered
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    query_coords_aligned = query_coords_centered @ R

    # Calculate RMSD
    diff = ref_coords_centered - query_coords_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    if verbose:
        print(f"MCS atoms: {len(ref_match)}/{ref_mol.GetNumHeavyAtoms()} (ref), {len(query_match)}/{query_mol.GetNumHeavyAtoms()} (query)")
        print(f"RMSD: {rmsd:.3f} Å")

    return rmsd, len(ref_match)


def main():
    parser = argparse.ArgumentParser(description="Calculate RMSD between crystal and docked ligands")
    parser.add_argument("-r", "--reference", required=True, help="Reference PDB file (crystal ligand)")
    parser.add_argument("-q", "--query", required=True, help="Query SDF file (docked poses)")
    parser.add_argument("-n", "--top_n", type=int, default=10, help="Number of top poses to analyze")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    # Load reference
    print(f"Loading reference from {args.reference}...")
    ref_mol = Chem.MolFromPDBFile(args.reference, removeHs=False)
    if ref_mol is None:
        print(f"Error: Could not load reference from {args.reference}")
        return

    print(f"  Reference: {ref_mol.GetNumAtoms()} atoms, {ref_mol.GetNumHeavyAtoms()} heavy atoms")

    # Load docked poses
    print(f"\nLoading docked poses from {args.query}...")
    suppl = Chem.SDMolSupplier(args.query, removeHs=False)
    poses = [mol for mol in suppl if mol is not None]

    if len(poses) == 0:
        print(f"Error: No valid molecules in {args.query}")
        return

    print(f"  Loaded {len(poses)} poses")
    print(f"  Pose format: {poses[0].GetNumAtoms()} atoms, {poses[0].GetNumHeavyAtoms()} heavy atoms")

    # Calculate RMSD for top N poses
    print(f"\nCalculating RMSD for top {min(args.top_n, len(poses))} poses:")
    print(f"{'Rank':<6} {'Energy':<12} {'RMSD':<10} {'MCS atoms':<12}")
    print("-" * 50)

    for i in range(min(args.top_n, len(poses))):
        pose = poses[i]

        # Get energy
        try:
            energy = float(pose.GetProp('Energy'))
        except:
            energy = 0.0

        # Calculate RMSD
        result = calculate_mcs_rmsd(ref_mol, pose, verbose=args.verbose)

        if result is not None:
            rmsd, mcs_atoms = result
            print(f"{i+1:<6} {energy:<12.3f} {rmsd:<10.3f} {mcs_atoms:<12}")
        else:
            print(f"{i+1:<6} {energy:<12.3f} {'FAILED':<10} {'-':<12}")

    # Find best RMSD
    print("\n" + "=" * 50)
    rmsds = []
    for pose in poses[:args.top_n]:
        result = calculate_mcs_rmsd(ref_mol, pose, verbose=False)
        if result is not None:
            rmsds.append(result[0])

    if rmsds:
        best_rmsd = min(rmsds)
        best_idx = rmsds.index(best_rmsd)
        print(f"Best RMSD: {best_rmsd:.3f} Å (Pose {best_idx + 1})")

        # Get energy of best RMSD pose
        try:
            best_energy = float(poses[best_idx].GetProp('Energy'))
            print(f"  Energy: {best_energy:.3f} kcal/mol")
        except:
            pass


if __name__ == "__main__":
    main()
