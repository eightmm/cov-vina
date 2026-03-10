#!/usr/bin/env python3
"""Extract ligand from PDB file and save as SDF."""

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem


def extract_ligand(pdb_path, output_sdf, ligand_name=None):
    """
    Extract ligand from PDB file and save to SDF.

    Args:
        pdb_path: Path to PDB file
        output_sdf: Output SDF file path
        ligand_name: Specific ligand residue name (e.g., '02J'). If None, extracts first HETATM.
    """
    # Load PDB
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)

    if mol is None:
        raise ValueError(f"Failed to load PDB: {pdb_path}")

    # Find ligand atoms
    ligand_atoms = []
    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info:
            res_name = pdb_info.GetResidueName().strip()
            is_hetatm = pdb_info.GetIsHeteroAtom()

            # Skip water and common ions
            if res_name in ['HOH', 'WAT', 'NA', 'CL', 'MG', 'ZN', 'CA', 'K']:
                continue

            if is_hetatm:
                if ligand_name is None or res_name == ligand_name:
                    ligand_atoms.append(atom.GetIdx())

    if not ligand_atoms:
        raise ValueError(f"No ligand atoms found (ligand_name={ligand_name})")

    print(f"Found {len(ligand_atoms)} ligand atoms")

    # Extract ligand as new molecule
    ligand_mol = Chem.RWMol()
    old_to_new = {}

    for old_idx in ligand_atoms:
        atom = mol.GetAtomWithIdx(old_idx)
        new_idx = ligand_mol.AddAtom(atom)
        old_to_new[old_idx] = new_idx

    # Copy bonds
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx in old_to_new and end_idx in old_to_new:
            ligand_mol.AddBond(
                old_to_new[begin_idx],
                old_to_new[end_idx],
                bond.GetBondType()
            )

    # Copy conformer
    conf = mol.GetConformer()
    new_conf = Chem.Conformer(len(ligand_atoms))
    for old_idx, new_idx in old_to_new.items():
        pos = conf.GetAtomPosition(old_idx)
        new_conf.SetAtomPosition(new_idx, pos)

    ligand_mol = ligand_mol.GetMol()
    ligand_mol.AddConformer(new_conf)

    # Sanitize
    try:
        Chem.SanitizeMol(ligand_mol)
    except Exception as e:
        print(f"Warning: Sanitization failed: {e}")
        print("Attempting to save without sanitization...")

    # Save to SDF
    writer = Chem.SDWriter(output_sdf)
    writer.write(ligand_mol)
    writer.close()

    print(f"Saved ligand to: {output_sdf}")
    print(f"SMILES: {Chem.MolToSmiles(ligand_mol)}")


def main():
    parser = argparse.ArgumentParser(description="Extract ligand from PDB and save as SDF")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-o", "--output", required=True, help="Output SDF file")
    parser.add_argument("-l", "--ligand", help="Ligand residue name (e.g., '02J')")

    args = parser.parse_args()
    extract_ligand(args.input, args.output, args.ligand)


if __name__ == "__main__":
    main()
