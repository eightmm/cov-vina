"""
Reconstruct the original (pre-reaction) ligand from covalent adduct in crystal structure.

For Michael addition: S-C-C → C=C (remove S, add double bond)
For SN2: S-C → Cl-C (remove S, add leaving group)
"""

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def extract_pje_from_pdb(pdb_file):
    """Extract PJE residue from PDB file."""
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)

    # Find PJE atoms
    pje_atoms = []
    conf = mol.GetConformer()

    for atom_idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atom_idx)
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() == "PJE":
            pje_atoms.append(atom_idx)

    if not pje_atoms:
        return None

    # Build PJE molecule
    pje_mol = Chem.RWMol()
    atom_map = {}

    for old_idx in pje_atoms:
        atom = mol.GetAtomWithIdx(old_idx)
        new_atom = Chem.Atom(atom.GetSymbol())
        new_idx = pje_mol.AddAtom(new_atom)
        atom_map[old_idx] = new_idx

    # Add bonds
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in pje_atoms and end in pje_atoms:
            pje_mol.AddBond(atom_map[begin], atom_map[end], bond.GetBondType())

    pje_mol = pje_mol.GetMol()

    # Add conformer
    new_conf = Chem.Conformer(len(pje_atoms))
    for old_idx, new_idx in atom_map.items():
        pos = conf.GetAtomPosition(old_idx)
        new_conf.SetAtomPosition(new_idx, pos)
    pje_mol.AddConformer(new_conf)

    return pje_mol


def find_bonded_carbon_near_sg(pje_mol, sg_pos, cutoff=2.5):
    """
    Find the carbon atom in PJE that is bonded to CYS SG.

    Returns (carbon_idx, neighbor_idx) where neighbor is likely the other carbon in C=C
    """
    conf = pje_mol.GetConformer()

    # Find carbon closest to SG
    best_idx = None
    best_dist = 999.0

    for atom_idx in range(pje_mol.GetNumAtoms()):
        atom = pje_mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            continue

        pos = conf.GetAtomPosition(atom_idx)
        pos_arr = np.array([pos.x, pos.y, pos.z])
        dist = np.linalg.norm(pos_arr - sg_pos)

        if dist < best_dist:
            best_dist = dist
            best_idx = atom_idx

    if best_idx is None or best_dist > cutoff:
        return None, None

    # Find neighbor carbon (likely the other end of original C=C)
    atom = pje_mol.GetAtomWithIdx(best_idx)
    for bond in atom.GetBonds():
        neighbor_idx = bond.GetOtherAtomIdx(best_idx)
        neighbor = pje_mol.GetAtomWithIdx(neighbor_idx)
        if neighbor.GetSymbol() == 'C':
            # Check if this carbon has reasonable geometry for vinyl
            return best_idx, neighbor_idx

    return best_idx, None


def reconstruct_vinyl_from_adduct(pje_mol, c1_idx, c2_idx):
    """
    Reconstruct C=C from S-C-C adduct.

    Simply changes C-C single bond to C=C double bond.
    """
    if c1_idx is None or c2_idx is None:
        print("Cannot reconstruct: missing carbon indices")
        return None

    editable = Chem.RWMol(pje_mol)

    # Change C-C bond to C=C
    bond = editable.GetBondBetweenAtoms(c1_idx, c2_idx)
    if bond:
        editable.RemoveBond(c1_idx, c2_idx)
        editable.AddBond(c1_idx, c2_idx, Chem.BondType.DOUBLE)

    result = editable.GetMol()

    try:
        Chem.SanitizeMol(result)
    except Exception as e:
        print(f"Warning: Sanitization failed: {e}")

    return result


def main():
    parser = argparse.ArgumentParser(description="Reconstruct ligand from covalent adduct")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-o", "--output", help="Output SDF file for reconstructed ligand")

    args = parser.parse_args()

    print(f"Loading PDB: {args.input}")

    # Extract PJE
    pje_mol = extract_pje_from_pdb(args.input)
    if pje_mol is None:
        print("ERROR: Could not extract PJE from PDB")
        return

    print(f"Extracted PJE: {pje_mol.GetNumAtoms()} atoms")

    # CYS145 SG position (from 6lu7.pdb)
    sg_pos = np.array([-14.043, 17.445, 66.228])

    # Find carbons involved in covalent bond
    c1_idx, c2_idx = find_bonded_carbon_near_sg(pje_mol, sg_pos)

    if c1_idx is not None:
        print(f"Found carbon bonded to SG: atom {c1_idx}")
        if c2_idx is not None:
            print(f"Found neighbor carbon: atom {c2_idx}")
            print("Reconstructing C=C double bond...")

            # Reconstruct
            reconstructed = reconstruct_vinyl_from_adduct(pje_mol, c1_idx, c2_idx)

            if reconstructed:
                smiles = Chem.MolToSmiles(reconstructed)
                print(f"Reconstructed SMILES: {smiles}")

                if args.output:
                    writer = Chem.SDWriter(args.output)
                    writer.write(reconstructed)
                    writer.close()
                    print(f"Saved to: {args.output}")
            else:
                print("ERROR: Reconstruction failed")
        else:
            print("ERROR: Could not find neighbor carbon for vinyl reconstruction")
    else:
        print("ERROR: Could not find carbon bonded to SG")

    # Show current structure
    try:
        Chem.SanitizeMol(pje_mol)
        adduct_smiles = Chem.MolToSmiles(pje_mol)
        print(f"\nAdduct SMILES (with S bonded): {adduct_smiles}")
    except:
        print("Could not generate SMILES for adduct")


if __name__ == "__main__":
    main()
