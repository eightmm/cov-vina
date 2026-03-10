"""
Create protein-ligand covalent adduct by removing leaving group and forming bond.
"""

from typing import Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from .anchor import WarheadHit, AnchorPoint


# Leaving group patterns for each warhead type
# Maps warhead_type -> list of atom indices to remove (relative to matched_atoms)
LEAVING_GROUPS = {
    # Alpha-halo carbonyls: remove halogen
    "chloroacetamide": [0],   # Remove Cl
    "bromoacetamide": [0],    # Remove Br
    "iodoacetamide": [0],     # Remove I
    "fluoroacetamide": [0],   # Remove F
    "chlorofluoroacetamide": [0, 1],  # Remove Cl, F

    # Michael acceptors: no leaving group (addition reaction)
    "acrylamide": [],
    "acrylic_acid": [],
    "acrylate": [],
    "enone": [],
    "vinyl_sulfonamide": [],
    "vinyl_sulfone": [],
    "maleimide": [],

    # Epoxide/aziridine/thiirane: ring opens, no atom removed
    "epoxide": [],
    "aziridine": [],
    "thiirane": [],

    # Nitriles: no leaving group (addition to C≡N)
    "aryl_nitrile": [],
    "alkyl_nitrile": [],
    "propiolamide": [],
    "propargylamide": [],
    "cyanoacrylamide": [],

    # Others: no leaving group
    "disulfide": [],
    "aldehyde": [],
    "isothiocyanate": [],
    "boronic_acid": [],
    "phosphonate": [],

    # Activated esters: leaving group removed (NHS, TFE)
    "nhs_ester": [],  # NHS stays attached, complex mechanism
    "tfe_ester": [],  # TFE group leaves, but complex

    # Acyl fluoride & sulfonyl fluoride: remove F
    "acyl_fluoride": [0],
    "sulfonyl_fluoride": [0],
}


def create_adduct_template(
    ligand_mol: Chem.Mol,
    warhead: WarheadHit,
    anchor: AnchorPoint,
) -> tuple[Chem.Mol, Optional[int], Optional[int], int]:
    """
    Create adduct template molecule (topology only, no conformers).

    This is called BEFORE conformer generation to:
    1. Remove leaving group from ligand
    2. Add protein anchor atoms (Cβ and S)
    3. Form covalent bonds

    The resulting molecule will be used for conformer generation with CB-S coordmap.

    Args:
        ligand_mol: Query ligand with warhead (NO conformers required)
        warhead: Detected warhead information
        anchor: Protein anchor point (contains both Cβ and S coordinates)

    Returns:
        adduct_mol: Modified ligand with leaving group removed and anchor atoms added
        cb_atom_idx: Index of added Cβ atom, or None if not added
        s_atom_idx: Index of added S atom, or None if not added
        new_reactive_idx: Updated reactive atom index after leaving group removal
    """
    # Step 1: Identify leaving group atoms to remove
    leaving_group_pattern = LEAVING_GROUPS.get(warhead.warhead_type, [])

    # Step 2: Get leaving group atom indices in the full molecule
    matched_atom_list = list(warhead.matched_atoms) if warhead.matched_atoms else []
    leaving_atoms_to_remove = [matched_atom_list[i] for i in leaving_group_pattern] if leaving_group_pattern and matched_atom_list else []

    # Step 3: Create editable copy
    editable_mol = Chem.RWMol(ligand_mol)

    # Step 4: Handle Michael addition (no leaving group)
    if not leaving_group_pattern:
        # Find the double bond adjacent to reactive atom (beta carbon)
        reactive_atom = editable_mol.GetAtomWithIdx(warhead.reactive_atom_idx)
        double_bond_neighbor = None

        for bond in reactive_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other_idx = bond.GetOtherAtomIdx(warhead.reactive_atom_idx)
                double_bond_neighbor = other_idx
                break

        if double_bond_neighbor is None:
            raise ValueError("No double bond found for Michael addition mechanism")

        # Change double bond (C=C) to single bond (C-C)
        editable_mol.RemoveBond(warhead.reactive_atom_idx, double_bond_neighbor)
        editable_mol.AddBond(warhead.reactive_atom_idx, double_bond_neighbor, Chem.BondType.SINGLE)
    else:
        # Step 5: Remove leaving group atoms (sort descending to avoid index shift)
        for atom_idx in sorted(leaving_atoms_to_remove, reverse=True):
            editable_mol.RemoveAtom(atom_idx)

    # Step 6: Add protein anchor atoms (Cβ and nucleophilic atom)
    cb_idx = None
    if anchor.cb_coord is not None:
        cb_idx = editable_mol.AddAtom(Chem.Atom(6))  # Carbon (Cβ)

    # Determine nucleophilic atom type from residue
    # Map atom name (from PDB) to atomic number
    nucleophile_atom_map = {
        "SG": 16,   # CYS: Sulfur (thiol)
        "OG": 8,    # SER: Oxygen (hydroxyl)
        "OG1": 8,   # THR: Oxygen (hydroxyl)
        "OH": 8,    # TYR: Oxygen (phenolic)
        "NZ": 7,    # LYS: Nitrogen (amine)
        "NE2": 7,   # HIS: Nitrogen (imidazole)
    }
    nuc_atomic_num = nucleophile_atom_map.get(anchor.atom_name, 16)  # Default to S
    nuc_idx = editable_mol.AddAtom(Chem.Atom(nuc_atomic_num))

    # Step 7: Form bonds: Cβ-Nuc-C(reactive)
    if cb_idx is not None:
        editable_mol.AddBond(cb_idx, nuc_idx, Chem.BondType.SINGLE)  # Cβ-Nuc bond

    # Calculate new reactive atom index after leaving group removal
    shift = sum(1 for idx in leaving_atoms_to_remove if idx < warhead.reactive_atom_idx)
    new_reactive_idx = warhead.reactive_atom_idx - shift

    editable_mol.AddBond(new_reactive_idx, nuc_idx, Chem.BondType.SINGLE)  # Nuc-C bond

    # Step 8: Get modified molecule
    adduct_mol = editable_mol.GetMol()

    # Step 9: Sanitize
    try:
        Chem.SanitizeMol(adduct_mol)
    except Exception as e:
        print(f"Warning: Sanitization failed after creating adduct template: {e}")

    return adduct_mol, cb_idx, nuc_idx, new_reactive_idx


def create_covalent_adduct(
    ligand_mol: Chem.Mol,
    warhead: WarheadHit,
    anchor: AnchorPoint,
    protein_mol: Optional[Chem.Mol] = None,
    add_anchor_atom: bool = True
) -> tuple[Chem.Mol, Optional[int], Optional[int]]:
    """
    DEPRECATED: Use create_adduct_template() before conformer generation instead.

    Create protein-ligand covalent adduct by:
    1. Removing leaving group from ligand
    2. Adding protein anchor atoms (Cβ and S) and forming covalent bonds

    Args:
        ligand_mol: Query ligand with warhead (WITH conformer)
        warhead: Detected warhead information
        anchor: Protein anchor point (contains both Cβ and S coordinates)
        protein_mol: Protein molecule (currently unused, for future extensions)
        add_anchor_atom: If True, add protein Cβ and S atoms to create true adduct

    Returns:
        adduct_mol: Modified ligand with leaving group removed and anchor atoms added
        cb_atom_idx: Index of added Cβ atom, or None if not added
        s_atom_idx: Index of added S atom, or None if not added
    """

    # Step 1: Identify leaving group atoms to remove
    leaving_group_pattern = LEAVING_GROUPS.get(warhead.warhead_type, [])

    if not leaving_group_pattern:
        # No leaving group - Michael addition mechanism
        # For Michael acceptors: C=C-C(=O) + S⁻ → C-C(S)-C(=O)
        if add_anchor_atom:
            editable_mol = Chem.RWMol(ligand_mol)
            has_conf = ligand_mol.GetNumConformers() > 0
            if has_conf:
                old_conf = ligand_mol.GetConformer()

            # Find the double bond adjacent to reactive atom (beta carbon)
            reactive_atom = editable_mol.GetAtomWithIdx(warhead.reactive_atom_idx)
            double_bond_neighbor = None

            for bond in reactive_atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other_idx = bond.GetOtherAtomIdx(warhead.reactive_atom_idx)
                    double_bond_neighbor = other_idx
                    break

            if double_bond_neighbor is None:
                print("Warning: No double bond found for Michael addition, returning copy")
                adduct_mol = Chem.Mol(ligand_mol)
                for conf_id in range(ligand_mol.GetNumConformers()):
                    conf = ligand_mol.GetConformer(conf_id)
                    new_conf = Chem.Conformer(conf)
                    adduct_mol.AddConformer(new_conf, assignId=True)
                return adduct_mol, None, None

            # Change double bond (C=C) to single bond (C-C)
            editable_mol.RemoveBond(warhead.reactive_atom_idx, double_bond_neighbor)
            editable_mol.AddBond(warhead.reactive_atom_idx, double_bond_neighbor, Chem.BondType.SINGLE)

            # Add Cβ and S atoms from protein (Cβ-S-C chain)
            # Add Cβ first (this will be the anchor - frozen)
            cb_idx = None
            if anchor.cb_coord is not None:
                cb_idx = editable_mol.AddAtom(Chem.Atom(6))  # Carbon (Cβ)

            # Add S atom
            s_idx = editable_mol.AddAtom(Chem.Atom(16))  # Sulfur

            # Form bonds: Cβ-S-C(reactive)
            if cb_idx is not None:
                editable_mol.AddBond(cb_idx, s_idx, Chem.BondType.SINGLE)  # Cβ-S bond
            editable_mol.AddBond(s_idx, warhead.reactive_atom_idx, Chem.BondType.SINGLE)  # S-C bond

            # Get molecule
            adduct_mol = editable_mol.GetMol()

            # Sanitize
            try:
                Chem.SanitizeMol(adduct_mol)
            except Exception as e:
                print(f"Warning: Sanitization failed after Michael addition: {e}")

            # Update conformer
            if has_conf:
                new_conf = Chem.Conformer(adduct_mol.GetNumAtoms())

                # Copy all original atom positions
                for idx in range(ligand_mol.GetNumAtoms()):
                    pos = old_conf.GetAtomPosition(idx)
                    new_conf.SetAtomPosition(idx, pos)

                # Add Cβ atom position if present
                if cb_idx is not None and anchor.cb_coord is not None:
                    new_conf.SetAtomPosition(cb_idx,
                                            Point3D(float(anchor.cb_coord[0]),
                                                   float(anchor.cb_coord[1]),
                                                   float(anchor.cb_coord[2])))

                # Add S atom position
                new_conf.SetAtomPosition(s_idx,
                                        Point3D(float(anchor.coord[0]),
                                               float(anchor.coord[1]),
                                               float(anchor.coord[2])))

                adduct_mol.RemoveAllConformers()
                adduct_mol.AddConformer(new_conf)

            return adduct_mol, cb_idx, s_idx
        else:
            # Don't add anchor atom - just return copy
            adduct_mol = Chem.Mol(ligand_mol)
            for conf_id in range(ligand_mol.GetNumConformers()):
                conf = ligand_mol.GetConformer(conf_id)
                new_conf = Chem.Conformer(conf)
                adduct_mol.AddConformer(new_conf, assignId=True)
            return adduct_mol, None, None

    # Step 2: Get leaving group atom indices in the full molecule
    matched_atom_list = list(warhead.matched_atoms)
    leaving_atoms_to_remove = [matched_atom_list[i] for i in leaving_group_pattern] if leaving_group_pattern else []

    # Step 3: Create editable copy
    editable_mol = Chem.RWMol(ligand_mol)
    has_conf = ligand_mol.GetNumConformers() > 0
    if has_conf:
        old_conf = ligand_mol.GetConformer()

    # Step 4: Remove leaving group atoms (sort descending to avoid index shift)
    for atom_idx in sorted(leaving_atoms_to_remove, reverse=True):
        editable_mol.RemoveAtom(atom_idx)

    # Step 5: Add protein anchor atoms (Cβ and S) if requested
    cb_idx = None
    s_idx = None
    if add_anchor_atom:
        # Add Cβ first (this will be the anchor - frozen)
        if anchor.cb_coord is not None:
            cb_idx = editable_mol.AddAtom(Chem.Atom(6))  # Carbon (Cβ)

        # Add sulfur atom (for CYS)
        s_idx = editable_mol.AddAtom(Chem.Atom(16))  # Atomic number 16 = Sulfur

        # Form bonds: Cβ-S-C(reactive)
        if cb_idx is not None:
            editable_mol.AddBond(cb_idx, s_idx, Chem.BondType.SINGLE)  # Cβ-S bond

        # Form bond between reactive atom and sulfur
        # After removing leaving groups, need to recalculate reactive atom index
        shift = sum(1 for idx in leaving_atoms_to_remove if idx < warhead.reactive_atom_idx)
        reactive_idx_after_removal = warhead.reactive_atom_idx - shift

        editable_mol.AddBond(reactive_idx_after_removal, s_idx, Chem.BondType.SINGLE)  # S-C bond

    # Step 6: Get modified molecule
    adduct_mol = editable_mol.GetMol()

    # Step 7: Sanitize
    try:
        Chem.SanitizeMol(adduct_mol)
    except Exception as e:
        print(f"Warning: Sanitization failed after creating adduct: {e}")
        # Try to continue anyway

    # Step 8: Update conformer coordinates
    if has_conf:
        new_conf = Chem.Conformer(adduct_mol.GetNumAtoms())

        # Map old indices to new indices (after removing leaving group)
        new_idx = 0
        for old_idx in range(ligand_mol.GetNumAtoms()):
            if old_idx not in leaving_atoms_to_remove:
                pos = old_conf.GetAtomPosition(old_idx)
                new_conf.SetAtomPosition(new_idx, pos)
                new_idx += 1

        # Add anchor atom positions (Cβ and S)
        if add_anchor_atom:
            # Add Cβ position if present
            if cb_idx is not None and anchor.cb_coord is not None:
                new_conf.SetAtomPosition(cb_idx,
                                         Point3D(float(anchor.cb_coord[0]),
                                                float(anchor.cb_coord[1]),
                                                float(anchor.cb_coord[2])))

            # Add S position
            if s_idx is not None:
                new_conf.SetAtomPosition(s_idx,
                                         Point3D(float(anchor.coord[0]),
                                                float(anchor.coord[1]),
                                                float(anchor.coord[2])))

        adduct_mol.RemoveAllConformers()
        adduct_mol.AddConformer(new_conf)

    return adduct_mol, cb_idx, s_idx


def get_covalent_exclusion_indices(
    ligand_mol: Chem.Mol,
    warhead: WarheadHit,
    n_hop_exclude: int = 2
) -> set[int]:
    """
    Get ligand atom indices that should be excluded from intermolecular scoring
    with the protein (because they're involved in the covalent bond).

    Args:
        ligand_mol: Ligand molecule (AFTER leaving group removal)
        warhead: Warhead hit (with reactive_atom_idx from BEFORE removal)
        n_hop_exclude: How many bonds away from reactive atom to exclude

    Returns:
        Set of atom indices to exclude from ligand-protein scoring
    """

    # Map old reactive_atom_idx to new index after leaving group removal
    leaving_group_pattern = LEAVING_GROUPS.get(warhead.warhead_type, [])
    matched_atom_list = list(warhead.matched_atoms)
    leaving_atoms = [matched_atom_list[i] for i in leaving_group_pattern]

    # Calculate how many atoms before reactive_atom were removed
    shift = sum(1 for idx in leaving_atoms if idx < warhead.reactive_atom_idx)
    new_reactive_idx = warhead.reactive_atom_idx - shift

    # Get all atoms within n_hop_exclude bonds
    exclude_indices = {new_reactive_idx}

    for hop in range(n_hop_exclude):
        current_layer = set(exclude_indices)
        for idx in current_layer:
            if idx < ligand_mol.GetNumAtoms():
                atom = ligand_mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    exclude_indices.add(neighbor.GetIdx())

    return exclude_indices


def get_protein_exclusion_residues(
    anchor: AnchorPoint,
    n_residue_exclude: int = 1
) -> list[str]:
    """
    Get protein residues that should be excluded from intermolecular scoring.

    Args:
        anchor: Anchor point
        n_residue_exclude: How many residues to exclude (0=none, 1=anchor only, 2=anchor+neighbors)

    Returns:
        List of residue identifiers like "CYS145:A"
    """
    # For now, just return the anchor residue
    # Could be extended to exclude neighboring residues
    return [f"{anchor.residue_name}{anchor.residue_num}:{anchor.chain_id}"]


def create_intermolecular_exclusion_mask(
    ligand_mol: Chem.Mol,
    protein_mol: Chem.Mol,
    ligand_exclude_indices: set[int],
    protein_residue_excludes: list[str],
    device
) -> 'torch.Tensor':
    """
    Create exclusion mask for intermolecular scoring in covalent docking.

    Args:
        ligand_mol: Ligand molecule (adduct, after leaving group removal)
        protein_mol: Protein/pocket molecule
        ligand_exclude_indices: Ligand atom indices to exclude (near covalent bond)
        protein_residue_excludes: Protein residue IDs to exclude (e.g., ["CYS145:A"])
        device: torch device

    Returns:
        Boolean tensor (1, num_ligand_atoms, num_protein_atoms)
        True = exclude this pair from scoring
    """
    import torch

    num_ligand_atoms = ligand_mol.GetNumAtoms()
    num_protein_atoms = protein_mol.GetNumAtoms()

    # Initialize mask (all False = include everything)
    mask = torch.zeros(num_ligand_atoms, num_protein_atoms, dtype=torch.bool, device=device)

    # Mark ligand atoms for exclusion
    for lig_idx in ligand_exclude_indices:
        mask[lig_idx, :] = True

    # Mark protein residue atoms for exclusion
    protein_exclude_atom_indices = []
    for atom_idx in range(num_protein_atoms):
        atom = protein_mol.GetAtomWithIdx(atom_idx)
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info:
            res_name = pdb_info.GetResidueName().strip()
            res_num = pdb_info.GetResidueNumber()
            chain_id = pdb_info.GetChainId()
            res_id = f"{res_name}{res_num}:{chain_id}"

            if res_id in protein_residue_excludes:
                protein_exclude_atom_indices.append(atom_idx)

    # Exclude all ligand atoms from interacting with excluded protein atoms
    if protein_exclude_atom_indices:
        mask[:, protein_exclude_atom_indices] = True

    return mask.unsqueeze(0)  # Add batch dimension
