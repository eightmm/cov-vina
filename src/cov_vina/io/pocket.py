"""Cached loading utilities for repeatedly reused protein pockets."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Callable, Dict, Tuple, Optional

import torch
import numpy as np
from rdkit import Chem


@dataclass(frozen=True)
class PocketBundle:
    mol: Chem.Mol
    coords: torch.Tensor
    features: dict


_POCKET_CACHE: Dict[Tuple[str, int, int, str], PocketBundle] = {}


def _cache_key(protein_pdb: str, device: torch.device) -> Tuple[str, int, int, str]:
    abs_path = os.path.abspath(protein_pdb)
    stat = os.stat(abs_path)
    return abs_path, stat.st_mtime_ns, stat.st_size, str(device)


def load_pocket_bundle(
    protein_pdb: str,
    device: torch.device,
    feature_builder: Callable[[Chem.Mol], dict],
) -> PocketBundle:
    """
    Load and cache pocket molecule data for repeated runs on the same receptor.

    The cache is invalidated when the pocket file path, modified time, file size,
    or target device changes.
    """
    key = _cache_key(protein_pdb, device)
    cached = _POCKET_CACHE.get(key)
    if cached is not None:
        return cached

    pocket_mol = Chem.MolFromPDBFile(protein_pdb, sanitize=False, removeHs=True)
    if pocket_mol is None:
        raise ValueError(f"Failed to load protein pocket from {protein_pdb}")

    pocket_coords = torch.tensor(
        pocket_mol.GetConformer().GetPositions(),
        dtype=torch.float32,
        device=device,
    )
    pocket_features = feature_builder(pocket_mol)

    bundle = PocketBundle(
        mol=pocket_mol,
        coords=pocket_coords,
        features=pocket_features,
    )
    _POCKET_CACHE[key] = bundle
    return bundle


def clear_pocket_cache() -> None:
    """Clear all cached pocket bundles."""
    _POCKET_CACHE.clear()


def extract_pocket_around_residue(
    protein_mol: Chem.Mol,
    residue_spec: str,
    cutoff: float = 12.0
) -> Chem.Mol:
    """
    Extract pocket atoms within cutoff distance of specified residue.

    Args:
        protein_mol: Full protein RDKit Mol object
        residue_spec: Residue specifier, e.g. "CYS145" or "CYS145:A"
        cutoff: Distance cutoff in Angstroms (default: 12.0)
                - 12Å is standard for Vina docking (covers ~10Å effective range)
                - Smaller values (10Å) for small molecules
                - Larger values (15-20Å) for large ligands/peptides

    Returns:
        Pocket Mol with atoms within cutoff of the residue

    Raises:
        ValueError: If residue not found or protein has no conformer
    """
    # Parse residue spec
    parts = residue_spec.split(':')
    if len(parts) == 2:
        res_name_num, chain_id = parts
    else:
        res_name_num = parts[0]
        chain_id = None

    # Extract residue name and number
    import re
    match = re.match(r'([A-Z]+)(\d+)', res_name_num)
    if not match:
        raise ValueError(f"Invalid residue spec: {residue_spec}. Expected format: 'CYS145' or 'CYS145:A'")

    res_name = match.group(1)
    res_num = int(match.group(2))

    # Get protein conformer
    if protein_mol.GetNumConformers() == 0:
        raise ValueError("Protein molecule has no conformer")
    conf = protein_mol.GetConformer()

    # Find residue atoms
    residue_atom_indices = []
    for atom in protein_mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is None:
            continue

        atom_res_name = pdb_info.GetResidueName().strip()
        atom_res_num = pdb_info.GetResidueNumber()
        atom_chain = pdb_info.GetChainId().strip()

        # Match residue
        if atom_res_name == res_name and atom_res_num == res_num:
            if chain_id is None or atom_chain == chain_id:
                residue_atom_indices.append(atom.GetIdx())

    if not residue_atom_indices:
        raise ValueError(
            f"Residue {residue_spec} not found in protein. "
            f"Available residues can be checked with: "
            f"'SELECT * FROM protein_mol WHERE residue_name=\"{res_name}\"'"
        )

    # Get residue center (for distance calculation) - vectorized
    residue_coords = np.array([
        [conf.GetAtomPosition(idx).x, conf.GetAtomPosition(idx).y, conf.GetAtomPosition(idx).z]
        for idx in residue_atom_indices
    ])
    residue_center = residue_coords.mean(axis=0)

    # Get all protein coordinates - vectorized
    num_atoms = protein_mol.GetNumAtoms()
    all_coords = np.array([
        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
        for i in range(num_atoms)
    ])

    # Calculate all distances at once - vectorized
    distances = np.linalg.norm(all_coords - residue_center, axis=1)

    # Find atoms within cutoff - vectorized boolean indexing
    pocket_atom_indices = np.where(distances <= cutoff)[0].tolist()

    if not pocket_atom_indices:
        raise ValueError(f"No atoms found within {cutoff}Å of residue {residue_spec}")

    # Create pocket molecule
    pocket_mol = Chem.RWMol()
    old_to_new_idx = {}

    # Add atoms
    for old_idx in pocket_atom_indices:
        atom = protein_mol.GetAtomWithIdx(old_idx)
        new_idx = pocket_mol.AddAtom(atom)
        old_to_new_idx[old_idx] = new_idx

    # Add bonds
    for bond in protein_mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx in old_to_new_idx and end_idx in old_to_new_idx:
            pocket_mol.AddBond(
                old_to_new_idx[begin_idx],
                old_to_new_idx[end_idx],
                bond.GetBondType()
            )

    # Add conformer
    pocket_conf = Chem.Conformer(len(pocket_atom_indices))
    for old_idx, new_idx in old_to_new_idx.items():
        pos = conf.GetAtomPosition(old_idx)
        pocket_conf.SetAtomPosition(new_idx, pos)

    pocket_mol = pocket_mol.GetMol()
    pocket_mol.AddConformer(pocket_conf, assignId=True)

    return pocket_mol
