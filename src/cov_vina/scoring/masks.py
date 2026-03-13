"""Intramolecular interaction mask computation."""

import torch
from rdkit import Chem


def compute_intramolecular_mask(mol: Chem.Mol, device: torch.device,
                                exclude_atom_indices: set = None) -> torch.Tensor:
    """
    Computes a boolean mask [N, N] for intramolecular Vina scoring.
    True means the interaction between atom i and j SHOULD be calculated.
    False means they are too close (1-2, 1-3 pairs, or in the same ring) and and should be ignored.

    OPTIMIZED: Vectorized operations instead of nested loops.
    """
    num_atoms = mol.GetNumAtoms()

    # 1. Get adjacency matrix (1-2 pairs = 1 hop) - numpy array
    adj = Chem.GetDistanceMatrix(mol)

    # 2. Vectorized: mask starts True where distance > 2 (ignore 1-2 and 1-3 pairs)
    mask = torch.from_numpy(adj > 2).to(device)

    # 3. Ignore atoms in the same ring - vectorized with broadcasting
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_tensor = torch.tensor(list(ring), dtype=torch.long, device=device)
        # Create a meshgrid for the ring atoms
        mask[ring_tensor[:, None], ring_tensor] = False

    if exclude_atom_indices:
        for idx in exclude_atom_indices:
            mask[idx, :] = False
            mask[:, idx] = False

    return mask
