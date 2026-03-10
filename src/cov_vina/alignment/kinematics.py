"""Forward kinematics for ligand torsion angle optimization."""

import torch
import torch.nn as nn
from torch import sin, cos
from rdkit import Chem
from typing import Dict, List, Tuple


def get_batched_rotation_matrix(axis: torch.Tensor, theta: torch.Tensor) -> torch.Tensor:
    """
    Computes Rodrigues' rotation matrices for a batch of axes/angles.
    axis: [B, 3]
    theta: [B]
    Returns: [B, 3, 3]
    """
    axis = axis / torch.norm(axis, dim=1, keepdim=True).clamp_min(1e-12)
    a, b, c = axis[:, 0], axis[:, 1], axis[:, 2]

    sin_t = sin(theta)
    cos_t = cos(theta)
    one_minus_cos_t = 1 - cos_t

    return torch.stack([
        cos_t + a * a * one_minus_cos_t,
        a * b * one_minus_cos_t - c * sin_t,
        a * c * one_minus_cos_t + b * sin_t,
        a * b * one_minus_cos_t + c * sin_t,
        cos_t + b * b * one_minus_cos_t,
        b * c * one_minus_cos_t - a * sin_t,
        a * c * one_minus_cos_t - b * sin_t,
        b * c * one_minus_cos_t + a * sin_t,
        cos_t + c * c * one_minus_cos_t,
    ], dim=1).reshape(-1, 3, 3)


def _compute_descendants(tree: Dict[int, List[Tuple[int, int, int]]],
                         frames: List[List[int]],
                         frame_idx: int) -> List[int]:
    descendants = []
    for child_frame, _, _ in tree[frame_idx]:
        descendants.extend(frames[child_frame])
        descendants.extend(_compute_descendants(tree, frames, child_frame))
    return descendants


def _build_kinematic_topology(mol: Chem.Mol, ref_indices: List[int], freeze_anchor: bool):
    num_atoms = mol.GetNumAtoms()
    rot_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    rot_bonds = mol.GetSubstructMatches(rot_smarts)

    ref_set = set(ref_indices)
    rot_bond_pairs = []
    if freeze_anchor:
        for pair in rot_bonds:
            p_set = set(pair)
            if not p_set.issubset(ref_set):
                rot_bond_pairs.append(p_set)
    else:
        rot_bond_pairs = [set(pair) for pair in rot_bonds]

    adj_list = {i: [] for i in range(num_atoms)}
    for bond in mol.GetBonds():
        u, v = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if {u, v} not in rot_bond_pairs:
            adj_list[u].append(v)
            adj_list[v].append(u)

    visited = set()
    frames = []
    atom_to_frame = {}
    for i in range(num_atoms):
        if i not in visited:
            frame_atoms = []
            queue = [i]
            visited.add(i)
            while queue:
                curr = queue.pop(0)
                frame_atoms.append(curr)
                atom_to_frame[curr] = len(frames)
                for neighbor in adj_list[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            frames.append(frame_atoms)

    root_frame_idx = -1
    max_overlap = -1
    for f_idx, frame_atoms in enumerate(frames):
        overlap = len(set(frame_atoms).intersection(ref_set))
        if overlap > max_overlap:
            max_overlap = overlap
            root_frame_idx = f_idx

    tree = {i: [] for i in range(len(frames))}
    q = [root_frame_idx]
    visited_frames = {root_frame_idx}
    kinematic_edges = []

    while q:
        curr_frame = q.pop(0)
        for u in frames[curr_frame]:
            atom = mol.GetAtomWithIdx(u)
            for neighbor in atom.GetNeighbors():
                v = neighbor.GetIdx()
                neighbor_frame = atom_to_frame[v]
                if neighbor_frame != curr_frame and neighbor_frame not in visited_frames:
                    tree[curr_frame].append((neighbor_frame, u, v))
                    kinematic_edges.append((u, v, neighbor_frame))
                    visited_frames.add(neighbor_frame)
                    q.append(neighbor_frame)

    parent_atoms = [edge[0] for edge in kinematic_edges]
    child_atoms = [edge[1] for edge in kinematic_edges]
    child_frames = [edge[2] for edge in kinematic_edges]
    frame_descendants_cache = {
        frame_idx: _compute_descendants(tree, frames, frame_idx)
        for frame_idx in range(len(frames))
    }
    atoms_to_rotate = []
    for frame_child in child_frames:
        atoms_to_rotate.append(frames[frame_child] + frame_descendants_cache[frame_child])

    return {
        "num_atoms": num_atoms,
        "frames": frames,
        "tree": tree,
        "parent_atoms": parent_atoms,
        "child_atoms": child_atoms,
        "child_frames": child_frames,
        "atoms_to_rotate": atoms_to_rotate,
        "num_torsions": len(kinematic_edges),
    }


class LigandKinematics(nn.Module):
    def __init__(self, mol: Chem.Mol, ref_indices: List[int], init_coords: torch.Tensor, device: torch.device, freeze_anchor: bool = True):
        """
        Forward kinematics for ligand torsion angle optimization.

        Automatically handles both single and batched inputs:
        - init_coords: [N, 3] -> single pose optimization
        - init_coords: [B, N, 3] -> batched optimization

        Args:
            mol: RDKit molecule
            ref_indices: List of atom indices that should remain fixed (anchor atoms)
            init_coords: [N, 3] or [B, N, 3] initial coordinates
            device: torch device
            freeze_anchor: If True, rotatable bonds inside anchor are frozen (rigid)
        """
        super().__init__()
        self.device = device
        topology = _build_kinematic_topology(mol, ref_indices, freeze_anchor)
        self.num_atoms = topology["num_atoms"]
        self.tree = topology["tree"]
        self.num_torsions = topology["num_torsions"]
        self.frames = topology["frames"]
        self.parent_atoms = topology["parent_atoms"]
        self.child_atoms = topology["child_atoms"]
        self.child_frames = topology["child_frames"]
        self.atoms_to_rotate_tensors = [
            torch.tensor(atoms, dtype=torch.long, device=self.device)
            for atoms in topology["atoms_to_rotate"]
        ]

        # Detect if batched input
        init_coords = init_coords.to(self.device)
        if init_coords.dim() == 2:
            # Single pose: [N, 3] -> add batch dimension
            self.is_batched = False
            self.base_coords = init_coords.unsqueeze(0)  # [1, N, 3]
            self.batch_size = 1
        elif init_coords.dim() == 3:
            # Batched: [B, N, 3]
            self.is_batched = True
            self.base_coords = init_coords.clone()
            self.batch_size = init_coords.shape[0]
        else:
            raise ValueError(f"init_coords must be [N, 3] or [B, N, 3], got shape {init_coords.shape}")

        # Always use batched parameter shape [B, num_torsions]
        self.thetas = nn.Parameter(torch.zeros(self.batch_size, self.num_torsions, device=self.device))

    def forward(self) -> torch.Tensor:
        """
        Compute new coordinates from torsion angles.

        Returns:
            [N, 3] if initialized with single pose
            [B, N, 3] if initialized with batched poses
        """
        coords = self.base_coords.clone()  # [B, N, 3]

        for i in range(self.num_torsions):
            idx_parent = self.parent_atoms[i]
            idx_child = self.child_atoms[i]
            theta = self.thetas[:, i]  # [B]

            origin = coords[:, idx_parent, :]  # [B, 3]
            axis = coords[:, idx_child, :] - origin  # [B, 3]
            R = get_batched_rotation_matrix(axis, theta)  # [B, 3, 3]

            atoms_tensor = self.atoms_to_rotate_tensors[i]
            if len(atoms_tensor) == 0:
                continue

            vecs = coords[:, atoms_tensor, :] - origin.unsqueeze(1)  # [B, M, 3]
            new_vecs = torch.matmul(vecs, R.transpose(1, 2))  # [B, M, 3]

            next_coords = coords.clone()
            next_coords[:, atoms_tensor, :] = new_vecs + origin.unsqueeze(1)
            coords = next_coords

        # Remove batch dimension if input was single pose
        if not self.is_batched:
            coords = coords.squeeze(0)  # [N, 3]

        return coords
