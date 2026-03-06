"""Torsion angle optimization using gradient descent."""

import torch
from rdkit import Chem
from typing import List

from ..alignment.kinematics import LigandKinematics


def optimize_torsions_vina(mol: Chem.Mol,
                           ref_indices: List[int],
                           init_coords: torch.Tensor,
                           pocket_coords: torch.Tensor,
                           query_features: dict,
                           pocket_features: dict,
                           device: torch.device,
                           num_steps: int = 100,
                           lr: float = 0.1,
                           freeze_mcs: bool = True,
                           num_rotatable_bonds: int = None,
                           weight_preset: str = 'vina',
                           batch_size: int = 8,
                           optimizer: str = 'adam',
                           early_stopping: bool = True,
                           patience: int = 30,
                           min_delta: float = 1e-5) -> torch.Tensor:
    """
    Optimizes torsion angles to minimize Vina score.

    Handles both single pose and batched optimization automatically.

    Args:
        mol: RDKit molecule
        ref_indices: MCS atom indices to freeze
        init_coords: Initial coordinates [N_atoms, 3] or [N_poses, N_atoms, 3]
        pocket_coords: Pocket coordinates [N_pocket_atoms, 3]
        query_features: Query molecular features
        pocket_features: Pocket molecular features
        device: torch device
        num_steps: Number of optimization steps (default: 100)
        lr: Learning rate (default: 0.1 for adam/adamw, 1.0 for lbfgs)
        freeze_mcs: Whether to freeze MCS atoms (default: True)
        num_rotatable_bonds: Number of rotatable bonds for torsion penalty
        weight_preset: Vina weight preset ('vina', 'vina_lp', 'vinardo')
        batch_size: Batch size for multi-pose optimization (default: 8)
        optimizer: Optimizer type ('adam', 'adamw', 'lbfgs') (default: 'adam')
        early_stopping: Enable early stopping (default: True)
        patience: Steps without improvement before stopping (default: 30)
        min_delta: Minimum improvement to reset patience (default: 1e-5)

    Returns:
        optimized_coords: [N_atoms, 3] or [N_poses, N_atoms, 3] depending on input
    """
    from ..scoring import vina_scoring
    from ..scoring.masks import compute_intramolecular_mask

    # Auto-detect single vs batched input
    if init_coords.ndim == 2:
        # Single pose: [N_atoms, 3] → add batch dim
        init_coords = init_coords.unsqueeze(0)
        single_pose = True
    else:
        single_pose = False

    n_poses = init_coords.shape[0]
    n_atoms = init_coords.shape[1]

    # Precompute intramolecular mask once
    intra_mask = compute_intramolecular_mask(mol, device)

    # Check if molecule has rotatable bonds
    test_model = LigandKinematics(mol, ref_indices, init_coords[0], device, freeze_mcs=freeze_mcs)
    if test_model.num_torsions == 0:
        if n_poses > 1:
            print("No rotatable bonds found - returning initial coordinates")
        return init_coords[0] if single_pose else init_coords.clone()

    optimized_coords = torch.zeros_like(init_coords)

    # Process in batches
    n_batches = (n_poses + batch_size - 1) // batch_size

    if n_poses > 1:
        print(f"Optimizing {n_poses} poses in {n_batches} batches (batch_size={batch_size})...")

    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_poses)

        if n_poses > 1:
            print(f"  Batch {batch_idx + 1}/{n_batches}: Optimizing poses {start_idx}-{end_idx-1}...")

        # Create separate models for each pose in batch
        models = []
        optimizers_list = []
        for i in range(start_idx, end_idx):
            model = LigandKinematics(mol, ref_indices, init_coords[i], device, freeze_mcs=freeze_mcs)

            # Select optimizer
            if optimizer.lower() == 'adam':
                opt = torch.optim.Adam(model.parameters(), lr=lr)
            elif optimizer.lower() == 'adamw':
                opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
            elif optimizer.lower() == 'lbfgs':
                opt = torch.optim.LBFGS(model.parameters(), lr=lr, max_iter=20,
                                       history_size=10, line_search_fn='strong_wolfe')
            else:
                raise ValueError(f"Unknown optimizer: {optimizer}. Choose 'adam', 'adamw', or 'lbfgs'")

            models.append(model)
            optimizers_list.append(opt)

        # Optimization loop
        if optimizer.lower() == 'lbfgs':
            # LBFGS requires closure function with early stopping
            batch_len = end_idx - start_idx
            converged = [False] * batch_len
            patience_counts = [0] * batch_len
            best_losses = [float('inf')] * batch_len

            for step in range(num_steps):
                active_count = 0

                for i, (model, opt) in enumerate(zip(models, optimizers_list)):
                    if early_stopping and converged[i]:
                        continue  # Skip converged poses

                    active_count += 1

                    def closure():
                        opt.zero_grad()
                        coords = model()
                        loss = vina_scoring(coords.unsqueeze(0), pocket_coords, query_features,
                                          pocket_features, num_rotatable_bonds, weight_preset,
                                          intramolecular_mask=intra_mask)
                        loss = loss.sum()
                        loss.backward()
                        return loss

                    loss = opt.step(closure)

                    # Early stopping check
                    if early_stopping:
                        loss_val = loss.item() if isinstance(loss, torch.Tensor) else loss
                        if loss_val < best_losses[i] - min_delta:
                            best_losses[i] = loss_val
                            patience_counts[i] = 0
                        else:
                            patience_counts[i] += 1

                        if patience_counts[i] >= patience:
                            converged[i] = True

                # Check if all poses converged
                if early_stopping and active_count == 0:
                    if n_poses > 1:
                        print(f"    All {batch_len} poses converged at step {step + 1}")
                    break
        else:
            # Adam/AdamW standard loop with early stopping
            batch_len = end_idx - start_idx
            converged = [False] * batch_len
            patience_counts = [0] * batch_len
            best_losses = [float('inf')] * batch_len

            for step in range(num_steps):
                active_count = 0

                for i, (model, opt) in enumerate(zip(models, optimizers_list)):
                    if early_stopping and converged[i]:
                        continue  # Skip converged poses

                    active_count += 1
                    opt.zero_grad()
                    coords = model()
                    loss = vina_scoring(coords.unsqueeze(0), pocket_coords, query_features, pocket_features,
                                       num_rotatable_bonds, weight_preset, intramolecular_mask=intra_mask)
                    loss = loss.sum()
                    loss.backward()
                    # Gradient clipping for stability
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    opt.step()

                    # Early stopping check
                    if early_stopping:
                        loss_val = loss.item()
                        if loss_val < best_losses[i] - min_delta:
                            best_losses[i] = loss_val
                            patience_counts[i] = 0
                        else:
                            patience_counts[i] += 1

                        if patience_counts[i] >= patience:
                            converged[i] = True

                # Check if all poses converged
                if early_stopping and active_count == 0:
                    if n_poses > 1:
                        print(f"    All {batch_len} poses converged at step {step + 1}")
                    break

        # Extract final coordinates
        with torch.no_grad():
            for i, model in enumerate(models):
                optimized_coords[start_idx + i] = model()

    if n_poses > 1:
        print(f"✓ Optimization complete!")

    # Return single pose without batch dim if input was single
    return optimized_coords[0] if single_pose else optimized_coords
