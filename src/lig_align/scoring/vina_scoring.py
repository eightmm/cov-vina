"""AutoDock Vina scoring functions."""

import torch
from .vina_params import VINA_WEIGHTS


def precompute_interaction_matrices(query_features: dict, pocket_features: dict, device: torch.device) -> dict:
    """
    Precompute feature interaction matrices that don't depend on coordinates.
    This should be called once and reused across batches.

    Returns:
        dict with 'is_hydrophobic', 'is_hbond', 'R_ij' (vdW sum)
    """
    num_query_atoms = query_features['vdw'].shape[0]
    num_pocket_atoms = pocket_features['vdw'].shape[0]

    # Hydrophobic match: 1.0 if BOTH atoms are hydrophobic
    q_hydro = query_features['hydro'].unsqueeze(1).expand(-1, num_pocket_atoms)
    p_hydro = pocket_features['hydro'].unsqueeze(0).expand(num_query_atoms, -1)
    is_hydrophobic = (q_hydro * p_hydro)

    # H-Bond match: (Donor-Acceptor) OR (Acceptor-Donor)
    q_hbd = query_features['hbd'].unsqueeze(1).expand(-1, num_pocket_atoms)
    p_hba = pocket_features['hba'].unsqueeze(0).expand(num_query_atoms, -1)

    q_hba = query_features['hba'].unsqueeze(1).expand(-1, num_pocket_atoms)
    p_hbd = pocket_features['hbd'].unsqueeze(0).expand(num_query_atoms, -1)

    is_hbond = ((q_hbd * p_hba) + (q_hba * p_hbd) > 0).float()

    # VdW Sum R_ij
    R_ij = query_features['vdw'].unsqueeze(1) + pocket_features['vdw'].unsqueeze(0)

    return {
        'is_hydrophobic': is_hydrophobic,
        'is_hbond': is_hbond,
        'R_ij': R_ij
    }


def vina_scoring(aligned_query_coords: torch.Tensor,
                 pocket_coords: torch.Tensor,
                 query_features: dict,
                 pocket_features: dict,
                 num_rotatable_bonds: int = None,
                 weight_preset: str = 'vina',
                 intramolecular_mask: torch.Tensor = None,
                 precomputed_matrices: dict = None) -> torch.Tensor:
    """
    Differentiable Vina Scoring logic (Full 5-term version matching DeepRMSD/Vinardo).
    Calculates Intermolecular Energy, and optionally Intramolecular Energy if mask is provided.

    OPTIMIZED: Accepts precomputed_matrices to avoid recomputing feature interactions.
    """
    batch_size = aligned_query_coords.shape[0]
    num_query_atoms = aligned_query_coords.shape[1]
    num_pocket_atoms = pocket_coords.shape[0]

    # Expand pocket coords for batch
    pocket_expanded = pocket_coords.unsqueeze(0).expand(batch_size, -1, -1)

    # 1. Coordinate Distance d_ij
    d_ij = torch.cdist(aligned_query_coords, pocket_expanded) # Shape: (batch_size, num_query_atoms, num_pocket_atoms)

    # 2. Get feature matrices (precomputed or compute on-the-fly)
    if precomputed_matrices is not None:
        R_ij = precomputed_matrices['R_ij'].unsqueeze(0).expand(batch_size, -1, -1)
        is_hydrophobic = precomputed_matrices['is_hydrophobic'].unsqueeze(0).expand(batch_size, -1, -1)
        is_hbond = precomputed_matrices['is_hbond'].unsqueeze(0).expand(batch_size, -1, -1)
    else:
        # Fallback: compute on-the-fly (backward compatible)
        R_ij = query_features['vdw'].unsqueeze(1) + pocket_features['vdw'].unsqueeze(0)
        R_ij = R_ij.unsqueeze(0).expand(batch_size, -1, -1)

        q_hydro = query_features['hydro'].unsqueeze(1).expand(-1, num_pocket_atoms)
        p_hydro = pocket_features['hydro'].unsqueeze(0).expand(num_query_atoms, -1)
        is_hydrophobic = (q_hydro * p_hydro).unsqueeze(0).expand(batch_size, -1, -1)

        q_hbd = query_features['hbd'].unsqueeze(1).expand(-1, num_pocket_atoms)
        p_hba = pocket_features['hba'].unsqueeze(0).expand(num_query_atoms, -1)

        q_hba = query_features['hba'].unsqueeze(1).expand(-1, num_pocket_atoms)
        p_hbd = pocket_features['hbd'].unsqueeze(0).expand(num_query_atoms, -1)

        is_hbond = ((q_hbd * p_hba) + (q_hba * p_hbd) > 0).float()
        is_hbond = is_hbond.unsqueeze(0).expand(batch_size, -1, -1)

    # 3. Surface Distance
    delta_d = d_ij - R_ij

    # --- Vina Terms Calculation ---
    if weight_preset == 'vinardo':
        # 1. Gauss 1 (Vinardo uses width 0.8)
        gauss1 = torch.exp(-((delta_d / 0.8) ** 2))
        gauss2 = torch.zeros_like(delta_d)

        # 3. Repulsion - OPTIMIZED: use torch.where instead of masking
        repulsion_mask = delta_d < 0
        repulsion = torch.where(repulsion_mask, delta_d ** 2, torch.zeros_like(delta_d))

        # 4. Hydrophobic (Vinardo shifts cutoff from 1.5 to 2.5)
        hydro_1 = is_hydrophobic * (delta_d <= 0.0) * 1.0
        hydro_2_cond = is_hydrophobic * (delta_d > 0.0) * (delta_d < 2.5) * 1.0
        hydro_2 = hydro_2_cond * (1.0 - (delta_d / 2.5))
        hydrophobic_term = hydro_1 + hydro_2

        # 5. HBonding (Vinardo intercept is -0.6)
        hbond_1 = is_hbond * (delta_d <= -0.6) * 1.0
        hbond_2_cond = is_hbond * (delta_d < 0) * (delta_d > -0.6) * 1.0
        hbond_2 = hbond_2_cond * (-delta_d) / 0.6
        hbond_term = hbond_1 + hbond_2

    else:
        # Vina & Vina_LP Constants
        # 1 & 2. Gauss 1 and Gauss 2
        gauss1 = torch.exp(-((delta_d / 0.5) ** 2))

        # OPTIMIZED: use torch.where for gauss2
        gauss2_mask = delta_d != 0
        gauss2 = torch.where(gauss2_mask, torch.exp(-(((delta_d - 3.0) / 2.0) ** 2)), torch.zeros_like(delta_d))

        # 3. Repulsion - OPTIMIZED: use torch.where
        repulsion_mask = delta_d < 0
        repulsion = torch.where(repulsion_mask, delta_d ** 2, torch.zeros_like(delta_d))

        # 4. Hydrophobic
        hydro_1 = is_hydrophobic * (delta_d <= 0.5) * 1.0
        hydro_2_cond = is_hydrophobic * (delta_d > 0.5) * (delta_d < 1.5) * 1.0
        hydro_2 = 1.5 * hydro_2_cond - (hydro_2_cond * delta_d)
        hydrophobic_term = hydro_1 + hydro_2

        # 5. HBonding
        hbond_1 = is_hbond * (delta_d <= -0.7) * 1.0
        hbond_2_cond = is_hbond * (delta_d < 0) * (delta_d > -0.7) * 1.0
        hbond_2 = hbond_2_cond * (-delta_d) / 0.7
        hbond_term = hbond_1 + hbond_2

    # --- Combine into Total Energy ---
    w = VINA_WEIGHTS[weight_preset]

    energy_matrix = (w['gauss1'] * gauss1) + \
                    (w['gauss2'] * gauss2) + \
                    (w['repulsion'] * repulsion) + \
                    (w['hydrophobic'] * hydrophobic_term) + \
                    (w['hbond'] * hbond_term)

    # Sum over all query-pocket atom pairs for each conformer
    inter_energy = energy_matrix.sum(dim=(1, 2))
    total_energy = inter_energy

    # --- Intramolecular Energy Calculation (Optional) ---
    if intramolecular_mask is not None:
        # Distance between query atoms themselves
        intra_d_ij = torch.cdist(aligned_query_coords, aligned_query_coords) # [Batch, Num_Atoms, Num_Atoms]

        # VdW Sum
        intra_R_ij = query_features['vdw'].unsqueeze(1) + query_features['vdw'].unsqueeze(0)
        intra_R_ij = intra_R_ij.unsqueeze(0).expand(batch_size, -1, -1)
        intra_delta_d = intra_d_ij - intra_R_ij
        # Intramolecular Gauss and Repulsion - OPTIMIZED: use torch.where
        if weight_preset == 'vinardo':
            intra_gauss1 = torch.exp(-((intra_delta_d / 0.8) ** 2))
            intra_gauss2 = torch.zeros_like(intra_delta_d)
        else:
            intra_gauss1 = torch.exp(-((intra_delta_d / 0.5) ** 2))
            intra_gauss2_mask = intra_delta_d != 0
            intra_gauss2 = torch.where(intra_gauss2_mask, torch.exp(-(((intra_delta_d - 3.0) / 2.0) ** 2)), torch.zeros_like(intra_delta_d))

        intra_repulsion_mask = intra_delta_d < 0
        intra_repulsion = torch.where(intra_repulsion_mask, intra_delta_d ** 2, torch.zeros_like(intra_delta_d))

        # Only compute steric clash (Gauss + Repulsion) for intramolecular to save ops, which is the main stability factor
        intra_energy_matrix = (w['gauss1'] * intra_gauss1) + \
                              (w['gauss2'] * intra_gauss2) + \
                              (w['repulsion'] * intra_repulsion)

        # Apply mask and sum
        masked_intra = intra_energy_matrix * intramolecular_mask.unsqueeze(0).float()
        # Divide by 2 to prevent double counting
        intra_energy_sum = masked_intra.sum(dim=(1, 2)) / 2.0

        total_energy = total_energy + intra_energy_sum

    # Apply Vina Torsional Entropy Penalty if requested
    if num_rotatable_bonds is not None:
        total_energy = total_energy / (1.0 + w['rot'] * num_rotatable_bonds)

    return total_energy
