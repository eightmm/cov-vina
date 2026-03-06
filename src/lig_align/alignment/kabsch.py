import torch
from typing import List, Tuple

def batched_kabsch_alignment(ref_coords: torch.Tensor, 
                             query_ensemble_coords: torch.Tensor, 
                             mapping: List[Tuple[int, int]],
                             device: torch.device) -> torch.Tensor:
    """
    Perform batched Kabsch algorithm aligning query conformers onto the reference ligand based on MCS coordinates.
    """
    ref_indices = [m[0] for m in mapping]
    query_indices = [m[1] for m in mapping]
    
    P = ref_coords[ref_indices].to(device).float()
    Q = query_ensemble_coords[:, query_indices, :].to(device).float()
    
    batch_size = Q.shape[0]
    
    # 1. Centering
    P_centroid = P.mean(dim=0, keepdim=True) # [1, 3]
    Q_centroid = Q.mean(dim=1, keepdim=True) # [Batch, 1, 3]
    
    P_centered = P - P_centroid # [N_mcs, 3]
    Q_centered = Q - Q_centroid # [Batch, N_mcs, 3]
    
    # 2. Covariance Matrix H = Q^T * P
    H = torch.bmm(Q_centered.transpose(1, 2), P_centered.expand(batch_size, -1, -1)) # [Batch, 3, 3]
    
    # 3. SVD
    U, S, Vh = torch.linalg.svd(H) 
    
    # 4. Compute Rotation Matrix R = V * U^T
    V = Vh.transpose(1, 2)
    Ut = U.transpose(1, 2)
    R = torch.bmm(V, Ut) # [Batch, 3, 3]
    
    # 5. Handle Reflection
    det = torch.linalg.det(R)
    reflection_mask = det < 0
    
    if reflection_mask.any():
        V_corr = V.clone()
        V_corr[reflection_mask, :, 2] *= -1
        R_corr = torch.bmm(V_corr, Ut)
        R[reflection_mask] = R_corr[reflection_mask]
        
    # 6. Translation vector
    t = P_centroid.unsqueeze(-1) - torch.bmm(R, Q_centroid.transpose(1, 2))
    
    # 7. Apply transformation to FULL query coordinates
    query_full = query_ensemble_coords.to(device).float()
    query_full_col = query_full.unsqueeze(-1)
    
    R_expand = R.unsqueeze(1)
    mapped_coords = torch.matmul(R_expand, query_full_col)
    
    t_expand = t.unsqueeze(1)
    mapped_coords = mapped_coords + t_expand
    
    return mapped_coords.squeeze(-1)
