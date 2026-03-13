"""
Utility functions for CovVina optimization.
"""
import torch


def warmup_gpu(device: torch.device, verbose: bool = False):
    """
    Pre-compile CUDA kernels to avoid warmup overhead on first ligand.

    This triggers JIT compilation of common operations used in the pipeline:
    - Distance matrix computation (RMSD clustering)
    - Optimizer initialization (Adam/LBFGS)
    - Tensor operations

    Args:
        device: torch device (cuda or cpu)
        verbose: print progress

    Returns:
        float: warmup time in seconds
    """
    import time

    if device.type == 'cpu':
        return 0.0  # No warmup needed for CPU

    start = time.time()

    if verbose:
        print("Warming up GPU kernels...")

    # Trigger RMSD kernel compilation
    dummy = torch.randn(100, 3, device=device)
    _ = torch.cdist(dummy, dummy)

    # Trigger optimizer compilation
    dummy_param = torch.randn(10, 3, device=device, requires_grad=True)
    optimizer = torch.optim.Adam([dummy_param], lr=0.1)
    loss = dummy_param.sum()
    loss.backward()
    optimizer.step()

    # Trigger common tensor ops
    _ = torch.einsum('ij,jk->ik', dummy[:10, :], dummy[:3, :].T)
    _ = torch.norm(dummy, dim=1)

    # Synchronize to ensure all kernels are compiled
    if device.type == 'cuda':
        torch.cuda.synchronize()

    elapsed = time.time() - start

    if verbose:
        print(f"  ✓ GPU warmed up in {elapsed:.2f}s")

    return elapsed

