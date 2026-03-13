"""
Batch docking API - automatically handles pocket caching.

This module provides a simple interface for batch docking that
automatically optimizes performance by caching pocket features.
"""
from typing import List, Dict, Union, Optional, Literal
from pathlib import Path

from .pipeline import run_covalent_pipeline, load_pocket_for_caching


def run_batch_docking(
    protein_pdb: str,
    ligands: Union[List[str], str],  # Single SMILES, List of SMILES, or path to .smi file
    reactive_residue: Optional[str] = None,
    output_dir: str = "results",
    # Pocket extraction
    pocket_cutoff: float = 12.0,
    # Conformer generation
    num_confs: int = 200,
    rmsd_threshold: float = 1.0,
    # Optimization
    optimize: bool = True,
    optimizer: Literal["adam", "adamw", "lbfgs"] = "adam",
    opt_steps: int = 100,
    opt_lr: float = 0.05,
    opt_batch_size: int = 128,
    # Scoring
    weight_preset: Literal["vina", "vina_lp", "vinardo"] = "vina",
    torsion_penalty: bool = True,
    # Output
    save_all_poses: Optional[bool] = None,
    top_k: Optional[int] = None,
    # Device
    device: Optional[str] = None,
    verbose: bool = True,
) -> List[Dict]:
    """
    Run covalent docking on multiple ligands with automatic pocket caching.

    This function automatically handles pocket caching for optimal performance.
    Just pass a list of SMILES and it takes care of the rest!

    Args:
        protein_pdb: Path to protein PDB file (full protein or pocket)
        ligands: Either:
            - Single SMILES string: "C=CC(=O)NC"
            - List of SMILES strings: ["C=CC(=O)NC", "O=CCc1ccccc1", ...]
            - Path to .smi file: "ligands.smi"
        reactive_residue: e.g., "CYS145" or None for auto-detect
        output_dir: Directory for results (single ligand) or base directory (multiple ligands)

        [... all other parameters same as run_covalent_pipeline ...]

        verbose: Print progress for each ligand

    Returns:
        List of result dicts, one per ligand:
            {
                'ligand': SMILES string,
                'name': ligand name (if from .smi) or index,
                'success': True/False,
                'output_file': path to final_poses.sdf (if successful),
                'best_score': float (if successful),
                'warhead_type': str (if successful),
                'num_poses': int (if successful),
                'runtime': float in seconds,
                'error': str (if failed),
            }

    Example:
        >>> from cov_vina import run_batch_docking
        >>>
        >>> # Single ligand (same interface!)
        >>> results = run_batch_docking(
        ...     protein_pdb="pocket.pdb",
        ...     ligands="C=CC(=O)NC",  # Single SMILES works!
        ...     reactive_residue="CYS145",
        ... )
        >>>
        >>> # Multiple ligands with list
        >>> results = run_batch_docking(
        ...     protein_pdb="pocket.pdb",
        ...     ligands=["C=CC(=O)NC", "O=CCc1ccccc1"],
        ...     reactive_residue="CYS145",
        ... )
        >>>
        >>> # Or with .smi file
        >>> results = run_batch_docking(
        ...     protein_pdb="pocket.pdb",
        ...     ligands="ligands.smi",
        ...     reactive_residue="CYS145",
        ... )
        >>>
        >>> # Check results
        >>> for r in results:
        ...     if r['success']:
        ...         print(f"{r['name']}: {r['best_score']:.3f} kcal/mol")
    """
    import os
    import time

    # Parse ligands input
    if isinstance(ligands, str):
        # Check if it's a file path or single SMILES
        if ligands.endswith('.smi') or ligands.endswith('.smiles') or '\n' not in ligands and ' ' in ligands:
            # Looks like a .smi file path
            ligand_list = _parse_smi_file(ligands)
        else:
            # Single SMILES string
            ligand_list = [(ligands, "ligand_1")]
    else:
        # List of SMILES
        ligand_list = [(smiles, f"ligand_{i+1}") for i, smiles in enumerate(ligands)]

    if verbose:
        print(f"Processing {len(ligand_list)} ligands...")
        print(f"Output directory: {output_dir}/")

    # Pre-load pocket once (automatic caching!)
    if verbose:
        print(f"\nLoading protein pocket: {protein_pdb}")
        print(f"Reactive residue: {reactive_residue or 'auto-detect'}")

    cache_start = time.time()
    cached_pocket = load_pocket_for_caching(
        protein_pdb=protein_pdb,
        reactive_residue=reactive_residue,
        pocket_cutoff=pocket_cutoff,
        device=device,
        verbose=verbose,
    )
    cache_time = time.time() - cache_start

    if verbose:
        print(f"✓ Pocket loaded in {cache_time:.2f}s")
        print(f"  Pocket atoms: {cached_pocket['pocket_bundle'].mol.GetNumAtoms()}")
        print(f"  Device: {cached_pocket['device']}")

    # GPU warmup to avoid overhead on first ligand
    if len(ligand_list) > 1:  # Only worthwhile for batch
        from .utils import warmup_gpu
        warmup_time = warmup_gpu(cached_pocket['device'], verbose=verbose)
        if verbose:
            print()
    else:
        warmup_time = 0.0

    # Process each ligand
    results = []
    total_start = time.time()

    for i, (smiles, name) in enumerate(ligand_list, 1):
        if verbose:
            print(f"[{i}/{len(ligand_list)}] {name}: {smiles}")

        # Create output directory
        mol_out_dir = os.path.join(output_dir, name)
        os.makedirs(mol_out_dir, exist_ok=True)

        mol_start = time.time()
        result = {
            'ligand': smiles,
            'name': name,
            'success': False,
            'runtime': 0.0,
        }

        try:
            # Run docking with cached pocket
            docking_result = run_covalent_pipeline(
                protein_pdb=protein_pdb,
                query_ligand=smiles,
                reactive_residue=reactive_residue,
                output_dir=mol_out_dir,
                pocket_cutoff=pocket_cutoff,
                num_confs=num_confs,
                rmsd_threshold=rmsd_threshold,
                optimize=optimize,
                optimizer=optimizer,
                opt_steps=opt_steps,
                opt_lr=opt_lr,
                opt_batch_size=opt_batch_size,
                weight_preset=weight_preset,
                torsion_penalty=torsion_penalty,
                save_all_poses=save_all_poses,
                top_k=top_k,
                device=device,
                verbose=False,  # Reduce clutter
                _cached_pocket=cached_pocket,  # Use cached pocket!
            )

            # Success
            result['success'] = True
            result['output_file'] = docking_result['output_file']
            result['best_score'] = docking_result['best_score']
            result['warhead_type'] = docking_result['warhead_type']
            result['num_poses'] = docking_result['num_poses']
            result['runtime'] = time.time() - mol_start

            if verbose:
                print(f"  ✓ Score: {result['best_score']:.3f} kcal/mol "
                      f"({result['num_poses']} poses, {result['runtime']:.2f}s)")

        except Exception as e:
            # Failed
            result['error'] = str(e)
            result['runtime'] = time.time() - mol_start

            if verbose:
                print(f"  ✗ Failed: {str(e)[:60]} ({result['runtime']:.2f}s)")

            # Save error log
            error_log = os.path.join(mol_out_dir, "error.log")
            with open(error_log, 'w') as f:
                f.write(f"SMILES: {smiles}\n")
                f.write(f"Error: {str(e)}\n")

        results.append(result)

    total_time = time.time() - total_start
    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful

    if verbose:
        print()
        print("=" * 60)
        print("Batch docking complete!")
        print(f"  Successful: {successful}/{len(results)}")
        print(f"  Failed: {failed}/{len(results)}")
        print(f"  Total time: {total_time:.2f}s")
        print(f"  Avg time/ligand: {total_time/len(results):.2f}s")
        print(f"  (Pocket load: {cache_time:.2f}s, Docking: {total_time-cache_time:.2f}s)")
        print("=" * 60)

    return results


def _parse_smi_file(smi_path: str) -> List[tuple]:
    """
    Parse .smi file and return list of (SMILES, name) tuples.

    Format:
        SMILES_string    molecule_name
        # Comments starting with # are ignored
    """
    molecules = []
    with open(smi_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                smiles = parts[0]
                name = parts[1]
                molecules.append((smiles, name))
            elif len(parts) == 1:
                # No name provided, use index
                smiles = parts[0]
                name = f"molecule_{len(molecules)+1}"
                molecules.append((smiles, name))

    return molecules
