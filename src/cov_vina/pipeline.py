"""
High-level pipeline API for covalent docking.

This module provides the main entry-point:
  - run_covalent_pipeline : Covalent docking with warhead anchor
"""
import os
import time
import torch
from typing import Optional, Union, Literal
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from rdkit.Geometry import Point3D

from .molecular import generate_conformers_and_cluster, compute_vina_features
from .scoring import compute_intramolecular_mask, vina_scoring
from .optimization import optimize_torsions_vina
from .io import load_pocket_bundle, process_query_ligand, PocketBundle, final_selection, extract_pocket_around_residue
from .molecular.relax import relax_pose_with_fixed_core
from .molecular.anchor import (
    detect_warheads,
    find_reactive_residues,
    create_covalent_coordmap,
    check_warhead_residue_compatibility,
)
from .molecular.adduct import (
    create_adduct_template,
    get_covalent_exclusion_indices,
    get_protein_exclusion_residues,
    create_intermolecular_exclusion_mask,
)


def load_pocket_for_caching(
    protein_pdb: str,
    reactive_residue: Optional[str] = None,
    pocket_cutoff: float = 12.0,
    device: Optional[str] = None,
    verbose: bool = True,
):
    """
    Pre-load pocket for batch docking (caching optimization).

    This allows reusing pocket features across multiple ligands,
    avoiding redundant computation.

    Args:
        protein_pdb: Path to protein PDB file
        reactive_residue: e.g., "CYS145" or None for auto-detect
        pocket_cutoff: Pocket extraction radius (Å)
        device: 'cuda' or 'cpu'
        verbose: Print progress

    Returns:
        dict with keys:
            - pocket_bundle: PocketBundle with mol, coords, features
            - anchor: CovalentAnchor with residue info
            - residue_spec_str: Formatted residue string
            - device: torch device used
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    device = torch.device(device)

    # Load protein directly
    protein_mol = Chem.MolFromPDBFile(protein_pdb, sanitize=False, removeHs=True)
    if protein_mol is None:
        raise ValueError(f"Failed to load protein from {protein_pdb}")

    # Find anchor residue
    anchors = find_reactive_residues(protein_mol, residue_spec=reactive_residue)
    if not anchors:
        if reactive_residue:
            raise ValueError(f"Could not find reactive residue {reactive_residue}")
        else:
            raise ValueError("No reactive residues found. Specify with reactive_residue='CYS145'")

    anchor = anchors[0]
    residue_spec_str = f"{anchor.residue_name}{anchor.residue_num}"
    if anchor.chain_id.strip():
        residue_spec_str += f":{anchor.chain_id}"

    if verbose:
        print(f"Anchor: {anchor.residue_name}{anchor.residue_num}"
              f":{anchor.chain_id} atom {anchor.atom_name} "
              f"(bond length {anchor.bond_length:.2f} Å)")
        print(f"Extracting pocket within {pocket_cutoff}Å of {residue_spec_str}...")

    # Extract pocket
    pocket_mol = extract_pocket_around_residue(
        protein_mol,
        residue_spec_str,
        cutoff=pocket_cutoff
    )

    if verbose:
        print(f"  Pocket: {pocket_mol.GetNumAtoms()} atoms "
              f"(full protein: {protein_mol.GetNumAtoms()} atoms)")

    # Compute pocket features once
    pocket_coords = torch.tensor(
        pocket_mol.GetConformer().GetPositions(),
        dtype=torch.float32,
        device=device,
    )
    pocket_features = compute_vina_features(pocket_mol, device)
    pocket_bundle = PocketBundle(
        mol=pocket_mol,
        coords=pocket_coords,
        features=pocket_features,
    )

    # Re-find anchor in pocket
    anchors_pocket = find_reactive_residues(pocket_mol, residue_spec=residue_spec_str)
    if not anchors_pocket:
        raise RuntimeError(f"Failed to re-locate anchor {residue_spec_str} in extracted pocket")
    anchor = anchors_pocket[0]

    return {
        "pocket_bundle": pocket_bundle,
        "anchor": anchor,
        "residue_spec_str": residue_spec_str,
        "device": device,
    }


def run_covalent_pipeline(
    protein_pdb: str,
    query_ligand: str,
    reactive_residue: Optional[str] = None,
    output_dir: str = "output_predictions",
    # Pocket extraction
    pocket_cutoff: float = 12.0,
    # Cached pocket (for batch docking optimization)
    _cached_pocket: Optional[dict] = None,
    # Conformer generation
    num_confs: int = 1000,
    rmsd_threshold: float = 1.0,
    # Force field
    mmff_optimize: bool = True,
    # Optimization
    optimize: bool = False,
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
) -> dict:
    """
    Run covalent docking pipeline.

    Instead of MCS alignment to a reference ligand, this places the ligand's
    reactive warhead at the covalent bond distance from a protein nucleophile
    (e.g. CYS SG), then scores and optimises with Vina.

    Args:
        protein_pdb: Path to protein PDB file (full protein or pocket).
        query_ligand: SMILES string or path to SDF file.
        reactive_residue: Residue specifier, e.g. "CYS145", "CYS145:A".
            If None, auto-detect the first supported reactive residue.
        output_dir: Directory to save results.

        pocket_cutoff: Distance cutoff (Å) for pocket extraction around residue.
            Default 12Å covers Vina's effective scoring range (~10Å).
            Use 10Å for small molecules, 15-20Å for large ligands/peptides.

        num_confs: Number of conformers to generate.
        rmsd_threshold: RMSD threshold for clustering (Å).

        mmff_optimize: Apply MMFF94 relaxation after anchor placement.

        optimize: Enable gradient-based torsion optimization.
        optimizer: Optimizer type.
        opt_steps: Optimization steps.
        opt_lr: Learning rate.
        opt_batch_size: Batch size for optimization.

        weight_preset: Vina scoring weights.
        torsion_penalty: Apply torsional entropy penalty.

        save_all_poses: Save all poses or just top-k.
        top_k: Number of top poses to save.

        device: Device to use.
        verbose: Print progress messages.

    Returns:
        dict with output_file, num_poses, best_score, runtime, warhead_type,
        anchor_residue, etc.
    """
    t0 = time.time()

    if not verbose:
        RDLogger.DisableLog('rdApp.warning')

    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(device)

    if verbose:
        print(f"Using device: {device}")

    # ------------------------------------------------------------------ #
    # 1-2. Load protein and extract pocket (or use cached)
    # ------------------------------------------------------------------ #
    if _cached_pocket is not None:
        # Use pre-loaded pocket (batch docking optimization)
        pocket_bundle = _cached_pocket["pocket_bundle"]
        anchor = _cached_pocket["anchor"]
        residue_spec_str = _cached_pocket["residue_spec_str"]
        device = _cached_pocket["device"]
        pocket_mol = pocket_bundle.mol  # Extract mol from bundle

        if verbose:
            print(f"Using cached pocket: {pocket_mol.GetNumAtoms()} atoms")
            print(f"Anchor: {anchor.residue_name}{anchor.residue_num}"
                  f":{anchor.chain_id} atom {anchor.atom_name}")
    else:
        # Load and extract pocket from scratch
        if verbose:
            print(f"Loading protein from {protein_pdb}...")

        from .io import extract_pocket_around_residue
        protein_mol = Chem.MolFromPDBFile(protein_pdb, sanitize=False, removeHs=True)
        if protein_mol is None:
            raise ValueError(f"Failed to load protein from {protein_pdb}")

        # Find reactive residue
        anchors = find_reactive_residues(protein_mol, residue_spec=reactive_residue)
        if not anchors:
            raise ValueError(
                f"No reactive residue found in protein"
                + (f" matching '{reactive_residue}'" if reactive_residue else "")
                + ". Supported residue types: " + ", ".join(sorted(
                    f"{k} ({v.atom_name})"
                    for k, v in __import__('cov_vina.molecular.anchor',
                                           fromlist=['REACTIVE_RESIDUES']).REACTIVE_RESIDUES.items()
                ))
            )
        anchor = anchors[0]

        residue_spec_str = f"{anchor.residue_name}{anchor.residue_num}"
        if anchor.chain_id:
            residue_spec_str += f":{anchor.chain_id}"

        if verbose:
            print(f"Anchor: {anchor.residue_name}{anchor.residue_num}"
                  f":{anchor.chain_id} atom {anchor.atom_name} "
                  f"(bond length {anchor.bond_length:.2f} Å)")
            print(f"Extracting pocket within {pocket_cutoff}Å of {residue_spec_str}...")

        pocket_mol = extract_pocket_around_residue(
            protein_mol,
            residue_spec_str,
            cutoff=pocket_cutoff
        )

        if verbose:
            print(f"  Pocket: {pocket_mol.GetNumAtoms()} atoms "
                  f"(full protein: {protein_mol.GetNumAtoms()} atoms)")

        # Create pocket bundle
        pocket_coords = torch.tensor(
            pocket_mol.GetConformer().GetPositions(),
            dtype=torch.float32,
            device=device,
        )
        pocket_features = compute_vina_features(pocket_mol, device)
        pocket_bundle = PocketBundle(
            mol=pocket_mol,
            coords=pocket_coords,
            features=pocket_features,
        )

        # Re-find anchor in pocket
        anchors_pocket = find_reactive_residues(pocket_mol, residue_spec=residue_spec_str)
        if not anchors_pocket:
            raise RuntimeError(f"Failed to re-locate anchor {residue_spec_str} in extracted pocket")
        anchor = anchors_pocket[0]

    # ------------------------------------------------------------------ #
    # 3. Load & detect warhead on query ligand
    # ------------------------------------------------------------------ #
    query_mol, canonical_smiles = process_query_ligand(query_ligand)
    if verbose:
        print(f"Query ligand: {canonical_smiles}")

    # Detect on H-removed copy (atom indices will match after RemoveHs in conformer step)
    query_no_h = Chem.RemoveHs(Chem.MolFromSmiles(canonical_smiles))
    warhead_hits = detect_warheads(query_no_h)
    if not warhead_hits:
        raise ValueError("No reactive warhead detected on query ligand. "
                         "Cannot perform covalent docking.")

    warhead = warhead_hits[0]
    if verbose:
        print(f"Warhead: {warhead.warhead_type} "
              f"(reactive atom idx {warhead.reactive_atom_idx})")
        if len(warhead_hits) > 1:
            others = ", ".join(h.warhead_type for h in warhead_hits[1:])
            print(f"  Also detected: {others}")

    # Validate warhead-residue compatibility
    is_compatible, compat_msg = check_warhead_residue_compatibility(
        warhead.warhead_type, anchor.residue_name, strict=False
    )
    if not is_compatible:
        raise ValueError(f"Incompatible warhead-residue combination: {compat_msg}")
    if verbose and "Warning" in compat_msg:
        print(f"  {compat_msg}")

    # ------------------------------------------------------------------ #
    # 4. Create adduct template (BEFORE conformer generation)
    # ------------------------------------------------------------------ #
    if verbose:
        print("Creating adduct template (removing leaving group, adding CB-S)...")

    # Remove hydrogens first to match warhead detection
    query_mol = Chem.RemoveHs(query_mol)

    # Create adduct template (topology only, no conformers)
    adduct_mol, cb_atom_idx, nuc_atom_idx, new_reactive_idx = create_adduct_template(
        query_mol, warhead, anchor
    )

    if verbose:
        print(f"  Original atoms: {query_mol.GetNumAtoms()}")
        print(f"  Adduct atoms: {adduct_mol.GetNumAtoms()}")
        if cb_atom_idx is not None:
            print(f"  Added CB atom at index: {cb_atom_idx}")
        if nuc_atom_idx is not None:
            nuc_name = anchor.atom_name  # e.g., "SG", "OG", "NZ"
            print(f"  Added {nuc_name} atom at index: {nuc_atom_idx}")
        print(f"  Reactive atom index: {warhead.reactive_atom_idx} → {new_reactive_idx}")

    # Update warhead info with new reactive atom index
    from .molecular.anchor import WarheadHit
    warhead = WarheadHit(
        warhead_type=warhead.warhead_type,
        reactive_atom_idx=new_reactive_idx,
        matched_atoms=tuple()  # No longer meaningful after modification
    )

    # ------------------------------------------------------------------ #
    # 5. Generate conformers with CB-S coordmap constraint
    # ------------------------------------------------------------------ #
    if verbose:
        print(f"Generating conformers with CB-S fixed at anchor position...")

    # Create coordmap with CB and nucleophile fixed (NOT reactive atom - maximize diversity)
    coord_map = create_covalent_coordmap(cb_atom_idx, nuc_atom_idx, anchor)

    if verbose:
        print(f"  Fixed atoms in coordmap: {list(coord_map.keys())}")

    # Generate conformers on adduct molecule with CB-Nuc constraint
    adduct_mol, rep_cids = generate_conformers_and_cluster(
        adduct_mol, device, num_confs, rmsd_threshold, coordMap=coord_map
    )
    if len(rep_cids) == 0:
        raise RuntimeError("Failed to generate conformers")

    if verbose:
        nuc_name = anchor.atom_name
        print(f"Generated {len(rep_cids)} representative conformers (Butina clustering with CB-{nuc_name})")

    # Atoms to fix during MMFF: CB and nucleophile
    anchor_query_indices = {nuc_atom_idx}
    if cb_atom_idx is not None:
        anchor_query_indices.add(cb_atom_idx)

    # ------------------------------------------------------------------ #
    # 6. Optional MMFF relaxation with CB-Nuc fixed
    # ------------------------------------------------------------------ #
    if mmff_optimize and verbose:
        nuc_name = anchor.atom_name
        print(f"Relaxing conformers via MMFF94 with CB-{nuc_name} fixed...")

    batch_size = len(rep_cids)
    num_atoms = adduct_mol.GetNumAtoms()
    aligned_coords = torch.zeros((batch_size, num_atoms, 3))

    for j, cid in enumerate(rep_cids):
        conf = adduct_mol.GetConformer(cid)

        # Relax non-anchor atoms
        if mmff_optimize:
            applied, message = relax_pose_with_fixed_core(
                adduct_mol, cid, anchor_query_indices, max_iters=500,
            )
            if verbose:
                print(f"  Conformer {cid}: {message}")

        aligned_coords[j] = torch.tensor(conf.GetPositions(), dtype=torch.float32)

    aligned_coords = aligned_coords.to(device)

    # Replace query_mol with adduct
    query_mol = adduct_mol

    # Store metadata
    query_mol.SetProp("CovVina_Warhead_Type", warhead.warhead_type)
    query_mol.SetProp("CovVina_Reactive_Atom_Idx", str(warhead.reactive_atom_idx))
    query_mol.SetProp("CovVina_Anchor_Residue",
                       f"{anchor.residue_name}{anchor.residue_num}:{anchor.chain_id}")
    query_mol.SetProp("CovVina_Bond_Length", f"{anchor.bond_length:.2f}")

    # ------------------------------------------------------------------ #
    # 7. Vina Scoring with exclusion mask for covalent region
    # ------------------------------------------------------------------ #
    if verbose:
        print(f"Scoring with '{weight_preset}' weights "
              f"(excluding covalent bond region from intermolecular scoring)...")

    pocket_coords = pocket_bundle.coords
    pocket_features = pocket_bundle.features

    # Build query features (no masking - we use explicit exclusion mask instead)
    query_features = compute_vina_features(query_mol, device)

    # Create exclusion mask for covalent region
    ligand_exclude_indices = get_covalent_exclusion_indices(
        query_mol, warhead, n_hop_exclude=2
    )
    protein_exclude_residues = get_protein_exclusion_residues(anchor)

    intermolecular_exclusion_mask = create_intermolecular_exclusion_mask(
        query_mol, pocket_mol, ligand_exclude_indices,
        protein_exclude_residues, device
    )

    if verbose:
        num_excluded_pairs = intermolecular_exclusion_mask.sum().item()
        total_pairs = query_mol.GetNumAtoms() * pocket_mol.GetNumAtoms()
        print(f"  Excluding {num_excluded_pairs}/{total_pairs} atom pairs "
              f"({num_excluded_pairs/total_pairs*100:.1f}%)")

    num_rotatable_bonds = None
    if torsion_penalty:
        from rdkit.Chem import rdMolDescriptors
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(query_mol)

    intra_mask = compute_intramolecular_mask(query_mol, device)
    scores = vina_scoring(
        aligned_coords, pocket_coords, query_features, pocket_features,
        num_rotatable_bonds, weight_preset, intramolecular_mask=intra_mask,
        intermolecular_exclusion_mask=intermolecular_exclusion_mask,
    )
    initial_scores = scores.clone()

    # ------------------------------------------------------------------ #
    # 8. Optional gradient optimization
    # ------------------------------------------------------------------ #
    if optimize:
        if verbose:
            print(f"\n--- Gradient-Based Torsion Optimization "
                  f"({len(rep_cids)} poses) ---")

        # For covalent docking, freeze the anchor atom (protein Cβ)
        # This allows both Cβ-S and S-C bonds to rotate freely
        ref_indices = [cb_atom_idx] if cb_atom_idx is not None else [nuc_atom_idx]

        aligned_coords = optimize_torsions_vina(
            mol=query_mol,
            ref_indices=ref_indices,
            init_coords=aligned_coords,
            pocket_coords=pocket_coords,
            query_features=query_features,
            pocket_features=pocket_features,
            device=device,
            num_steps=opt_steps,
            lr=opt_lr,
            freeze_anchor=True,  # Always freeze anchor (protein Cβ atom)
            num_rotatable_bonds=num_rotatable_bonds,
            weight_preset=weight_preset,
            batch_size=opt_batch_size,
            optimizer=optimizer,
            intermolecular_exclusion_mask=intermolecular_exclusion_mask,
        )

        new_scores = vina_scoring(
            aligned_coords, pocket_coords, query_features, pocket_features,
            num_rotatable_bonds, weight_preset,
            intramolecular_mask=intra_mask,
            intermolecular_exclusion_mask=intermolecular_exclusion_mask,
        )

        score_diffs = new_scores - scores
        best_idx = torch.argmin(new_scores).item()
        if verbose:
            print(f"  Best pose: {best_idx} "
                  f"score {new_scores[best_idx]:.3f} kcal/mol "
                  f"(delta {score_diffs[best_idx]:.3f})")
            print(f"  Avg improvement: {score_diffs.mean():.3f} kcal/mol")
        scores = new_scores

        query_mol.SetProp("CovVina_Gradient_Optimized", "True")

    # ------------------------------------------------------------------ #
    # 9. Save results
    # ------------------------------------------------------------------ #
    if save_all_poses is None:
        save_all_poses = optimize
    if top_k is None:
        top_k = None if save_all_poses else 3

    os.makedirs(output_dir, exist_ok=True)

    if top_k is None:
        out_sdf = os.path.join(output_dir, "covalent_poses_all.sdf")
        final_selection(
            query_mol, rep_cids, aligned_coords, scores,
            initial_scores=initial_scores, top_k=None, output_path=out_sdf,
        )
        num_saved = len(rep_cids)
    else:
        out_sdf = os.path.join(output_dir, f"covalent_pose_top{top_k}.sdf")
        final_selection(
            query_mol, rep_cids, aligned_coords, scores,
            initial_scores=initial_scores, top_k=top_k, output_path=out_sdf,
        )
        num_saved = min(top_k, len(rep_cids))

    t1 = time.time()
    runtime = t1 - t0
    best_score = float(torch.min(scores).item())

    if verbose:
        print(f"\n-> Covalent docking completed in {runtime:.2f}s")
        print(f"-> Results saved to: {out_sdf}")

    return {
        "output_file": out_sdf,
        "num_poses": num_saved,
        "best_score": best_score,
        "runtime": runtime,
        "num_conformers": num_confs,
        "num_representatives": len(rep_cids),
        "warhead_type": warhead.warhead_type,
        "anchor_residue": f"{anchor.residue_name}{anchor.residue_num}",
        "anchor_atom": anchor.atom_name,
        "canonical_smiles": canonical_smiles,
        "device": str(device),
    }


def _mask_anchor_atom_features(query_features: dict, anchor_idx: int) -> None:
    """
    Zero-out Vina features for the anchor atom so it is excluded from
    intermolecular scoring.  The atom's VdW radius is set to 0 which
    effectively removes all distance-based terms for that pair.
    """
    query_features['vdw'][anchor_idx] = 0.0
    query_features['hydro'][anchor_idx] = 0.0
    query_features['hbd'][anchor_idx] = 0.0
    query_features['hba'][anchor_idx] = 0.0
