"""
Temporary script to extract and refactor the core of pipeline.py steps 4-6.
This will be manually copied back to pipeline.py after verification.
"""

# NEW STEPS 4-6 for pipeline.py (to replace lines ~202-379)

STEPS_4_TO_6 = '''
    # ------------------------------------------------------------------ #
    # 4. Create adduct template (BEFORE conformer generation)
    # ------------------------------------------------------------------ #
    if verbose:
        print("Creating adduct template (removing leaving group, adding CB-S)...")

    # Remove hydrogens first to match warhead detection
    query_mol = Chem.RemoveHs(query_mol)

    # Create adduct template (topology only, no conformers)
    adduct_mol, cb_atom_idx, s_atom_idx, new_reactive_idx = create_adduct_template(
        query_mol, warhead, anchor
    )

    if verbose:
        print(f"  Original atoms: {query_mol.GetNumAtoms()}")
        print(f"  Adduct atoms: {adduct_mol.GetNumAtoms()}")
        if cb_atom_idx is not None:
            print(f"  Added CB atom at index: {cb_atom_idx}")
        if s_atom_idx is not None:
            print(f"  Added S atom at index: {s_atom_idx}")
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

    # Create coordmap with CB and S fixed (NOT reactive atom - maximize diversity)
    coord_map = create_covalent_coordmap(cb_atom_idx, s_atom_idx, anchor)

    if verbose:
        print(f"  Fixed atoms in coordmap: {list(coord_map.keys())}")

    # Generate conformers on adduct molecule with CB-S constraint
    adduct_mol, rep_cids = aligner.step1_generate_conformers(
        adduct_mol, num_confs=num_confs, rmsd_threshold=rmsd_threshold,
        coordMap=coord_map,
    )
    if len(rep_cids) == 0:
        raise RuntimeError("Failed to generate conformers")

    if verbose:
        print(f"Generated {len(rep_cids)} representative conformers (Butina clustering with CB-S)")

    # ------------------------------------------------------------------ #
    # 6. Optional MMFF relaxation with CB-S fixed
    # ------------------------------------------------------------------ #
    if mmff_optimize and verbose:
        print("Relaxing conformers via MMFF94 with CB-S fixed...")

    batch_size = len(rep_cids)
    num_atoms = adduct_mol.GetNumAtoms()
    aligned_coords = torch.zeros((batch_size, num_atoms, 3))

    # Atoms to fix during MMFF: CB and S
    anchor_query_indices = {s_atom_idx}
    if cb_atom_idx is not None:
        anchor_query_indices.add(cb_atom_idx)

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

    aligned_coords = aligned_coords.to(aligner.device)

    # Replace query_mol with adduct
    query_mol = adduct_mol
'''

print("New steps 4-6 defined. Manual integration required.")
print("Lines to replace in pipeline.py: approximately 202-379")
print(STEPS_4_TO_6)
