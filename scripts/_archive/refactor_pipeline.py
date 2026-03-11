#!/usr/bin/env python
"""Refactor pipeline.py to use CB-S fixed conformer generation."""

import re

# Read original file
with open("src/cov_vina/pipeline.py", "r") as f:
    content = f.read()

# Find the section to replace (from step 4 to before "# Store metadata")
# We need to replace from "# 4. Generate conformers" to just before "# Store metadata"

# Pattern: find section starting with "# 4. Generate conformers"
# and ending just before "# Store metadata"
pattern = re.compile(
    r'(    # -{70}\n'
    r'    # 4\. Generate conformers.*?\n'
    r'    # -{70}\n)'
    r'.*?'
    r'(    # Store metadata)',
    re.DOTALL
)

replacement = r'''\1    # ------------------------------------------------------------------ #
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

    \2'''

# Apply replacement
new_content = pattern.sub(replacement, content)

if new_content == content:
    print("ERROR: Pattern not found! No changes made.")
    exit(1)

# Write modified file
with open("src/cov_vina/pipeline.py", "w") as f:
    f.write(new_content)

print("Successfully refactored pipeline.py")
print("Old file backed up to pipeline.py.backup")
