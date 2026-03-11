#!/usr/bin/env python
"""Delete old sections 5-5.5 from pipeline.py and replace with simple MMFF section."""

# Read file
with open("src/cov_vina/pipeline.py", "r") as f:
    lines = f.readlines()

# Find START: "    # 5. Exact coordinate placement"
# Find END: "    # Store metadata" (after anchor_query_indices = {...})
start_idx = None
end_idx = None

for i, line in enumerate(lines):
    if "# 5. Exact coordinate placement" in line:
        start_idx = i - 1  # Include the comment separator line above
    if start_idx is not None and "# Store metadata" in line:
        end_idx = i
        break

if start_idx is None or end_idx is None:
    print(f"ERROR: Could not find section boundaries")
    print(f"start_idx={start_idx}, end_idx={end_idx}")
    exit(1)

print(f"Deleting lines {start_idx+1} to {end_idx} (0-indexed: {start_idx}-{end_idx-1})")

# New section to insert
new_section = """    # ------------------------------------------------------------------ #
    # 6. Optional MMFF relaxation with CB-S fixed
    # ------------------------------------------------------------------ #
    if mmff_optimize and verbose:
        print("Relaxing conformers via MMFF94 with CB-S fixed...")

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

    aligned_coords = aligned_coords.to(aligner.device)

    # Replace query_mol with adduct
    query_mol = adduct_mol

    # Store metadata
"""

# Reconstruct file
new_lines = lines[:start_idx] + [new_section] + lines[end_idx:]

# Write back
with open("src/cov_vina/pipeline.py", "w") as f:
    f.writelines(new_lines)

print(f"Successfully replaced lines {start_idx+1}-{end_idx}")
print(f"Inserted {len(new_section.splitlines())} new lines")
