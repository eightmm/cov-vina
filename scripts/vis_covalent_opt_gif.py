#!/usr/bin/env python3
"""
Generate 3D animated GIF visualization for covalent docking with live optimization.

Adapts vis_opt_gif.py for covalent docking - runs optimization and tracks trajectory in real-time.
"""

import argparse
import io
import os
import numpy as np
import torch
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger

from cov_vina.molecular import compute_vina_features
from cov_vina.alignment import LigandKinematics
from cov_vina.scoring import vina_scoring, compute_intramolecular_mask
from cov_vina.io import load_pocket_bundle, process_query_ligand
from cov_vina.io.visualization import draw_molecule_3d
from cov_vina.molecular.anchor import detect_warheads, find_reactive_residues, create_covalent_coordmap
from cov_vina.molecular.relax import relax_pose_with_fixed_core


def draw_covalent_adduct_3d(ax, mol, coords, anchor_idx, reactive_idx, label, alpha=1.0, cb_idx=None, s_idx=None):
    """
    Draw covalent adduct with color-coded atoms and bonds (ChimeraX-style).

    Color scheme (by element, similar to ChimeraX 'color byhetero'):
    - C: Gray (#909090)
    - N: Blue (#3050F8)
    - O: Red (#FF0D0D)
    - S: Yellow (#FFFF30)
    - F: Light green (#90E050)
    - Cl: Green (#1FF01F)
    - Br: Dark red (#A62929)
    - Cβ-S bond: Cyan (#00FFFF) - protein side chain
    - S-C bond: Magenta (#FF00FF) - covalent attachment
    - Rotatable bonds: Green (#40B040)
    - Heavy atoms only (no H), no element labels

    Highlights:
    - Cβ anchor: Orange with red edge
    - S atom: Gold with red edge
    """
    # ChimeraX-like element colors
    ELEMENT_COLORS = {
        'C': '#909090',   # Gray
        'N': '#3050F8',   # Blue
        'O': '#FF0D0D',   # Red
        'S': '#FFFF30',   # Yellow
        'F': '#90E050',   # Light green
        'Cl': '#1FF01F',  # Green
        'Br': '#A62929',  # Dark red
        'I': '#940094',   # Purple
        'P': '#FF8000',   # Orange
    }

    # Skip hydrogen atoms
    heavy_indices = [i for i in range(mol.GetNumAtoms())
                     if mol.GetAtomWithIdx(i).GetSymbol() != 'H']

    # Determine rotatable bonds
    from rdkit import Chem
    rot_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
    rotatable_bonds = set(mol.GetSubstructMatches(rot_smarts))

    # Draw bonds first (behind atoms)
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        # Skip if involves hydrogen
        if i not in heavy_indices or j not in heavy_indices:
            continue

        # Determine bond color
        if cb_idx is not None and s_idx is not None:
            # New scheme: Cβ-S-C chain
            if (i == cb_idx and j == s_idx) or (i == s_idx and j == cb_idx):
                # Cβ-S bond (protein side chain)
                bond_color = '#00FFFF'  # Cyan
                linewidth = 5
            elif (i == s_idx and j == reactive_idx) or (i == reactive_idx and j == s_idx):
                # S-C bond (covalent attachment to ligand)
                bond_color = '#FF00FF'  # Magenta
                linewidth = 5
            elif {i, j} in rotatable_bonds:
                # Rotatable bond
                bond_color = '#40B040'  # Green
                linewidth = 3.5
            else:
                # Normal bond
                bond_color = '#606060'  # Dark gray
                linewidth = 2.5
        else:
            # Old scheme (fallback)
            if (i == anchor_idx and j == reactive_idx) or (i == reactive_idx and j == anchor_idx):
                # C-S covalent bond - highlight
                bond_color = '#FF00FF'  # Magenta
                linewidth = 5
            elif {i, j} in rotatable_bonds:
                # Rotatable bond
                bond_color = '#40B040'  # Green
                linewidth = 3.5
            else:
                # Normal bond
                bond_color = '#606060'  # Dark gray
                linewidth = 2.5

        ax.plot([coords[i, 0], coords[j, 0]],
                [coords[i, 1], coords[j, 1]],
                [coords[i, 2], coords[j, 2]],
                color=bond_color, alpha=alpha, linewidth=linewidth)

    # Draw atoms (no labels, color by element)
    for i in heavy_indices:
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()

        # Get element color
        if cb_idx is not None and i == cb_idx:
            # Anchor Cβ atom - orange highlight
            color = '#FF8C00'  # Dark orange
            size = 150
            edgecolor = '#FF0000'  # Red edge
            edgewidth = 2.5
        elif s_idx is not None and i == s_idx:
            # S atom (part of covalent bond) - gold highlight
            color = '#FFD700'  # Gold
            size = 150
            edgecolor = '#FF6347'  # Tomato red edge
            edgewidth = 2.5
        elif i == anchor_idx and symbol == 'S':
            # Fallback: old anchor S highlighting
            color = '#FFD700'  # Gold
            size = 150
            edgecolor = '#FF6347'  # Tomato red edge
            edgewidth = 2.5
        else:
            # Normal atom - color by element
            color = ELEMENT_COLORS.get(symbol, '#FF1493')  # Default: hot pink
            # Size by element
            if symbol == 'C':
                size = 100
            elif symbol in ['N', 'O', 'S']:
                size = 120
            else:
                size = 130
            edgecolor = 'black'
            edgewidth = 1.0

        ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2],
                   c=color, s=size, alpha=alpha, edgecolors=edgecolor, linewidth=edgewidth, zorder=10)


def render_2d_warhead(mol: Chem.Mol, warhead_atoms: list[int], leaving_group_atoms: list[int] = None, width=400, height=400):
    """
    Render 2D structure with warhead highlighted.

    Args:
        mol: Original molecule (WITH leaving group)
        warhead_atoms: Warhead atoms to highlight in red
        leaving_group_atoms: Leaving group atoms to highlight in orange/yellow
    """
    mol_copy = Chem.RWMol(Chem.Mol(mol))
    mol_copy.RemoveAllConformers()
    rdDepictor.Compute2DCoords(mol_copy)

    d2d = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = d2d.drawOptions()
    opts.clearBackground = False

    # Red for warhead (reactive region)
    warhead_color = (0.85, 0.15, 0.15, 0.6)
    # Orange/Yellow for leaving group
    leaving_color = (1.0, 0.65, 0.0, 0.7)

    atom_colors = {}
    for i in warhead_atoms:
        atom_colors[i] = warhead_color

    if leaving_group_atoms:
        for i in leaving_group_atoms:
            atom_colors[i] = leaving_color

    highlight_bonds = []
    bond_colors = {}
    for bond in mol_copy.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        # Warhead bonds (red)
        if i in warhead_atoms and j in warhead_atoms:
            idx = bond.GetIdx()
            highlight_bonds.append(idx)
            bond_colors[idx] = warhead_color
        # Leaving group bonds (orange)
        elif leaving_group_atoms and i in leaving_group_atoms and j in leaving_group_atoms:
            idx = bond.GetIdx()
            highlight_bonds.append(idx)
            bond_colors[idx] = leaving_color
        # Bond between warhead and leaving group (gradient)
        elif leaving_group_atoms and ((i in warhead_atoms and j in leaving_group_atoms) or (j in warhead_atoms and i in leaving_group_atoms)):
            idx = bond.GetIdx()
            highlight_bonds.append(idx)
            bond_colors[idx] = (0.92, 0.40, 0.08, 0.65)  # Mix of red and orange

    rdMolDraw2D.PrepareAndDrawMolecule(
        d2d, mol_copy,
        highlightAtoms=list(atom_colors.keys()),
        highlightAtomColors=atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=bond_colors,
    )
    d2d.FinishDrawing()
    return Image.open(io.BytesIO(d2d.GetDrawingText()))


def main():
    parser = argparse.ArgumentParser(description="Visualize covalent docking optimization with 3D GIF")
    parser.add_argument("-q", "--query", required=True, help="Query SMILES with warhead")
    parser.add_argument("-o", "--output", default=None, help="Output GIF filename (saved to protein_dir/visualizations/)")
    parser.add_argument("-p", "--protein", required=True, help="Protein PDB path")
    parser.add_argument("-r", "--reactive_residue", default=None,
                       help="Reactive residue (e.g. 'CYS145' or 'CYS145:A'). Auto-detect if omitted.")
    parser.add_argument("-t", "--title", default="Covalent Docking", help="Plot title")
    parser.add_argument("--steps", type=int, default=200, help="Optimization steps")
    parser.add_argument("--free_anchor", action="store_true", help="Allow anchor atom to move")
    parser.add_argument("--weight_preset", type=str, choices=["vina", "vina_lp", "vinardo"],
                       default="vina", help="Vina weight preset")
    parser.set_defaults(torsion_penalty=True)
    torsion_group = parser.add_mutually_exclusive_group()
    torsion_group.add_argument("--torsion_penalty", dest="torsion_penalty", action="store_true",
                              help="Include torsional entropy penalty (default)")
    torsion_group.add_argument("--no_torsion_penalty", dest="torsion_penalty", action="store_false",
                              help="Disable torsional entropy penalty")

    args = parser.parse_args()

    # Auto-generate output path if not provided
    if args.output is None:
        # Extract warhead type from SMILES for filename
        protein_dir = os.path.dirname(os.path.abspath(args.protein))
        smiles_clean = args.query.replace("=", "").replace("(", "").replace(")", "")[:20]
        args.output = os.path.join(protein_dir, "visualizations", f"{smiles_clean}_covalent.gif")
    elif not os.path.isabs(args.output):
        # If relative path, put in protein_dir/visualizations/
        protein_dir = os.path.dirname(os.path.abspath(args.protein))
        args.output = os.path.join(protein_dir, "visualizations", args.output)

    # Create output directory
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    RDLogger.DisableLog('rdApp.warning')

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # 1. Load protein pocket
    print(f"Loading protein pocket from {args.protein}...")
    pocket_bundle = load_pocket_bundle(args.protein, device, lambda mol: compute_vina_features(mol, device))

    # 2. Process query ligand
    query_mol, canonical_smiles = process_query_ligand(args.query)
    print(f"Query ligand: {canonical_smiles}")

    # 3. Detect warhead
    query_mol_h = Chem.AddHs(query_mol)
    warhead_hits = detect_warheads(query_mol_h)

    if not warhead_hits:
        raise ValueError("No reactive warhead detected on query ligand. Cannot perform covalent docking.")

    warhead = warhead_hits[0]
    print(f"Warhead: {warhead.warhead_type} (reactive atom idx {warhead.reactive_atom_idx})")

    # Map back to heavy atoms
    heavy_map = {}
    heavy_i = 0
    for ai in range(query_mol_h.GetNumAtoms()):
        if query_mol_h.GetAtomWithIdx(ai).GetAtomicNum() != 1:
            heavy_map[ai] = heavy_i
            heavy_i += 1

    reactive_atom_idx_heavy = heavy_map.get(warhead.reactive_atom_idx)
    if reactive_atom_idx_heavy is None:
        raise ValueError("Reactive atom is hydrogen - cannot anchor")

    warhead_heavy_atoms = [heavy_map[ai] for ai in warhead.matched_atoms if ai in heavy_map]

    # 4. Find reactive residue
    reactive_residues = find_reactive_residues(pocket_bundle.mol, args.reactive_residue)
    if not reactive_residues:
        raise ValueError(f"No reactive residue found matching '{args.reactive_residue}'")

    anchor_res = reactive_residues[0]
    print(f"Anchor: {anchor_res.residue_name}{anchor_res.residue_num}:{anchor_res.chain_id} atom {anchor_res.atom_name}")

    # 5. Create coordMap for anchor constraint (use HEAVY atom index)
    # warhead.reactive_atom_idx is for H-含む molecule, need to map to heavy-only
    from cov_vina.molecular.anchor import WarheadHit
    from rdkit.Geometry import Point3D

    # Create coordMap with correct heavy-atom index
    target_pos = anchor_res.coord + anchor_res.bond_vector * anchor_res.bond_length
    coordMap = {
        reactive_atom_idx_heavy: Point3D(
            float(target_pos[0]),
            float(target_pos[1]),
            float(target_pos[2]),
        )
    }

    # 6. Generate conformer (coordMap often doesn't work, so we'll set manually)
    print(f"Generating conformer...")
    AllChem.EmbedMolecule(query_mol, randomSeed=42)

    if query_mol.GetNumConformers() == 0:
        raise ValueError("Failed to generate conformer")

    # Manually set anchor atom position
    conf = query_mol.GetConformer(0)
    conf.SetAtomPosition(
        reactive_atom_idx_heavy,
        Point3D(float(target_pos[0]), float(target_pos[1]), float(target_pos[2]))
    )

    # Verify
    anchor_atom_pos = conf.GetAtomPosition(reactive_atom_idx_heavy)
    print(f"  Anchor atom set to: ({anchor_atom_pos.x:.3f}, {anchor_atom_pos.y:.3f}, {anchor_atom_pos.z:.3f})")
    print(f"  Target (anchor+bond): ({target_pos[0]:.2f}, {target_pos[1]:.2f}, {target_pos[2]:.2f})")

    # 7. Relax with MMFF94 (freeze anchor)
    print("Relaxing with MMFF94 (anchor fixed)...")
    applied, msg = relax_pose_with_fixed_core(query_mol, 0, {reactive_atom_idx_heavy}, max_iters=500)
    print(f"  {msg}")

    init_coords = torch.tensor(query_mol.GetConformer().GetPositions(), dtype=torch.float32, device=device)

    # 7.5. Create adduct (remove leaving group)
    print("Creating covalent adduct (removing leaving group)...")
    from cov_vina.molecular.adduct import (
        create_covalent_adduct,
        get_covalent_exclusion_indices,
        get_protein_exclusion_residues,
        create_intermolecular_exclusion_mask,
        LEAVING_GROUPS
    )

    warhead_before_adduct = warhead
    query_mol_original = Chem.Mol(query_mol)  # Save original molecule for 2D rendering

    query_mol_adduct, cb_atom_idx, s_atom_idx = create_covalent_adduct(query_mol, warhead, anchor_res, add_anchor_atom=True)

    print(f"  Original atoms: {query_mol.GetNumAtoms()}")
    print(f"  Adduct atoms: {query_mol_adduct.GetNumAtoms()}")
    if cb_atom_idx is not None:
        print(f"  Added anchor atom (C\u03b2) at index: {cb_atom_idx}")
    if s_atom_idx is not None:
        print(f"  Added anchor atom (S) at index: {s_atom_idx}")

    # Replace query_mol with adduct (adduct already has conformer from create_covalent_adduct)
    query_mol = query_mol_adduct

    # Get coordinates from adduct molecule's conformer
    if query_mol_adduct.GetNumConformers() > 0:
        init_coords = torch.tensor(
            query_mol_adduct.GetConformer(0).GetPositions(),
            dtype=torch.float32,
            device=device
        )
    else:
        raise RuntimeError("Adduct molecule has no conformer!")

    # Determine leaving atoms for visualization
    leaving_pattern = LEAVING_GROUPS.get(warhead.warhead_type, [])
    matched_list = list(warhead.matched_atoms)
    leaving_atoms_heavy = [heavy_map[matched_list[i]] for i in leaving_pattern if matched_list[i] in heavy_map]

    # Calculate new reactive atom index after leaving group removal
    shift = sum(1 for idx in leaving_atoms_heavy if idx < warhead_before_adduct.reactive_atom_idx)
    new_reactive_idx = warhead_before_adduct.reactive_atom_idx - shift
    print(f"  Reactive atom index: {warhead_before_adduct.reactive_atom_idx} → {new_reactive_idx}")

    reactive_atom_idx_heavy = new_reactive_idx

    # Update warhead_heavy_atoms for visualization (approximate - first few atoms)
    warhead_heavy_atoms = list(range(min(3, query_mol.GetNumAtoms())))

    # 8. Setup kinematics
    # Use protein Cβ atom as anchor (freeze it to keep protein fixed)
    # This allows both Cβ-S and S-C bonds to rotate freely
    if cb_atom_idx is not None:
        anchor_indices = [cb_atom_idx]  # Freeze protein Cβ atom
        freeze_anchor = True
    elif s_atom_idx is not None:
        anchor_indices = [s_atom_idx]  # Fallback: freeze S if Cβ not available
        freeze_anchor = True
    else:
        # Fallback to old behavior
        freeze_anchor = not args.free_anchor
        anchor_indices = [] if not freeze_anchor else [reactive_atom_idx_heavy]

    model = LigandKinematics(
        query_mol,
        anchor_indices,
        init_coords,
        device,
        freeze_anchor=freeze_anchor
    )

    if cb_atom_idx is not None:
        print(f"  Freezing anchor atom (Cβ) at index {cb_atom_idx}, Cβ-S-C chain can rotate freely")
    elif s_atom_idx is not None:
        print(f"  Freezing anchor atom (S) at index {s_atom_idx}, all ligand atoms free to rotate")

    if model.num_torsions == 0:
        print("No rotatable bonds found. Skipping optimization.")
        return

    # 9. Prepare scoring
    num_rotatable_bonds = None
    if args.torsion_penalty:
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(query_mol)

    intra_mask = compute_intramolecular_mask(query_mol, device)
    pocket_coords = pocket_bundle.coords
    query_feat = compute_vina_features(query_mol, device)
    pocket_feat = pocket_bundle.features

    # Create exclusion mask for covalent region
    ligand_exclude_indices = get_covalent_exclusion_indices(
        query_mol, warhead_before_adduct, n_hop_exclude=2
    )
    protein_exclude_residues = get_protein_exclusion_residues(anchor_res)
    intermolecular_exclusion_mask = create_intermolecular_exclusion_mask(
        query_mol, pocket_bundle.mol, ligand_exclude_indices,
        protein_exclude_residues, device
    )
    print(f"  Exclusion mask: {intermolecular_exclusion_mask.sum().item()} pairs excluded")

    # 10. Optimize and track trajectory
    optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

    losses = []
    coords_history = []

    print(f"Optimizing for {args.steps} steps...")
    for step in range(args.steps):
        optimizer.zero_grad()
        coords = model()  # (N, 3)

        loss = vina_scoring(
            coords.unsqueeze(0), pocket_coords, query_feat, pocket_feat,
            num_rotatable_bonds, args.weight_preset,
            intramolecular_mask=intra_mask.unsqueeze(0),
            intermolecular_exclusion_mask=intermolecular_exclusion_mask
        )

        loss = loss.sum()
        loss.backward()
        optimizer.step()

        losses.append(loss.item())
        coords_history.append(coords.detach().cpu().numpy())

        if (step + 1) % 100 == 0:
            print(f"  Step {step+1}/{args.steps}: {loss.item():.3f} kcal/mol")

    init_numpy = init_coords.cpu().numpy()

    # 11. Generate 2x2 GIF animation
    print("Generating 2x2 GIF animation...")
    fig = plt.figure(figsize=(16, 12))

    ax_topleft = fig.add_subplot(2, 2, 1)
    ax_topright = fig.add_subplot(2, 2, 2)
    ax_botleft = fig.add_subplot(2, 2, 3)
    ax_botright = fig.add_subplot(2, 2, 4, projection='3d')

    # Render 2D structures: Original (left) and Adduct (right)
    # Original molecule with leaving group
    warhead_img_original = render_2d_warhead(query_mol_original, warhead_heavy_atoms, leaving_atoms_heavy, width=400, height=400)

    # Adduct molecule with S atom (no leaving group highlight)
    warhead_img_adduct = render_2d_warhead(query_mol, warhead_heavy_atoms, None, width=400, height=400)

    # Combine images side by side
    import numpy as np
    from PIL import Image as PILImage

    # Convert to numpy arrays
    img1_arr = np.array(warhead_img_original)
    img2_arr = np.array(warhead_img_adduct)

    # Concatenate horizontally
    combined_img = np.concatenate([img1_arr, img2_arr], axis=1)

    ax_topleft.imshow(combined_img)
    ax_topleft.axis('off')
    ax_topleft.set_title(f"Original: {canonical_smiles}  |  Adduct: +S (Anchor at {anchor_res.residue_name}{anchor_res.residue_num})",
                        fontsize=11, fontweight='bold', color='#2c3e50')

    # Metadata panel
    ax_topright.axis('off')
    metadata_text = f"""
COVALENT DOCKING

Warhead:       {warhead.warhead_type}
Reactive Atom: {reactive_atom_idx_heavy}
Anchor:        {anchor_res.residue_name}{anchor_res.residue_num}:{anchor_res.chain_id}
Anchor Atom:   {anchor_res.atom_name}

MOLECULE
Heavy Atoms:      {query_mol.GetNumAtoms()}
Rotatable Bonds:  {num_rotatable_bonds if num_rotatable_bonds else 0}
Flexible Torsions: {model.num_torsions}

OPTIMIZATION
Steps:         {args.steps}
Freeze Anchor: {'Yes' if freeze_anchor else 'No'}
Initial Score: {losses[0]:.2f} kcal/mol
Final Score:   {losses[-1]:.2f} kcal/mol
Delta:         {losses[-1] - losses[0]:.2f} kcal/mol
"""

    ax_topright.text(0.5, 0.5, metadata_text, fontsize=11, ha='center', va='center',
                    family='monospace',
                    bbox=dict(boxstyle='round,pad=1', facecolor='lightgray', alpha=0.3))

    # 3D bounds
    pad = 2.0
    all_coords = init_numpy
    min_x, max_x = all_coords[:, 0].min() - pad, all_coords[:, 0].max() + pad
    min_y, max_y = all_coords[:, 1].min() - pad, all_coords[:, 1].max() + pad
    min_z, max_z = all_coords[:, 2].min() - pad, all_coords[:, 2].max() + pad

    def update(frame):
        ax_botleft.clear()
        ax_botright.clear()

        # Energy landscape
        ax_botleft.plot(losses[:frame+1], color='#e74c3c', linewidth=3)
        ax_botleft.set_xlim(0, args.steps)
        padding = 5
        ax_botleft.set_ylim(min(losses) - padding, max(losses) + padding)
        ax_botleft.set_title(f"Energy Minimization (Step {frame+1}/{args.steps})",
                            fontsize=14, fontweight='bold')
        ax_botleft.set_xlabel("Optimization Step", fontsize=12)
        ax_botleft.set_ylabel("Vina Energy Score (kcal/mol)", fontsize=12)
        ax_botleft.grid(True, linestyle='--', alpha=0.7)
        ax_botleft.scatter([frame], [losses[frame]], color='black', s=100, zorder=5)
        ax_botleft.text(frame + 2, losses[frame], f"{losses[frame]:.2f}",
                       fontsize=12, fontweight='bold')

        # 3D structure overlay
        ax_botright.set_title(f"3D Optimization Trajectory\n({model.num_torsions} Flexible Torsions, Heavy Atoms Only)",
                             fontsize=13, fontweight='bold')
        ax_botright.set_xlim(min_x, max_x)
        ax_botright.set_ylim(min_y, max_y)
        ax_botright.set_zlim(min_z, max_z)

        # Initial (gray, faint)
        draw_covalent_adduct_3d(ax_botright, query_mol, init_numpy,
                                anchor_idx=cb_atom_idx if cb_atom_idx is not None else (s_atom_idx if s_atom_idx is not None else reactive_atom_idx_heavy),
                                reactive_idx=reactive_atom_idx_heavy,
                                label="Initial", alpha=0.2,
                                cb_idx=cb_atom_idx, s_idx=s_atom_idx)

        # Current (colored, solid with legend)
        draw_covalent_adduct_3d(ax_botright, query_mol, coords_history[frame],
                                anchor_idx=cb_atom_idx if cb_atom_idx is not None else (s_atom_idx if s_atom_idx is not None else reactive_atom_idx_heavy),
                                reactive_idx=reactive_atom_idx_heavy,
                                label="Current", alpha=1.0,
                                cb_idx=cb_atom_idx, s_idx=s_atom_idx)

        # Custom legend with color codes
        from matplotlib.patches import Patch
        if cb_atom_idx is not None and s_atom_idx is not None:
            legend_elements = [
                Patch(facecolor='#FF8C00', label='Cβ Atom (Anchor)'),
                Patch(facecolor='#FFD700', label='S Atom'),
                Patch(facecolor='#00FFFF', label='Cβ-S Bond'),
                Patch(facecolor='#FF00FF', label='S-C Bond'),
                Patch(facecolor='#40B040', label='Rotatable Bonds'),
                Patch(facecolor='#606060', label='Other Bonds')
            ]
        else:
            legend_elements = [
                Patch(facecolor='#FFD700', label='S Atom (Anchor)'),
                Patch(facecolor='#00CED1', label='Reactive C'),
                Patch(facecolor='#FF00FF', label='C-S Bond'),
                Patch(facecolor='#32CD32', label='Rotatable Bonds'),
                Patch(facecolor='#808080', label='Other Bonds')
            ]
        ax_botright.legend(handles=legend_elements, loc='upper right', fontsize=7, framealpha=0.9)
        ax_botright.set_xticks([])
        ax_botright.set_yticks([])
        ax_botright.set_zticks([])

        return ax_botleft, ax_botright

    ani = FuncAnimation(fig, update, frames=len(coords_history), interval=100, blit=False)
    ani.save(args.output, writer=PillowWriter(fps=10))
    print(f"✓ GIF saved to {args.output}")


if __name__ == "__main__":
    main()
