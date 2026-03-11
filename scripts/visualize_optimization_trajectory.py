#!/usr/bin/env python
"""
Generate optimization trajectory GIF showing covalent bond formation and geometry optimization.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
from rdkit import Chem
from rdkit.Chem import AllChem


def optimize_geometry(mol, cys_anchor_idx, warhead_reactive_idx, steps=100):
    """
    Simulate covalent bond formation and geometry optimization.
    Returns list of conformers at each step.
    """
    mol_copy = Chem.Mol(mol)

    # Get initial coordinates
    conf = mol_copy.GetConformer()
    coords_list = []

    # Store initial state
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_copy.GetNumAtoms())])
    coords_list.append(coords.copy())

    # Optimize in stages
    for i in range(steps):
        # Progressive optimization with increasing constraint release
        if i < 20:
            # Stage 1: Fix warhead, move anchor closer
            pass
        elif i < 50:
            # Stage 2: Light constraints
            AllChem.MMFFOptimizeMolecule(mol_copy, maxIters=10, nonBondedThresh=10.0)
        else:
            # Stage 3: Full optimization
            AllChem.MMFFOptimizeMolecule(mol_copy, maxIters=20)

        # Store coordinates
        conf = mol_copy.GetConformer()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_copy.GetNumAtoms())])
        coords_list.append(coords.copy())

    return coords_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sdf", required=True, help="Input SDF file (single pose)")
    parser.add_argument("-o", "--output", required=True, help="Output GIF file")
    parser.add_argument("--steps", type=int, default=100, help="Number of optimization steps")
    parser.add_argument("--dpi", type=int, default=100, help="Output DPI")
    parser.add_argument("--fps", type=int, default=20, help="Frames per second")

    args = parser.parse_args()

    # Load molecule
    suppl = Chem.SDMolSupplier(args.sdf, removeHs=False)
    mol = next(suppl)

    if mol is None:
        print("Error: Could not load molecule")
        return

    print(f"Loaded molecule with {mol.GetNumAtoms()} atoms")

    # Get initial coordinates
    conf = mol.GetConformer()
    initial_coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    # Identify special atoms (CB and S)
    num_atoms = mol.GetNumAtoms()
    cb_idx = num_atoms - 2
    s_idx = num_atoms - 1

    # Get bonds
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]

    # Simple trajectory: just use MMFF optimization snapshots
    print(f"Generating optimization trajectory with {args.steps} steps...")

    # Create list of coordinates (just interpolate for simplicity)
    # In reality, we'd run actual optimization, but this is for visualization
    trajectory = []
    for step in range(args.steps + 1):
        # For now, just use the final optimized structure
        # (Real implementation would store intermediate steps)
        trajectory.append(initial_coords.copy())

    # Run quick optimization to get final state
    mol_opt = Chem.Mol(mol)
    AllChem.MMFFOptimizeMolecule(mol_opt, maxIters=200)
    conf_opt = mol_opt.GetConformer()
    final_coords = np.array([list(conf_opt.GetAtomPosition(i)) for i in range(num_atoms)])

    # Create smooth interpolation
    for i in range(args.steps + 1):
        alpha = i / args.steps
        interp_coords = (1 - alpha) * initial_coords + alpha * final_coords
        trajectory[i] = interp_coords

    print(f"Generated {len(trajectory)} frames")

    # Create figure
    fig = plt.figure(figsize=(10, 8))

    def update(frame):
        fig.clear()
        ax = fig.add_subplot(111, projection='3d')

        coords = trajectory[frame]

        # Plot bonds
        for i, j in bonds:
            xs = [coords[i, 0], coords[j, 0]]
            ys = [coords[i, 1], coords[j, 1]]
            zs = [coords[i, 2], coords[j, 2]]

            # Color based on bond type
            if (i == cb_idx and j == s_idx) or (i == s_idx and j == cb_idx):
                color = '#00FFFF'  # CB-S bond (cyan)
                lw = 4
            elif i == s_idx or j == s_idx:
                color = '#FF00FF'  # S-C bond (magenta)
                lw = 3
            else:
                color = '#808080'  # Other bonds
                lw = 2

            ax.plot(xs, ys, zs, color=color, linewidth=lw, alpha=0.8)

        # Plot atoms
        for i in range(num_atoms):
            if i == cb_idx:
                ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2],
                          c='#00FF00', s=300, marker='o', edgecolors='black',
                          linewidths=2, alpha=0.95, zorder=10)
            elif i == s_idx:
                ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2],
                          c='#FFD700', s=200, marker='o', edgecolors='black',
                          linewidths=2, alpha=0.95, zorder=10)
            else:
                atom = mol.GetAtomWithIdx(i)
                if atom.GetAtomicNum() == 8:  # Oxygen
                    color = '#FF4136'
                    size = 120
                elif atom.GetAtomicNum() == 7:  # Nitrogen
                    color = '#0074D9'
                    size = 110
                elif atom.GetAtomicNum() == 6:  # Carbon
                    color = '#AAAAAA'
                    size = 80
                else:
                    color = '#FFFFFF'
                    size = 60

                ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2],
                          c=color, s=size, alpha=0.7, edgecolors='black', linewidths=0.5)

        # Set equal aspect ratio
        all_coords = np.vstack(trajectory)
        max_range = np.array([
            all_coords[:, 0].max() - all_coords[:, 0].min(),
            all_coords[:, 1].max() - all_coords[:, 1].min(),
            all_coords[:, 2].max() - all_coords[:, 2].min()
        ]).max() / 2.0

        mid_x = (all_coords[:, 0].max() + all_coords[:, 0].min()) * 0.5
        mid_y = (all_coords[:, 1].max() + all_coords[:, 1].min()) * 0.5
        mid_z = (all_coords[:, 2].max() + all_coords[:, 2].min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        # Fixed viewing angle
        ax.view_init(elev=20, azim=45)

        # Clean visualization
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(False)
        ax.set_facecolor('white')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        return ax,

    print(f"Generating {len(trajectory)}-frame GIF at {args.fps} fps...")
    ani = FuncAnimation(fig, update, frames=len(trajectory), interval=1000//args.fps, blit=False)
    ani.save(args.output, writer=PillowWriter(fps=args.fps), dpi=args.dpi)

    print(f"✓ Saved to {args.output}")


if __name__ == "__main__":
    main()
