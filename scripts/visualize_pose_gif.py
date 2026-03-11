#!/usr/bin/env python
"""
Generate rotating 3D GIF visualization for a docked pose.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sdf", required=True, help="Input SDF file")
    parser.add_argument("-o", "--output", required=True, help="Output GIF file")
    parser.add_argument("-r", "--reference", help="Reference crystal PDB file")
    parser.add_argument("--dpi", type=int, default=100, help="Output DPI")
    parser.add_argument("--frames", type=int, default=72, help="Number of frames")
    parser.add_argument("--fps", type=int, default=15, help="Frames per second")

    args = parser.parse_args()

    # Load molecule
    suppl = Chem.SDMolSupplier(args.sdf, removeHs=False)
    mol = next(suppl)

    if mol is None:
        print("Error: Could not load molecule")
        return

    print(f"Loaded molecule with {mol.GetNumAtoms()} atoms")

    # Get properties
    try:
        score = float(mol.GetProp('Vina_Score_Final'))
        warhead = mol.GetProp('CovVina_Warhead_Type')
        rank_str = mol.GetProp('Rank')
        print(f"  Rank: {rank_str}")
        print(f"  Score: {score:.3f} kcal/mol")
        print(f"  Warhead: {warhead}")
    except:
        score = 0.0
        warhead = "unknown"
        rank_str = "?"

    # Extract coordinates
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

    # Load reference crystal if provided
    crystal_coords = None
    crystal_bonds = None
    if args.reference:
        crystal = Chem.MolFromPDBFile(args.reference, removeHs=False)
        if crystal:
            print(f"Loaded reference crystal with {crystal.GetNumAtoms()} atoms")
            crystal_conf = crystal.GetConformer()
            crystal_coords = np.array([list(crystal_conf.GetAtomPosition(i)) for i in range(crystal.GetNumAtoms())])
            crystal_bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in crystal.GetBonds()]
        else:
            print("Warning: Could not load reference crystal")

    # Identify special atoms (last 2 are CB and S)
    num_atoms = mol.GetNumAtoms()
    cb_idx = num_atoms - 2
    s_idx = num_atoms - 1

    # Build bond list
    bonds = []
    for bond in mol.GetBonds():
        bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    # Create figure
    fig = plt.figure(figsize=(10, 8))

    def update(frame):
        fig.clear()
        ax = fig.add_subplot(111, projection='3d')

        # Rotation angle
        angle = frame * (360 / args.frames)
        ax.view_init(elev=20, azim=angle)

        # Plot crystal structure first (if available)
        if crystal_coords is not None and crystal_bonds is not None:
            # Plot crystal bonds (thin green lines)
            for i, j in crystal_bonds:
                xs = [crystal_coords[i, 0], crystal_coords[j, 0]]
                ys = [crystal_coords[i, 1], crystal_coords[j, 1]]
                zs = [crystal_coords[i, 2], crystal_coords[j, 2]]
                ax.plot(xs, ys, zs, color='#00FF00', linewidth=1, alpha=0.3)

            # Plot crystal atoms (small green spheres)
            ax.scatter(crystal_coords[:, 0], crystal_coords[:, 1], crystal_coords[:, 2],
                      c='#00FF00', s=40, alpha=0.3, label='Crystal', zorder=1)

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
                          linewidths=2, alpha=0.95, label='CB (anchor)', zorder=10)
            elif i == s_idx:
                ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2],
                          c='#FFD700', s=200, marker='o', edgecolors='black',
                          linewidths=2, alpha=0.95, label='S (covalent)', zorder=10)
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
        max_range = np.array([
            coords[:, 0].max() - coords[:, 0].min(),
            coords[:, 1].max() - coords[:, 1].min(),
            coords[:, 2].max() - coords[:, 2].min()
        ]).max() / 2.0

        mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
        mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
        mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        # Remove axis labels for cleaner look
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

    print(f"Generating {args.frames}-frame GIF at {args.fps} fps...")
    ani = FuncAnimation(fig, update, frames=args.frames, interval=1000//args.fps, blit=False)
    ani.save(args.output, writer=PillowWriter(fps=args.fps), dpi=args.dpi)

    print(f"✓ Saved to {args.output}")


if __name__ == "__main__":
    main()
