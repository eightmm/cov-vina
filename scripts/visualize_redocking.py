#!/usr/bin/env python
"""
Visualize crystal ligand redocking results: success vs failure cases.
"""
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch


def kabsch_align(P, Q):
    """Align Q to P using Kabsch algorithm."""
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    H = Q_c.T @ P_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    Q_aligned = Q_c @ R + P.mean(axis=0)
    return Q_aligned


def kabsch_rmsd(P, Q):
    """Calculate RMSD after Kabsch alignment."""
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    H = Q_c.T @ P_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return np.sqrt(np.mean(np.sum((P_c - Q_c @ R)**2, axis=1)))


def plot_overlay(ax, crystal_coords, pose_coords, title, rmsd, score):
    """Plot 3D overlay of crystal and docked pose."""
    # Plot crystal (green)
    ax.plot(crystal_coords[:, 0], crystal_coords[:, 1], crystal_coords[:, 2],
            'o-', color='#2ECC40', markersize=6, linewidth=1.5, label='Crystal', alpha=0.8)

    # Plot docked (orange/red depending on RMSD)
    color = '#FF4136' if rmsd > 3.0 else '#FF851B'
    ax.plot(pose_coords[:, 0], pose_coords[:, 1], pose_coords[:, 2],
            'o-', color=color, markersize=6, linewidth=1.5, label='Docked', alpha=0.8)

    # Highlight key atoms
    ax.scatter(crystal_coords[0, 0], crystal_coords[0, 1], crystal_coords[0, 2],
               s=150, c='#2ECC40', marker='*', edgecolors='black', linewidths=1.5, zorder=10)
    ax.scatter(pose_coords[0, 0], pose_coords[0, 1], pose_coords[0, 2],
               s=150, c=color, marker='*', edgecolors='black', linewidths=1.5, zorder=10)

    ax.set_title(f'{title}\nRMSD: {rmsd:.2f} Å, Score: {score:.2f} kcal/mol',
                 fontsize=11, fontweight='bold')
    ax.set_xlabel('X (Å)', fontsize=9)
    ax.set_ylabel('Y (Å)', fontsize=9)
    ax.set_zlabel('Z (Å)', fontsize=9)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    # Set equal aspect ratio
    max_range = np.array([
        crystal_coords[:, 0].max() - crystal_coords[:, 0].min(),
        crystal_coords[:, 1].max() - crystal_coords[:, 1].min(),
        crystal_coords[:, 2].max() - crystal_coords[:, 2].min()
    ]).max() / 2.0

    mid_x = (crystal_coords[:, 0].max() + crystal_coords[:, 0].min()) * 0.5
    mid_y = (crystal_coords[:, 1].max() + crystal_coords[:, 1].min()) * 0.5
    mid_z = (crystal_coords[:, 2].max() + crystal_coords[:, 2].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)


def main():
    parser = argparse.ArgumentParser(description="Visualize redocking success and failure cases")
    parser.add_argument("-r", "--reference", required=True, help="Reference PDB (crystal)")
    parser.add_argument("-q", "--query", required=True, help="Docked poses SDF")
    parser.add_argument("-o", "--output", required=True, help="Output PNG file")
    parser.add_argument("--dpi", type=int, default=150, help="Output DPI")

    args = parser.parse_args()

    # Load crystal
    print(f"Loading crystal from {args.reference}...")
    crystal = Chem.MolFromPDBFile(args.reference, removeHs=False)

    # Load poses
    print(f"Loading docked poses from {args.query}...")
    suppl = Chem.SDMolSupplier(args.query, removeHs=False)
    poses = [mol for mol in suppl if mol is not None]
    print(f"  Loaded {len(poses)} poses")

    # Find MCS
    print("Finding maximum common substructure...")
    mcs = rdFMCS.FindMCS([crystal, poses[0]], timeout=5, completeRingsOnly=False)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    print(f"  MCS: {mcs.numAtoms} atoms")

    # Get crystal MCS coordinates
    crystal_match = crystal.GetSubstructMatch(mcs_mol)
    crystal_conf = crystal.GetConformer()
    crystal_coords = np.array([list(crystal_conf.GetAtomPosition(i)) for i in crystal_match])

    # Calculate RMSD for all poses
    print("Calculating RMSD for all poses...")
    rmsds = []
    scores = []
    for pose in poses:
        pose_match = pose.GetSubstructMatch(mcs_mol)
        pose_conf = pose.GetConformer()
        pose_coords = np.array([list(pose_conf.GetAtomPosition(i)) for i in pose_match])
        rmsd = kabsch_rmsd(crystal_coords, pose_coords)
        rmsds.append(rmsd)
        scores.append(float(pose.GetProp('Vina_Score_Final')))

    # Find best and worst cases
    best_rmsd_idx = np.argmin(rmsds)
    worst_rmsd_idx = np.argmax(rmsds)
    best_score_idx = 0  # Top ranked by score

    print(f"\nBest RMSD: Pose {best_rmsd_idx+1}, RMSD={rmsds[best_rmsd_idx]:.3f} Å, Score={scores[best_rmsd_idx]:.3f}")
    print(f"Best Score: Pose {best_score_idx+1}, RMSD={rmsds[best_score_idx]:.3f} Å, Score={scores[best_score_idx]:.3f}")
    print(f"Worst RMSD: Pose {worst_rmsd_idx+1}, RMSD={rmsds[worst_rmsd_idx]:.3f} Å, Score={scores[worst_rmsd_idx]:.3f}")

    # Create figure
    fig = plt.figure(figsize=(16, 5))

    # Plot 1: Best RMSD case (success)
    ax1 = fig.add_subplot(131, projection='3d')
    pose_match = poses[best_rmsd_idx].GetSubstructMatch(mcs_mol)
    pose_conf = poses[best_rmsd_idx].GetConformer()
    pose_coords = np.array([list(pose_conf.GetAtomPosition(i)) for i in pose_match])
    pose_coords_aligned = kabsch_align(crystal_coords, pose_coords)
    plot_overlay(ax1, crystal_coords, pose_coords_aligned,
                 f'Success: Best RMSD (Rank {best_rmsd_idx+1})',
                 rmsds[best_rmsd_idx], scores[best_rmsd_idx])

    # Plot 2: Best score case
    ax2 = fig.add_subplot(132, projection='3d')
    pose_match = poses[best_score_idx].GetSubstructMatch(mcs_mol)
    pose_conf = poses[best_score_idx].GetConformer()
    pose_coords = np.array([list(pose_conf.GetAtomPosition(i)) for i in pose_match])
    pose_coords_aligned = kabsch_align(crystal_coords, pose_coords)
    quality = "Success" if rmsds[best_score_idx] < 2.5 else "Moderate"
    plot_overlay(ax2, crystal_coords, pose_coords_aligned,
                 f'{quality}: Best Score (Rank 1)',
                 rmsds[best_score_idx], scores[best_score_idx])

    # Plot 3: Worst RMSD case (failure)
    ax3 = fig.add_subplot(133, projection='3d')
    pose_match = poses[worst_rmsd_idx].GetSubstructMatch(mcs_mol)
    pose_conf = poses[worst_rmsd_idx].GetConformer()
    pose_coords = np.array([list(pose_conf.GetAtomPosition(i)) for i in pose_match])
    pose_coords_aligned = kabsch_align(crystal_coords, pose_coords)
    plot_overlay(ax3, crystal_coords, pose_coords_aligned,
                 f'Failure: Worst RMSD (Rank {worst_rmsd_idx+1})',
                 rmsds[worst_rmsd_idx], scores[worst_rmsd_idx])

    # Add overall title
    fig.suptitle(f'7AEH Crystal Ligand Redocking: MCS Alignment ({mcs.numAtoms} atoms)\n'
                 f'Mean RMSD: {np.mean(rmsds):.2f} ± {np.std(rmsds):.2f} Å  |  '
                 f'Best RMSD: {min(rmsds):.2f} Å  |  '
                 f'{len([r for r in rmsds if r < 2.5])}/{len(rmsds)} poses < 2.5 Å',
                 fontsize=13, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    print(f"\n✓ Saved to {args.output}")


if __name__ == "__main__":
    main()
