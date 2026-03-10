# API Reference

## Covalent Docking API

```python
from cov_vina import run_covalent_pipeline

results = run_covalent_pipeline(
    protein_pdb="pocket.pdb",
    query_ligand="C=CC(=O)Nc1ccccc1",  # must contain a reactive warhead
    reactive_residue="CYS145",          # or None for auto-detect
    output_dir="output",
    num_confs=1000,
    optimize=True,
    optimizer="lbfgs",
)
```

Return keys:

```python
{
    "output_file": "output/covalent_pose_top3.sdf",
    "num_poses": 3,
    "best_score": -5.2,
    "runtime": 12.3,
    "num_conformers": 1000,
    "num_representatives": 18,
    "warhead_type": "acrylamide",
    "anchor_residue": "CYS145",
    "anchor_atom": "SG",
    "canonical_smiles": "C=CC(=O)Nc1ccccc1",
    "device": "cuda",
}
```

### Covalent-Specific Parameters

```text
protein_pdb        Protein PDB path
query_ligand       Query SMILES or SDF path (must have a warhead)
reactive_residue   Residue spec, e.g. "CYS145", "CYS145:A", or None (auto-detect)
freeze_anchor      Keep anchor atom fixed during optimization (default: True)
```

### Warhead Detection API

```python
from cov_vina.molecular.anchor import detect_warheads, find_reactive_residues
from rdkit import Chem

# Detect warheads on a ligand
mol = Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1")
hits = detect_warheads(mol)
# hits[0].warhead_type == "acrylamide"
# hits[0].reactive_atom_idx == 0

# Find reactive residues in a protein pocket
pocket = Chem.MolFromPDBFile("pocket.pdb", sanitize=False, removeHs=True)
anchors = find_reactive_residues(pocket)
# anchors[0].residue_name == "CYS"
# anchors[0].residue_num == 145
# anchors[0].coord == array([x, y, z])
```

## MCS-Guided Docking API

```python
from cov_vina import run_pipeline

results = run_pipeline(
    protein_pdb="protein.pdb",
    ref_ligand="ref.sdf",
    query_ligand="SMILES",
    output_dir="output",
    num_confs=1000,
    rmsd_threshold=1.0,
    optimize=True,
    optimizer="lbfgs",
    verbose=True,
)
```

Return keys:

```python
{
    "output_file": "output/predicted_pose_top3.sdf",
    "num_poses": 3,
    "best_score": -6.038,
    "runtime": 23.5,
    "num_conformers": 1000,
    "num_representatives": 22,
    "mcs_size": 10,
    "mcs_positions": 1,
    "canonical_smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "device": "cuda",
}
```

## Common Parameters

```text
protein_pdb        Protein PDB path
query_ligand       Query SMILES or SDF path
output_dir         Output directory
num_confs          Number of conformers to generate
rmsd_threshold     RMSD clustering threshold
optimize           Enable torsion optimization
optimizer          adam | adamw | lbfgs
opt_steps          Number of optimization steps
opt_lr             Optimization learning rate
opt_batch_size     Number of poses processed per optimization batch (default: 128)
weight_preset      vina | vina_lp | vinardo
torsion_penalty    Apply torsional entropy penalty (default: True)
verbose            Print progress
```

MCS-only parameters:

```text
ref_ligand         Reference SDF path
mcs_mode           auto | single | multi | cross
freeze_mcs         Keep MCS atoms fixed during optimization
```

Covalent-only parameters:

```text
reactive_residue   Residue spec or None
freeze_anchor      Keep anchor atom fixed during optimization
```

## Low-Level API

For stepwise control, use `LigandAligner`.

```python
from cov_vina import LigandAligner

aligner = LigandAligner(device="cuda")
mapping = aligner.step2_find_mcs(ref_mol, query_mol)
query_mol, rep_cids = aligner.step1_generate_conformers(
    query_mol,
    num_confs=1000,
    rmsd_threshold=1.0,
)
aligned = aligner.step3_batched_kabsch_alignment(ref_coords, query_coords, mapping)
scores = aligner.step4_vina_scoring(aligned, pocket_coords, query_feat, pocket_feat)
optimized = aligner.step6_refine_pose(
    query_mol,
    mcs_indices,
    aligned,
    pocket_coords,
    query_feat,
    pocket_feat,
    num_steps=100,
    batch_size=8,
)
```

## Script Inventory

- `scripts/run_covalent_pipeline.py`: covalent docking with warhead anchor
- `scripts/run_pipeline.py`: MCS-based end-to-end pose generation
- `scripts/optimize_pose.py`: optimize a single input pose
- `scripts/vis_comparison_grid.py`: generate comparison panels across optimization settings
- `scripts/vis_opt_gif.py`: make optimization animations
- `scripts/vis_ref_opt_gif.py`: compare reference-guided optimization trajectories
