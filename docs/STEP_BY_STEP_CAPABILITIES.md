# LigAlign: Step-by-Step Capabilities

Complete reference of all functions and parameters available at each pipeline step.

---

## Step 0: Initialization

### `LigandAligner(device=None)`

Initialize the aligner with device selection.

**Parameters:**
- `device`: `None`, `"cuda"`, or `"cpu"` (default: auto-detect)

**Example:**
```python
from lig_align import LigandAligner

aligner = LigandAligner()              # Auto-detect
aligner = LigandAligner(device='cuda') # Force GPU
aligner = LigandAligner(device='cpu')  # Force CPU
```

---

## Step 1: Conformer Generation

### `step1_generate_conformers(mol, num_confs, rmsd_threshold, coordMap)`

Generate conformers and cluster by RMSD.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mol` | `Chem.Mol` | Required | Query molecule (must have H atoms) |
| `num_confs` | `int` | `1000` | Number of conformers to generate |
| `rmsd_threshold` | `float` | `1.0` | RMSD clustering threshold (Å) |
| `coordMap` | `dict` | `None` | Optional coordinate constraints `{atom_idx: Point3D}` |

**Returns:**
- `mol`: Molecule with conformers added
- `representative_cids`: List of cluster centroid conformer IDs

**What it does:**
1. Generate N conformers using RDKit ETKDG
2. Optional MMFF minimize (if no coordMap)
3. Remove hydrogens
4. GPU-accelerated pairwise RMSD calculation (centered, no alignment)
5. Butina clustering with RMSD threshold
6. Return cluster representatives

**Example:**
```python
from rdkit import Chem
from rdkit.Chem import AllChem

query_mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
query_mol = Chem.AddHs(query_mol)

# Basic usage
mol, cids = aligner.step1_generate_conformers(
    query_mol,
    num_confs=1000,
    rmsd_threshold=1.0
)
print(f"Generated {len(cids)} cluster representatives")

# With coordinate constraints (for MCS alignment)
from rdkit.Geometry import Point3D
coordMap = {
    0: Point3D(1.0, 2.0, 3.0),
    1: Point3D(4.0, 5.0, 6.0)
}
mol, cids = aligner.step1_generate_conformers(
    query_mol,
    num_confs=1000,
    rmsd_threshold=1.0,
    coordMap=coordMap
)
```

**Tips:**
- Higher `num_confs` → better sampling (slower)
- Lower `rmsd_threshold` → more representatives (more diverse, slower)
- Typical: `num_confs=1000`, `rmsd_threshold=1.0` (robust default)
- For fast screening: `num_confs=100`, `rmsd_threshold=2.0`
- For high accuracy: `num_confs=3000`, `rmsd_threshold=0.5`

---

## Step 2: MCS Search

### `step2_find_mcs(ref_mol, query_mol, return_all_positions, cross_match, min_fragment_size, max_fragments)`

Find Maximum Common Substructure with three modes.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ref_mol` | `Chem.Mol` | Required | Reference molecule |
| `query_mol` | `Chem.Mol` | Required | Query molecule |
| `return_all_positions` | `bool` | `False` | Return all MCS positions (Mode 2) |
| `cross_match` | `bool` | `False` | Cross-matching mode (Mode 3) |
| `min_fragment_size` | `int` | `None` | Min atoms per fragment (Mode 3) |
| `max_fragments` | `int` | `3` | Max fragments to find (Mode 3) |

**Returns:**
- **Mode 1** (default): `[(ref_idx, query_idx), ...]` - Single mapping
- **Mode 2/3**: `[mapping1, mapping2, ...]` - List of mappings

**Three Modes:**

#### Mode 1: Single Position (1:1) - DEFAULT
```python
mapping = aligner.step2_find_mcs(ref_mol, query_mol)
# Returns: [(0, 0), (1, 1), (2, 2), ...]
# Single best MCS alignment
```

**Use when:**
- Reference is asymmetric
- You want fastest performance
- Default behavior (backward compatible)

#### Mode 2: Multi-Position (1:N) - Symmetric Reference
```python
mappings = aligner.step2_find_mcs(
    ref_mol,
    query_mol,
    return_all_positions=True
)
# Returns: [mapping1, mapping2, ...]
# All possible positions where query MCS matches reference
```

**Use when:**
- Reference is symmetric (e.g., dibenzyl with two phenyl rings)
- Query can match multiple positions
- You want to explore all alignment options

**Example:**
```python
# Reference: dibenzyl (Ph-CH2-CH2-Ph)
# Query: benzene
mappings = aligner.step2_find_mcs(ref_mol, query_mol, return_all_positions=True)
# Returns 2 mappings: one for each phenyl ring
```

#### Mode 3: Cross-Matching (N:M) - Both Symmetric
```python
mappings = aligner.step2_find_mcs(
    ref_mol,
    query_mol,
    cross_match=True,
    min_fragment_size=5,
    max_fragments=3
)
# Returns: [combo1, combo2, ...]
# All cross-combinations of multiple MCS fragments
```

**Use when:**
- Both reference AND query have multiple similar fragments
- Disconnected MCS regions
- Complex multi-fragment alignment

**Parameters for Mode 3:**
- `min_fragment_size`: Minimum atoms for a fragment (default: 5)
- `max_fragments`: Maximum fragments to find (default: 3)

**Example:**
```python
# Both molecules have multiple phenyl rings
mappings = aligner.step2_find_mcs(
    ref_mol,
    query_mol,
    cross_match=True,
    min_fragment_size=5,
    max_fragments=3
)
# Returns all valid cross-matching combinations
```

**Performance:**
- Mode 1: ~0.12 ms (<0.01% of pipeline)
- Mode 2: ~0.13 ms (<0.01% of pipeline)
- Mode 3: ~0.50 ms (<0.01% of pipeline)

---

## Step 3: Kabsch Alignment

### `step3_batched_kabsch_alignment(ref_coords, query_ensemble_coords, mapping)`

Align query conformer ensemble to reference using MCS mapping.

**Parameters:**
| Parameter | Type | Shape | Description |
|-----------|------|-------|-------------|
| `ref_coords` | `Tensor` | `[N_atoms, 3]` | Reference coordinates |
| `query_ensemble_coords` | `Tensor` | `[N_confs, N_atoms, 3]` | Query conformers |
| `mapping` | `List[Tuple]` | - | MCS mapping from step2 |

**Returns:**
- `aligned_coords`: `[N_confs, N_atoms, 3]` - Aligned coordinates

**What it does:**
1. Extract MCS atoms from both molecules
2. Compute optimal rotation/translation (SVD)
3. Apply transformation to all query atoms
4. GPU-accelerated batch processing

**Example:**
```python
import torch

# Get reference coordinates
ref_coords = torch.tensor(ref_mol.GetConformer().GetPositions())

# Get query ensemble coordinates
query_coords = torch.zeros((len(cids), query_mol.GetNumAtoms(), 3))
for i, cid in enumerate(cids):
    query_coords[i] = torch.tensor(query_mol.GetConformer(cid).GetPositions())

# Align
aligned = aligner.step3_batched_kabsch_alignment(
    ref_coords,
    query_coords,
    mapping
)
# aligned.shape: [N_representatives, N_atoms, 3]
```

---

## Step 4: Vina Scoring

### `step4_vina_scoring(aligned_query_coords, pocket_coords, query_features, pocket_features, num_rotatable_bonds, weight_preset, intramolecular_mask)`

Compute differentiable Vina score.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `aligned_query_coords` | `Tensor` | Required | Aligned query coordinates `[N_poses, N_atoms, 3]` |
| `pocket_coords` | `Tensor` | Required | Pocket coordinates `[N_pocket_atoms, 3]` |
| `query_features` | `dict` | Required | Query Vina features |
| `pocket_features` | `dict` | Required | Pocket Vina features |
| `num_rotatable_bonds` | `int` | `None` | For torsional penalty |
| `weight_preset` | `str` | `"vina"` | Scoring weights: `"vina"`, `"vina_lp"`, `"vinardo"` |
| `intramolecular_mask` | `Tensor` | `None` | Mask for intramolecular interactions |

**Returns:**
- `scores`: `[N_poses]` - Vina scores in kcal/mol

**Compute Vina Features:**
```python
query_features = aligner.compute_vina_features(query_mol)
pocket_features = aligner.compute_vina_features(pocket_mol)
```

**Feature Dictionary:**
```python
{
    'coords': Tensor,      # [N, 3]
    'vdw_radii': Tensor,   # [N]
    'atom_types': Tensor,  # [N] (0-4: C, N, O, S, other)
    'h_bond_donor': Tensor,    # [N] binary
    'h_bond_acceptor': Tensor, # [N] binary
    'hydrophobic': Tensor      # [N] binary
}
```

**Intramolecular Mask:**
```python
from lig_align.scoring import compute_intramolecular_mask

# Exclude bonded atoms (1-2, 1-3, 1-4 neighbors)
intra_mask = compute_intramolecular_mask(query_mol, device)
```

**Scoring Weights:**

| Preset | Description | Use Case |
|--------|-------------|----------|
| `"vina"` | Standard AutoDock Vina | General purpose (default) |
| `"vina_lp"` | Vina with local preference | Prioritize local interactions |
| `"vinardo"` | Vinardo scoring function | Alternative scoring |

**Example:**
```python
from lig_align.scoring import compute_intramolecular_mask

# Compute features
query_features = aligner.compute_vina_features(query_mol)
pocket_features = aligner.compute_vina_features(pocket_mol)

# Intramolecular mask
intra_mask = compute_intramolecular_mask(query_mol, aligner.device)

# Score all poses
scores = aligner.step4_vina_scoring(
    aligned_coords,
    pocket_coords,
    query_features,
    pocket_features,
    num_rotatable_bonds=None,
    weight_preset="vina",
    intramolecular_mask=intra_mask
)

# Find best pose
best_idx = torch.argmin(scores)
print(f"Best score: {scores[best_idx]:.2f} kcal/mol")
```

**With Torsional Penalty:**
```python
from rdkit.Chem import rdMolDescriptors

num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(query_mol)

scores = aligner.step4_vina_scoring(
    aligned_coords,
    pocket_coords,
    query_features,
    pocket_features,
    num_rotatable_bonds=num_rotatable,  # Add penalty
    weight_preset="vina"
)
```

---

## Step 5: Final Selection & Output

### `step5_final_selection(mol, representative_cids, aligned_coords, scores, top_k, output_path)`

Select and save top-k poses.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mol` | `Chem.Mol` | Required | Query molecule |
| `representative_cids` | `List[int]` | Required | Conformer IDs |
| `aligned_coords` | `Tensor` | Required | Coordinates `[N_poses, N_atoms, 3]` |
| `scores` | `Tensor` | Required | Vina scores `[N_poses]` |
| `top_k` | `int` | `None` | Number of poses to save (None = all) |
| `output_path` | `str` | `"output.sdf"` | Output SDF file path |

**Returns:**
- `selected_indices`: Indices of selected poses

**Example:**
```python
# Save top 3 poses
aligner.step5_final_selection(
    query_mol,
    representative_cids,
    aligned_coords,
    scores,
    top_k=3,
    output_path="top3_poses.sdf"
)

# Save all poses
aligner.step5_final_selection(
    query_mol,
    representative_cids,
    aligned_coords,
    scores,
    top_k=None,
    output_path="all_poses.sdf"
)
```

**Output SDF Properties:**
- `Vina_Score`: Energy in kcal/mol
- `MCS_Num_Atoms`: MCS size
- `MCS_Ref_Coverage`: % of reference covered
- `MCS_Query_Coverage`: % of query covered
- Additional pipeline parameters

---

## Step 6: Gradient Optimization

### `step6_refine_pose(mol, ref_indices, init_coords, pocket_coords, query_features, pocket_features, num_steps, lr, freeze_mcs, num_rotatable_bonds, weight_preset, batch_size, optimizer, early_stopping, patience, min_delta)`

Gradient-based torsion optimization.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mol` | `Chem.Mol` | Required | Query molecule |
| `ref_indices` | `List[int]` | Required | MCS atom indices (from query) |
| `init_coords` | `Tensor` | Required | Initial coordinates `[N, 3]` or `[B, N, 3]` |
| `pocket_coords` | `Tensor` | Required | Pocket coordinates |
| `query_features` | `dict` | Required | Query Vina features |
| `pocket_features` | `dict` | Required | Pocket Vina features |
| `num_steps` | `int` | `100` | Optimization steps |
| `lr` | `float` | `0.05` | Learning rate |
| `freeze_mcs` | `bool` | `True` | Keep MCS atoms fixed |
| `num_rotatable_bonds` | `int` | `None` | For torsional penalty |
| `weight_preset` | `str` | `"vina"` | Scoring weights |
| `batch_size` | `int` | `8` | Batch size for multi-pose optimization |
| `optimizer` | `str` | `"adam"` | Optimizer: `"adam"`, `"adamw"`, `"lbfgs"` |
| `early_stopping` | `bool` | `True` | Enable early stopping |
| `patience` | `int` | `30` | Steps without improvement |
| `min_delta` | `float` | `1e-5` | Minimum improvement threshold |

**Returns:**
- `optimized_coords`: Same shape as `init_coords`

**Optimizer Comparison:**

| Optimizer | Speed | Accuracy | Convergence | Best For |
|-----------|-------|----------|-------------|----------|
| **adam** | ⚡⚡⚡ Fast | ✓ Good | Fast, smooth | Default, speed |
| **adamw** | ⚡⚡⚡ Fast | ✓ Good | Stable | Large molecules |
| **lbfgs** | ⚡ Slower (2-3×) | ✓✓ Better | Best | High accuracy |

**Example - Single Pose:**
```python
# Optimize one pose
query_indices = [m[1] for m in mapping]  # MCS indices from query

optimized = aligner.step6_refine_pose(
    query_mol,
    ref_indices=query_indices,
    init_coords=aligned_coords[best_idx],  # [N_atoms, 3]
    pocket_coords=pocket_coords,
    query_features=query_features,
    pocket_features=pocket_features,
    num_steps=100,
    lr=0.05,
    optimizer="adam"
)
# optimized.shape: [N_atoms, 3]
```

**Example - Batch Optimization:**
```python
# Optimize all cluster representatives
optimized = aligner.step6_refine_pose(
    query_mol,
    ref_indices=query_indices,
    init_coords=aligned_coords,  # [N_poses, N_atoms, 3]
    pocket_coords=pocket_coords,
    query_features=query_features,
    pocket_features=pocket_features,
    num_steps=100,
    lr=0.05,
    batch_size=8,       # Process 8 poses at a time
    optimizer="lbfgs",
    early_stopping=True,
    patience=30
)
# optimized.shape: [N_poses, N_atoms, 3]
```

**With Different Optimizers:**
```python
# Fast but good
optimized = aligner.step6_refine_pose(..., optimizer="adam", num_steps=100)

# More stable
optimized = aligner.step6_refine_pose(..., optimizer="adamw", num_steps=100)

# Best accuracy (slower)
optimized = aligner.step6_refine_pose(..., optimizer="lbfgs", num_steps=50)
```

**Free MCS (allow MCS to move):**
```python
optimized = aligner.step6_refine_pose(
    query_mol,
    ref_indices=query_indices,
    init_coords=aligned_coords,
    pocket_coords=pocket_coords,
    query_features=query_features,
    pocket_features=pocket_features,
    freeze_mcs=False  # Allow MCS to optimize
)
```

**Early Stopping:**
```python
optimized = aligner.step6_refine_pose(
    ...,
    early_stopping=True,
    patience=30,      # Stop if no improvement for 30 steps
    min_delta=1e-5    # Minimum improvement threshold
)
```

---

## Complete Pipeline Example

```python
from lig_align import LigandAligner
from rdkit import Chem
from rdkit.Chem import AllChem
import torch

# Initialize
aligner = LigandAligner(device='cuda')

# Load molecules
ref_suppl = Chem.SDMolSupplier("ref_ligand.sdf")
ref_mol = ref_suppl[0]

query_mol = Chem.MolFromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
query_mol = Chem.AddHs(query_mol)

pocket_mol = Chem.MolFromPDBFile("pocket.pdb", sanitize=False, removeHs=True)

# Step 1: Generate conformers
query_mol, cids = aligner.step1_generate_conformers(
    query_mol,
    num_confs=1000,
    rmsd_threshold=1.0
)
print(f"Generated {len(cids)} representatives")

# Step 2: Find MCS
mapping = aligner.step2_find_mcs(ref_mol, query_mol)
print(f"MCS size: {len(mapping)} atoms")

# Step 3: Align
ref_coords = torch.tensor(ref_mol.GetConformer().GetPositions())
query_coords = torch.zeros((len(cids), query_mol.GetNumAtoms(), 3))
for i, cid in enumerate(cids):
    query_coords[i] = torch.tensor(query_mol.GetConformer(cid).GetPositions())

aligned = aligner.step3_batched_kabsch_alignment(
    ref_coords,
    query_coords,
    mapping
)

# Step 4: Score
pocket_coords = torch.tensor(pocket_mol.GetConformer().GetPositions())
query_features = aligner.compute_vina_features(query_mol)
pocket_features = aligner.compute_vina_features(pocket_mol)

from lig_align.scoring import compute_intramolecular_mask
intra_mask = compute_intramolecular_mask(query_mol, aligner.device)

scores = aligner.step4_vina_scoring(
    aligned,
    pocket_coords,
    query_features,
    pocket_features,
    weight_preset="vina",
    intramolecular_mask=intra_mask
)

print(f"Best initial score: {torch.min(scores):.2f} kcal/mol")

# Step 6: Optimize (all poses)
query_indices = [m[1] for m in mapping]

optimized = aligner.step6_refine_pose(
    query_mol,
    ref_indices=query_indices,
    init_coords=aligned,
    pocket_coords=pocket_coords,
    query_features=query_features,
    pocket_features=pocket_features,
    num_steps=100,
    lr=0.05,
    batch_size=8,
    optimizer="lbfgs",
    early_stopping=True
)

# Rescore optimized
new_scores = aligner.step4_vina_scoring(
    optimized,
    pocket_coords,
    query_features,
    pocket_features,
    intramolecular_mask=intra_mask
)

print(f"Best optimized score: {torch.min(new_scores):.2f} kcal/mol")
print(f"Average improvement: {(new_scores - scores).mean():.2f} kcal/mol")

# Step 5: Save
aligner.step5_final_selection(
    query_mol,
    cids,
    optimized,
    new_scores,
    top_k=None,
    output_path="optimized_poses.sdf"
)

print(f"Saved {len(cids)} poses to optimized_poses.sdf")
```

---

## Summary Table

| Step | Function | Key Parameters | Returns |
|------|----------|----------------|---------|
| **0** | `LigandAligner()` | `device` | aligner object |
| **1** | `step1_generate_conformers()` | `num_confs`, `rmsd_threshold`, `coordMap` | mol, cids |
| **2** | `step2_find_mcs()` | `return_all_positions`, `cross_match` | mapping(s) |
| **3** | `step3_batched_kabsch_alignment()` | `ref_coords`, `query_coords`, `mapping` | aligned coords |
| **4** | `step4_vina_scoring()` | `weight_preset`, `intramolecular_mask` | scores |
| **5** | `step5_final_selection()` | `top_k`, `output_path` | indices |
| **6** | `step6_refine_pose()` | `optimizer`, `batch_size`, `freeze_mcs` | optimized coords |

---

## Tips & Best Practices

### Conformer Generation
- **Default**: `num_confs=1000`, `rmsd_threshold=1.0`
- **Fast screening**: `num_confs=100`, `rmsd_threshold=2.0`
- **High accuracy**: `num_confs=3000`, `rmsd_threshold=0.5`

### MCS Search
- **Use Mode 1** for most cases (fastest)
- **Use Mode 2** for symmetric references
- **Use Mode 3** for complex multi-fragment matching

### Scoring
- **Use "vina"** for general purpose
- **Use "vina_lp"** if better results
- **Always use intramolecular_mask** for accuracy

### Optimization
- **Use "adam"** for speed (default)
- **Use "lbfgs"** for accuracy (slower)
- **Use "adamw"** for stability
- **Enable early_stopping** to save time
- **Adjust batch_size** based on GPU memory

### General
- Always add hydrogens before step 1
- Use GPU if available (20-50× faster)
- Monitor `num_representatives` / `num_confs` ratio
- Rescore after optimization
