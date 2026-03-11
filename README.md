# CovVina

GPU-accelerated covalent docking with PyTorch-based Vina scoring and gradient optimization.

## Features

- **Automatic warhead detection**: 29 warhead types (acrylamide, aldehyde, epoxide, nitrile, etc.)
- **Multi-residue support**: CYS, SER, THR, TYR, LYS, HIS
- **GPU-accelerated**: PyTorch-based Vina scoring and gradient optimization
- **Fast batch docking**: ~0.24s per ligand for virtual screening

## Installation

```bash
# Clone repository
git clone https://github.com/eightmm/cov-vina.git
cd cov-vina

# Install with uv (recommended)
uv venv
source .venv/bin/activate  # or `.venv\Scripts\activate` on Windows
uv sync

# Or with pip
pip install -e .
```

**Requirements**: Python 3.12+, CUDA-capable GPU (optional but recommended)

## Quick Start

### Single Ligand Docking

```bash
# Basic covalent docking
uv run python scripts/run_covalent_pipeline.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o output/ \
  --optimize

# Auto-detect residue
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "SMILES_STRING" \
  -o output/

# Custom parameters
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "SMILES_STRING" \
  -r CYS145 \
  --num_confs 1000 \
  --steps 200 \
  --optimize
```

### Python API (Simplified - 1 or many ligands!)

```python
from cov_vina import run_batch_docking

# Single ligand - same interface!
results = run_batch_docking(
    protein_pdb="examples/6lu7/6lu7_pocket.pdb",
    ligands="C=CC(=O)NC",  # Just one SMILES
    reactive_residue="CYS145",
)

# Multiple ligands - same interface!
results = run_batch_docking(
    protein_pdb="examples/6lu7/6lu7_pocket.pdb",
    ligands=["C=CC(=O)NC", "O=CCc1ccccc1"],  # List of SMILES
    reactive_residue="CYS145",
)

# Or from .smi file - same interface!
results = run_batch_docking(
    protein_pdb="examples/6lu7/6lu7_pocket.pdb",
    ligands="ligands.smi",  # File path
    reactive_residue="CYS145",
)

# Check results
for r in results:
    if r['success']:
        print(f"{r['name']}: {r['best_score']:.3f} kcal/mol")
```

See [docs/USAGE.md](docs/USAGE.md) for CLI usage and detailed examples.

### Advanced Options (for run_covalent_pipeline)

```python
from cov_vina import run_covalent_pipeline

result = run_covalent_pipeline(
    # Required
    protein_pdb="pocket.pdb",           # Protein PDB file
    query_ligand="SMILES_OR_SDF",       # SMILES string or SDF file path

    # Anchor residue
    reactive_residue="CYS145",          # e.g. "CYS145", "CYS145:A", or None (auto-detect)

    # Pocket extraction
    pocket_cutoff=12.0,                 # Pocket radius (Å) - 10Å for small, 15-20Å for peptides

    # Conformer generation
    num_confs=1000,                     # Number of conformers to generate (default: 1000)
    rmsd_threshold=1.0,                 # Butina clustering RMSD threshold (Å)

    # Force field relaxation
    mmff_optimize=True,                 # Apply MMFF94 relaxation (default: True)

    # Gradient optimization
    optimize=False,                     # Enable PyTorch gradient optimization (default: False)
    optimizer="adam",                   # Optimizer: "adam", "adamw", "lbfgs" (default: "adam")
    opt_steps=100,                      # Optimization steps (default: 100)
    opt_lr=0.05,                        # Learning rate (default: 0.05 for adam/adamw, 1.0 for lbfgs)
    opt_batch_size=128,                 # Batch size for GPU optimization (default: 128)

    # Scoring
    weight_preset="vina",               # Vina weights: "vina", "vinardo", "vina_lp" (default: "vina")
    torsion_penalty=True,               # Include torsional entropy penalty (default: True)

    # Output
    output_dir="output/",               # Output directory
    save_all_poses=None,                # Save all poses (True) or top-k (False/None)
    top_k=None,                         # Number of top poses to save (None = save all)

    # Runtime
    device=None,                        # "cuda", "cpu", or None (auto-detect)
    verbose=True,                       # Print progress messages (default: True)
)

# Result dictionary
print(f"Output: {result['output_file']}")
print(f"Poses: {result['num_poses']}")
print(f"Best score: {result['best_score']:.3f} kcal/mol")
print(f"Warhead: {result['warhead_type']}")
print(f"Anchor: {result['anchor_residue']}")
print(f"Runtime: {result['runtime']:.2f}s")
```

**Key Parameters:**

- **pocket_cutoff**: 12Å default covers Vina's ~10Å scoring range. Use 10Å for small molecules, 15-20Å for large peptides.
- **num_confs**: More conformers = better sampling but slower. 200-500 for quick tests, 1000+ for publication.
- **optimizer**: Adam (fast, stable), LBFGS (slower, better convergence), AdamW (with weight decay).
- **opt_batch_size**: GPU processes this many poses simultaneously. Reduce if GPU memory limited.
- **torsion_penalty**: Vina's entropy term. Disable for interaction-only scoring (not recommended).
- **Early stopping**: Automatically enabled (patience=30, min_delta=1e-5). Stops when no improvement.

## Example: SARS-CoV-2 Mpro

Dock vinyl-alanine to SARS-CoV-2 main protease (CYS145):

```bash
uv run python scripts/run_covalent_pipeline.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o examples/6lu7/output/ \
  --optimize \
  --num_confs 1000 \
  --steps 200
```

**Expected output** (~2 seconds on RTX PRO 6000):
```
Using device: cuda
Loading protein from examples/6lu7/6lu7_pocket.pdb...
Anchor: CYS145:A atom SG (bond length 1.82 Å)
Warhead: acrylamide (reactive atom idx 0)
Generated 7 representative conformers
Best score: -0.475 kcal/mol
✓ Saved 7 poses to output/covalent_poses_all.sdf
```

## Supported Warheads

| Category | Warhead Types |
|----------|---------------|
| **Michael acceptors** | acrylamide, acrylic_acid, acrylate, enone, vinyl_sulfonamide, vinyl_sulfone, maleimide |
| **Alkylating agents** | chloroacetamide, bromoacetamide, iodoacetamide, fluoroacetamide |
| **Ring opening** | epoxide, aziridine, thiirane |
| **Nitriles** | aryl_nitrile, alkyl_nitrile, propiolamide, propargylamide |
| **SER-specific** | boronic_acid, phosphonate |
| **Others** | disulfide, sulfonyl_fluoride, aldehyde, isothiocyanate, acyl_fluoride |

## Supported Residues

| Residue | Nucleophile | Bond Length | Primary Targets |
|---------|-------------|-------------|-----------------|
| **CYS** | SG (sulfur) | 1.82 Å | Michael acceptors, alkylating agents |
| **SER** | OG (oxygen) | 1.43 Å | Boronic acids, phosphonates |
| **LYS** | NZ (nitrogen) | 1.47 Å | Activated esters, sulfonyl fluorides |
| **THR** | OG1 (oxygen) | 1.43 Å | Similar to SER |
| **TYR** | OH (oxygen) | 1.43 Å | Similar to SER |
| **HIS** | NE2 (nitrogen) | 1.47 Å | Similar to LYS |

## Output Files

SDF files contain metadata:

```python
# Covalent docking metadata
mol.GetProp("CovVina_Warhead_Type")        # "acrylamide"
mol.GetProp("CovVina_Reactive_Atom_Idx")   # "0"
mol.GetProp("CovVina_Anchor_Residue")      # "CYS145:A"
mol.GetProp("CovVina_Bond_Length")         # "1.82"

# Scoring metadata
mol.GetProp("Vina_Score_Initial")          # "-0.125"
mol.GetProp("Vina_Score_Final")            # "-0.475"
mol.GetProp("Vina_Score_Delta")            # "-0.350"
mol.GetProp("Rank")                        # "1"
```

## Visualization

Generate optimization trajectory GIF:

```bash
uv run python scripts/vis_covalent_opt_gif.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o trajectory.gif \
  --steps 200
```

## Documentation

- **[docs/USAGE.md](docs/USAGE.md)**: Detailed usage guide and workflows
- **[docs/API_REFERENCE.md](docs/API_REFERENCE.md)**: Complete API documentation
- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)**: Technical architecture
- **[reports/](reports/)**: Results, benchmarks, and progress reports

## Performance

**Hardware**: NVIDIA RTX PRO 6000 (98GB VRAM), AMD Ryzen 9 7950X

**Single ligand** (first run, includes pocket loading):

| Task | Time | GPU Usage |
|------|------|-----------|
| Warhead detection | <0.01s | CPU |
| Conformer generation (1000) | 0.8s | CPU |
| Butina clustering | 0.2s | GPU |
| Vina scoring (7 poses) | 0.01s | GPU |
| Gradient optimization (200 steps) | 0.9s | GPU |
| **Total (full pipeline)** | **~2.0s** | Mixed |

**Batch docking** (pocket cached, ~0.24s per ligand):
- 10 ligands: 2.4s total
- 100 ligands: 24s total

## Citation

If you use CovVina in your research, please cite:

```bibtex
@software{covvina2026,
  title = {CovVina: GPU-Accelerated Covalent Docking},
  author = {CovVina Contributors},
  year = {2026},
  url = {https://github.com/eightmm/cov-vina}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- Vina scoring implementation inspired by AutoDock Vina
- RDKit for molecular manipulation
- PyTorch for GPU acceleration
