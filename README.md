# CovVina

GPU-accelerated covalent docking with PyTorch-based Vina scoring and gradient optimization.

## Features

- **Automatic warhead detection**: 28 warhead types (acrylamide, epoxide, nitrile, boronic acid, etc.)
- **Multi-residue support**: CYS, SER, THR, TYR, LYS, HIS
- **Adduct-first approach**: CB-S topology created before conformer generation
- **GPU-accelerated**: PyTorch-based scoring and optimization
- **Smart clustering**: Butina clustering on full adduct structure
- **Flexible optimization**: Gradient-based torsion refinement

## Installation

```bash
# Clone repository
git clone https://github.com/your-org/cov-vina.git
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

### Command Line

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

### Python API

```python
from cov_vina import run_covalent_pipeline

# Basic usage
results = run_covalent_pipeline(
    protein_pdb="examples/6lu7/6lu7_pocket.pdb",
    query_ligand="C=CC(=O)N[C@@H](C)C(=O)O",  # vinyl-alanine
    reactive_residue="CYS145",
    output_dir="output/",
    optimize=True,
)

# Check results
print(f"Best score: {results['best_score']:.3f} kcal/mol")
print(f"Warhead: {results['warhead_type']}")
print(f"Anchor: {results['anchor_residue']}")
print(f"Output: {results['output_file']}")
```

### Advanced Options

```python
results = run_covalent_pipeline(
    protein_pdb="pocket.pdb",
    query_ligand="SMILES",
    reactive_residue="CYS145",  # or None for auto-detect

    # Conformer generation
    num_confs=1000,              # Number of conformers
    rmsd_threshold=1.0,          # Butina clustering threshold (Å)

    # Optimization
    optimize=True,               # Enable gradient optimization
    opt_steps=200,               # Optimization steps
    opt_lr=0.05,                 # Learning rate
    optimizer="adam",            # adam, adamw, or lbfgs

    # Scoring
    weight_preset="vina",        # vina, vinardo, or vina_lp

    # Output
    output_dir="output/",
    save_all_poses=False,        # Save all or top-k only
    top_k=3,                     # Number of poses to save
    verbose=True,
)
```

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
✓ Saved 7 poses to examples/6lu7/output/covalent_poses_all.sdf
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

| Task | Time | GPU Usage |
|------|------|-----------|
| Warhead detection | <0.01s | CPU |
| Conformer generation (1000) | 0.8s | CPU |
| Butina clustering | 0.2s | GPU |
| Vina scoring (7 poses) | 0.01s | GPU |
| Gradient optimization (200 steps) | 0.9s | GPU |
| **Total (full pipeline)** | **~2.0s** | Mixed |

## Citation

If you use CovVina in your research, please cite:

```bibtex
@software{covvina2026,
  title = {CovVina: GPU-Accelerated Covalent Docking},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/your-org/cov-vina}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- Vina scoring implementation inspired by AutoDock Vina
- RDKit for molecular manipulation
- PyTorch for GPU acceleration
