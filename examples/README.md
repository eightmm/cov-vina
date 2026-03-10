# CovVina Examples

Pre-configured examples for testing CovVina covalent docking.

## Available Systems

### 6lu7 - SARS-CoV-2 Main Protease (Mpro)

**Target**: CYS145
**PDB**: 6LU7
**Organism**: SARS-CoV-2
**Function**: Main protease (Mpro, 3CLpro)

**Files**:
- `6lu7.pdb`: Full protein structure
- `6lu7_pocket.pdb`: 12Å pocket around CYS145
- `6lu7_ligand.sdf`: Crystal ligand (for reference)
- `PJE_ideal.sdf`: Idealized crystal ligand geometry
- `results/`: Docking results for test ligands
- `redocking/`: Redocking benchmarks
- `visualizations/`: Trajectory GIFs

#### Quick Test

```bash
cd examples/6lu7

# Basic docking
uv run python ../../scripts/run_covalent_pipeline.py \
  -p 6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o test_output/

# With optimization
uv run python ../../scripts/run_covalent_pipeline.py \
  -p 6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o test_output/ \
  --optimize \
  -n 500 \
  --opt_steps 200
```

#### Expected Results

- **Runtime**: ~1-2s (200 conformers), ~2-3s (500 conformers)
- **Best score**: -0.4 to -0.5 kcal/mol (optimized)
- **Warhead**: acrylamide
- **Representative poses**: 5-15 (depends on conformer count)

### 3poz - Cathepsin K

**Target**: CYS25
**PDB**: 3POZ
**Organism**: Human
**Function**: Lysosomal cysteine protease

**Files**:
- `3poz.pdb`: Full protein structure

**Usage**:
```bash
# Extract pocket first
uv run python scripts/extract_and_redock_crystal.py \
  --pdb examples/3poz/3poz.pdb \
  --residue CYS25
```

### 1m17 - Papain

**Target**: CYS25
**PDB**: 1M17
**Organism**: Carica papaya (papaya)
**Function**: Cysteine protease

**Files**:
- `1m17.pdb`: Full protein structure

**Usage**:
```bash
# Extract pocket first
uv run python scripts/extract_and_redock_crystal.py \
  --pdb examples/1m17/1m17.pdb \
  --residue CYS25
```

## Test Molecules

`test_molecules.txt` contains various warhead types for testing:

| Molecule | Warhead Type | Mechanism | SMILES |
|----------|--------------|-----------|--------|
| Vinyl-alanine | acrylamide | Michael addition | `C=CC(=O)N[C@@H](C)C(=O)O` |
| Acrylamide-phenyl | acrylamide | Michael addition | `C=CC(=O)Nc1ccccc1` |
| Chloroacetamide | chloroacetamide | SN2 displacement | `ClCC(=O)Nc1ccccc1` |
| Vinyl sulfonamide | vinyl_sulfonamide | Michael addition | `C=CS(=O)(=O)Nc1ccccc1` |
| Bromoacetamide | bromoacetamide | SN2 displacement | `BrCC(=O)Nc1ccccc1` |
| Acrylamide-methoxy | acrylamide | Michael addition | `C=CC(=O)Nc1ccc(OC)cc1` |

## Benchmark Results

See [../reports/results.md](../reports/results.md) for detailed benchmarking data.

### 6lu7 Summary

| Ligand | Warhead | Runtime | Best Score | RMSD | Poses |
|--------|---------|---------|------------|------|-------|
| Vinyl-alanine | acrylamide | 1.7s | -0.475 | <1.0Å | 11 |

*Hardware: NVIDIA RTX PRO 6000, 200 conformers, 100 optimization steps*

## Visualization

Generate optimization trajectory:

```bash
cd examples/6lu7

uv run python ../../scripts/vis_covalent_opt_gif.py \
  -p 6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o trajectory.gif \
  --steps 150
```

## Directory Structure

```
examples/
├── README.md                    # This file
├── test_molecules.txt           # Test ligand SMILES
├── 6lu7/                        # SARS-CoV-2 Mpro (ready to use)
│   ├── 6lu7.pdb                 # Full protein
│   ├── 6lu7_pocket.pdb          # Extracted pocket
│   ├── 6lu7_ligand.sdf          # Crystal ligand
│   ├── PJE_ideal.sdf            # Idealized ligand
│   ├── results/                 # Docking results
│   │   └── vinyl_alanine/       # Example result
│   ├── redocking/               # Redocking benchmarks
│   └── visualizations/          # Trajectory GIFs
├── 3poz/                        # Cathepsin K
│   └── 3poz.pdb
└── 1m17/                        # Papain
    └── 1m17.pdb
```

## Adding New Examples

1. Place PDB file in new directory: `examples/my_protein/`
2. Extract pocket around reactive residue
3. Test with a simple warhead
4. Document results in `README.md`

## Citation

If you use these examples:

```bibtex
@misc{covvina_examples2026,
  title = {CovVina Example Structures},
  author = {CovVina Contributors},
  year = {2026},
  url = {https://github.com/eightmm/cov-vina/tree/main/examples}
}
```

## PDB Sources

- **6LU7**: Zhang et al. (2020) *Science* - SARS-CoV-2 Mpro
- **3POZ**: Brak et al. (2011) - Cathepsin K inhibitor
- **1M17**: Schröder et al. (1993) - Papain
