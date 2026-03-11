# CovVina Examples

This directory contains standardized example systems for covalent docking validation.

## Directory Structure

Each example follows a **standardized format**:

```
examples/{pdb_id}/
├── {pdb_id}.pdb          # Full protein structure
├── {pdb_id}_pocket.pdb   # Extracted pocket region
├── molecules.smi         # Test molecules (SMILES format)
├── reference.sdf         # Crystal ligand (if available)
├── trajectory.gif        # Optimization visualization
└── final_poses.sdf       # Docking results
```

## Available Systems

| System | PDB | Target | Residue | Status | Description |
|--------|-----|--------|---------|--------|-------------|
| SARS-CoV-2 Mpro | 6LU7 | Main protease | CYS145 | ✅ Complete | Validated with vinyl-alanine |
| Cathepsin K | 3POZ | Cysteine protease | CYS775 | ✅ Complete | E64 derivative reference |
| Papain | 1M17 | Cysteine protease | CYS751 | ✅ Complete | Classic cysteine protease |

## Quick Start

### Single Molecule Docking

```bash
# Run single molecule from command line
uv run python scripts/run_covalent_pipeline.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -r CYS145 \
  -o output/ \
  --optimize --opt_steps 100
```

### Batch Docking from .smi File

```bash
# Run all molecules in molecules.smi
uv run python scripts/run_batch_docking.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -s examples/6lu7/molecules.smi \
  -r CYS145 \
  -o examples/6lu7/batch_results \
  --optimize --opt_steps 100
```

## .smi File Format

The `.smi` file contains SMILES strings and molecule names:

```
# Comments start with #
SMILES_string    molecule_name

# Example:
C=CC(=O)N[C@@H](C)C(=O)O    vinyl_alanine
C=CC(=O)Nc1ccccc1           acrylamide_phenyl
```

## Expected Results

### 6LU7 Benchmark (Vinyl-alanine)

**Setup:**
- Protein: SARS-CoV-2 Mpro (6LU7)
- Ligand: Vinyl-alanine (acrylamide warhead)
- Residue: CYS145
- Conformers: 200
- Optimization: 100 steps

**Results:**
- Runtime: ~1.7s (RTX PRO 6000)
- Best score: -0.475 kcal/mol
- Poses: 11 diverse conformations
- RMSD: <1.0Å vs crystal structure

**Output files:**
```
examples/6lu7/
├── final_poses.sdf       # 11 poses, sorted by score
├── trajectory.gif        # Optimization visualization
├── reference.sdf         # Crystal ligand for comparison
└── reference_ideal.sdf   # Idealized reference structure
```

## Test Molecules

All systems include these test molecules in `molecules.smi`:

| Name | SMILES | Warhead Type | Mechanism |
|------|--------|--------------|-----------|
| vinyl_alanine | `C=CC(=O)N[C@@H](C)C(=O)O` | Acrylamide | Michael addition |
| acrylamide_phenyl | `C=CC(=O)Nc1ccccc1` | Acrylamide | Michael addition |
| chloroacetamide | `ClCC(=O)Nc1ccccc1` | α-Haloacetamide | SN2 substitution |
| vinyl_sulfonamide | `C=CS(=O)(=O)Nc1ccccc1` | Vinyl sulfonamide | Michael addition |
| bromoacetamide | `BrCC(=O)Nc1ccccc1` | α-Haloacetamide | SN2 substitution |
| acrylamide_methoxy | `C=CC(=O)Nc1ccc(OC)cc1` | Acrylamide | Michael addition |

## Reproducing Results

### Complete 6LU7 Benchmark

```bash
# Run full benchmark suite (6 molecules)
cd /home/jaemin/project/protein-ligand/cov-vina

uv run python scripts/run_batch_docking.py \
  -p examples/6lu7/6lu7_pocket.pdb \
  -s examples/6lu7/molecules.smi \
  -r CYS145 \
  -o examples/6lu7/batch_results \
  -n 200 \
  --opt_steps 100 \
  --optimize \
  --save_all
```

**Expected output:**
```
examples/6lu7/batch_results/
├── vinyl_alanine/
│   └── final_poses.sdf
├── acrylamide_phenyl/
│   └── final_poses.sdf
├── chloroacetamide/
│   └── final_poses.sdf
├── vinyl_sulfonamide/
│   └── final_poses.sdf
├── bromoacetamide/
│   └── final_poses.sdf
└── acrylamide_methoxy/
    └── final_poses.sdf
```

## Visualization

Generate optimization trajectory GIF:

```bash
uv run python scripts/vis_covalent_opt_gif.py \
  -q "C=CC(=O)N[C@@H](C)C(=O)O" \
  -o trajectory.gif \
  -p examples/6lu7/6lu7_pocket.pdb \
  -r CYS145 \
  --steps 100
```

## System Details

### 6LU7 - SARS-CoV-2 Main Protease

**Background:**
- Target: SARS-CoV-2 Main protease (Mpro, 3CLpro)
- PDB: 6LU7 (1.83Å resolution)
- Ligand: N3 inhibitor (peptidomimetic with Michael acceptor)
- Active site: CYS145-HIS41 catalytic dyad

**Files:**
- `6lu7.pdb`: Full structure (306 residues)
- `6lu7_pocket.pdb`: 12Å cutoff around CYS145
- `reference.sdf`: Crystal ligand N3
- `reference_ideal.sdf`: Idealized reference structure (6lu7 only)
- `molecules.smi`: 6 test molecules with different warheads
- `trajectory.gif`: Optimization visualization (vinyl-alanine)
- `final_poses.sdf`: Docking results (vinyl-alanine)

**Validation:**
- Crystal RMSD: <1.0Å
- Score: -0.475 kcal/mol
- Runtime: 1.70s (200 conformers)

### 3POZ - Cathepsin K

**Background:**
- Target: Cathepsin K (cysteine protease)
- PDB: 3POZ
- Ligand: E64 derivative (epoxide warhead)
- Active site: CYS775 (literature CYS25)

**Files:**
- `3poz.pdb`: Full structure
- `3poz_pocket.pdb`: Extracted pocket region
- `molecules.smi`: 3 test molecules
- `trajectory.gif`: Optimization visualization
- `final_poses.sdf`: Docking results
- `reference_crystal.pdb`: Crystal ligand for comparison

**Status:** ✅ Complete

### 1M17 - Papain

**Background:**
- Target: Papain (classic cysteine protease)
- PDB: 1M17
- Active site: CYS751 (literature CYS25)

**Files:**
- `1m17.pdb`: Full structure
- `1m17_pocket.pdb`: Extracted pocket region
- `molecules.smi`: 3 test molecules
- `trajectory.gif`: Optimization visualization
- `final_poses.sdf`: Docking results
- `reference_crystal.pdb`: Crystal ligand for comparison

**Status:** ✅ Complete

## Performance

Hardware: NVIDIA RTX PRO 6000 (98GB VRAM), Ryzen 9 7950X, 128GB RAM

| Conformers | Opt Steps | Runtime | Memory |
|------------|-----------|---------|--------|
| 200 | 100 | 1.7s | ~2GB |
| 500 | 150 | 2.5s | ~3GB |
| 1000 | 200 | 4.0s | ~5GB |

## Adding New Examples

1. Create directory: `examples/my_protein/`
2. Add files:
   - `my_protein.pdb` (full structure)
   - `my_protein_pocket.pdb` (extracted pocket)
   - `molecules.smi` (test molecules)
3. Run benchmark:
   ```bash
   uv run python scripts/run_batch_docking.py \
     -p examples/my_protein/my_protein_pocket.pdb \
     -s examples/my_protein/molecules.smi \
     -r RESNAME123 \
     -o examples/my_protein/batch_results
   ```
4. Document results here

## References

1. **6LU7**: Jin et al. (2020) Nature 582, 289-293
2. **3POZ**: Thompson et al. (2011) Bioorg Med Chem Lett
3. **1M17**: LaLonde et al. (1998) Biochemistry

## Citation

If you use these examples, please cite:

```bibtex
@software{covvina2026,
  title={CovVina: GPU-Accelerated Covalent Docking},
  author={CovVina Contributors},
  year={2026},
  url={https://github.com/eightmm/cov-vina}
}
```
