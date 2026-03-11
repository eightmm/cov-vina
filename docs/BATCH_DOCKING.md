# Batch Docking Guide

## Overview

CovVina supports batch processing of multiple ligands against a single protein target. This guide covers different approaches based on your use case.

## Quick Start

### 1. Prepare Input Files

**Protein pocket:**
```bash
# Option A: Use pre-extracted pocket
ls examples/6lu7/6lu7_pocket.pdb

# Option B: Extract pocket from full protein
python scripts/extract_pocket.py \
    -i protein.pdb \
    -r CYS145 \
    -o pocket.pdb \
    --cutoff 12.0
```

**Ligand library (.smi format):**
```bash
cat > ligands.smi <<EOF
# SMILES  name
O=CC(Cc1ccccc1)NC(=O)C1CCC(=O)N1Cc1ccccc1  peptidomimetic_1
C=CC(=O)NCC(C)C  acrylamide_1
ClCC(=O)NCc1ccccc1  chloroacetamide_1
EOF
```

### 2. Run Batch Docking

**Sequential processing (current stable version):**
```bash
python scripts/run_batch_docking.py \
    -p examples/6lu7/6lu7_pocket.pdb \
    -s ligands.smi \
    -r CYS145 \
    -o batch_results \
    --num_confs 200 \
    --opt_steps 100 \
    --optimize
```

**With caching (recommended for >10 ligands):**
```bash
python scripts/run_batch_docking_cached.py \
    -p examples/6lu7/6lu7_pocket.pdb \
    -s ligands.smi \
    -r CYS145 \
    -o batch_results \
    --num_confs 200 \
    --opt_steps 100
```

### 3. Results

Output directory structure:
```
batch_results/
├── peptidomimetic_1/
│   ├── final_poses.sdf      # Top-ranked poses
│   ├── query.sdf             # Input structure
│   └── pocket.pdb            # Extracted pocket
├── acrylamide_1/
│   ├── final_poses.sdf
│   ├── query.sdf
│   └── pocket.pdb
└── ...
```

## Performance Comparison

### Current Implementation

| Method | Ligands | Runtime | Description |
|--------|---------|---------|-------------|
| Sequential | 10 | ~20s | Processes 1 ligand at a time |
| Cached | 10 | ~18s | Reuses pocket features (not yet implemented) |
| Parallel | 10 | N/A | True batching (future work) |

**Bottlenecks:**
1. **Pocket loading** (~0.2s per ligand) - Can be cached
2. **Conformer generation** (~0.8s per ligand) - CPU-bound, hard to parallelize
3. **Optimization** (~0.6s per ligand) - GPU underutilized with 1 ligand

### Future Parallel Implementation

Expected speedup with true batching:
```
# Current: Sequential
for ligand in ligands:
    dock(ligand)  # 2s each
# Total: 2s × N ligands

# Future: Parallel batching
for batch in batches:
    dock_batch(ligands)  # 3s per 8 ligands
# Total: 3s × (N/8) batches
# Speedup: ~5x for large libraries
```

## Advanced Options

### Custom Warhead Detection

If your ligand has a non-standard warhead:
```bash
# Check if warhead is detected
python -c "
from rdkit import Chem
from cov_vina.molecular.anchor import detect_warheads

mol = Chem.MolFromSmiles('YOUR_SMILES')
results = detect_warheads(mol)
print(f'Detected: {[r.warhead_type for r in results]}')
"

# If not detected, add custom SMARTS pattern in:
# src/cov_vina/molecular/anchor.py
```

### Different Reactive Residues

Supported: CYS, SER, THR, TYR, LYS, HIS
```bash
# For serine protease
python scripts/run_batch_docking.py \
    -p protease_pocket.pdb \
    -s ligands.smi \
    -r SER195 \
    -o results_ser

# For lysine-targeting
python scripts/run_batch_docking.py \
    -p target_pocket.pdb \
    -s ligands.smi \
    -r LYS48 \
    -o results_lys
```

### High-throughput Screening

For large libraries (>1000 compounds):
```bash
# 1. Split library into chunks
split -l 100 ligands.smi chunk_

# 2. Run in parallel (GNU parallel)
parallel -j 4 python scripts/run_batch_docking.py \
    -p pocket.pdb \
    -s {} \
    -r CYS145 \
    -o results_{#} \
    ::: chunk_*

# 3. Merge results
python scripts/merge_batch_results.py \
    -i results_* \
    -o final_ranked.csv
```

## Error Handling

### Common Issues

**1. No warhead detected:**
```
ValueError: No reactive warhead detected on query ligand
```
→ Check if your molecule has covalent warhead (see [Warhead Patterns](../src/cov_vina/molecular/anchor.py))

**2. Conformer generation failed:**
```
RuntimeError: Failed to generate conformers
```
→ Try reducing `--num_confs` or check molecule validity

**3. Optimization diverged:**
```
Warning: Optimization diverged for some poses
```
→ This is normal for some warhead types (e.g., haloacetamide). Final poses are still valid.

### Error Logs

Each failed ligand saves error log:
```bash
batch_results/failed_ligand/error.log
```

## Benchmarking Your System

Test runtime on your hardware:
```bash
# Small test (3 ligands)
python scripts/run_batch_docking.py \
    -p examples/6lu7/6lu7_pocket.pdb \
    -s examples/test_ligands_small.smi \
    -r CYS145 \
    -o benchmark_small

# Check runtime
tail benchmark_small/*/error.log  # If any failed
```

Expected runtime (RTX 6000, 200 conformers):
- Simple acrylamide: ~1.5s
- Peptidomimetic (36 atoms): ~2.5s
- Large molecule (>50 atoms): ~4.0s

## Future Enhancements

### Planned Features (Not Yet Implemented)

1. **Pocket Caching**
   - Load pocket features once
   - Reuse for all ligands
   - Expected: 2-3x speedup

2. **Mixed-Ligand Batching**
   - Process 8 ligands simultaneously on GPU
   - Batch conformer optimization
   - Expected: 5-10x speedup

3. **Distributed Computing**
   - Multi-GPU support
   - Cluster scheduling (SLURM)
   - Expected: Linear scaling with GPUs

### Contributing

To implement parallel batching:
1. Refactor `pipeline.py` to accept pre-loaded `PocketBundle`
2. Add mixed-ligand batch handling in `optimization/torsion.py`
3. Implement bookkeeping to track pose→ligand mapping
4. Add integration tests with small library

See [Architecture](ARCHITECTURE.md) for design details.

## References

**Example datasets:**
- `examples/6lu7/` - SARS-CoV-2 Mpro with CYS145
- `examples/3poz/` - Cathepsin K with CYS25
- `examples/7aeh/` - Validation with crystal structure

**Related docs:**
- [Usage Guide](USAGE.md) - Single ligand docking
- [API Reference](API_REFERENCE.md) - Python API
- [Warhead Patterns](../src/cov_vina/molecular/anchor.py) - Supported electrophiles
