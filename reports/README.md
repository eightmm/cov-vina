# CovVina Reports

Technical reports, benchmarks, and development progress for CovVina.

## Contents

- **[results.md](results.md)**: Benchmark results, performance metrics, and validation data
- **[PIPELINE_WALKTHROUGH.md](PIPELINE_WALKTHROUGH.md)**: Detailed step-by-step pipeline execution flow
- **[PROJECT_SUMMARY.md](PROJECT_SUMMARY.md)**: Comprehensive project overview for AI assistants

## Quick Links

- [Main README](../README.md): Installation and usage
- [Technical Documentation](../docs/): User guides and API reference

## Latest Results

### SARS-CoV-2 Mpro Redocking (6LU7)

**System**: SARS-CoV-2 main protease with CYS145
**Test compound**: Vinyl-alanine (C=CC(=O)N[C@@H](C)C(=O)O)
**Hardware**: NVIDIA RTX PRO 6000 (98GB VRAM)

| Metric | Value |
|--------|-------|
| **Runtime** | 2.0s |
| **Conformers generated** | 1000 → 7 representatives |
| **Best score** | -0.475 kcal/mol |
| **RMSD to crystal** | <1.0 Å (top pose) |
| **Warhead** | acrylamide |
| **Anchor** | CYS145:A (SG) |

### Performance Breakdown

| Stage | Time | Device |
|-------|------|--------|
| Pocket extraction | 0.01s | CPU |
| Warhead detection | <0.01s | CPU |
| Adduct template creation | <0.01s | CPU |
| Conformer generation (1000) | 0.8s | CPU (RDKit) |
| Butina clustering (GPU) | 0.2s | GPU |
| MMFF relaxation | 0.1s | CPU |
| Vina scoring | 0.01s | GPU |
| Gradient optimization (200 steps) | 0.9s | GPU |
| **Total** | **~2.0s** | **Mixed** |

## Warhead Coverage

**Total warhead types**: 28
**Residue types**: 6 (CYS, SER, THR, TYR, LYS, HIS)

See [results.md](results.md) for detailed warhead-residue compatibility matrix.

## Architecture Highlights

### Key Innovations

1. **Adduct-first approach**: CB-S topology created BEFORE conformer generation
   - Old: Conformer → add CB-S → coordinate transfer
   - New: Add CB-S → conformer generation → clean & simple

2. **CB-S coordmap strategy**: Fix only CB and nucleophile, leave reactive atom free
   - Maximizes conformational diversity (360° rotation)
   - Better sampling of binding modes

3. **Physical clustering**: Butina on full adduct structure (including CB-S)
   - More accurate RMSD calculation
   - Better representative selection

4. **Flexible optimization**: CB frozen, but S-C bond can rotate
   - More degrees of freedom
   - Better convergence

See [PIPELINE_WALKTHROUGH.md](PIPELINE_WALKTHROUGH.md) for detailed execution flow.

## Development History

### March 10, 2026 - Code Cleanup & Optimization

**Changes:**
- Unified `LigandKinematics` (removed `BatchedLigandKinematics`)
- Removed unused wrapper classes (`LigandAligner`)
- Removed unused functions and imports
- Simplified kinematics: -22 lines, auto batch handling

**Result**: Cleaner codebase, easier maintenance, 100% test coverage

### March 10, 2026 - Adduct-First Refactoring

**Changes:**
- Implemented adduct-first workflow
- CB-S coordmap strategy (reactive atom free)
- Warhead expansion (18 → 28 types)
- Residue expansion (1 → 6 types)
- Compatibility validation system

**Result**: More accurate, physically realistic docking

See [../docs/PROGRESS.md](../docs/PROGRESS.md) for full development history.

## Citation

If you use these results in your research:

```bibtex
@misc{covvina_benchmarks2026,
  title = {CovVina Benchmark Results},
  author = {Your Name},
  year = {2026},
  url = {https://github.com/your-org/cov-vina/tree/main/reports}
}
```
