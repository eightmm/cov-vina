# LigAlign Documentation

This directory contains technical documentation and design analysis for LigAlign.

## Main Documentation

- **[../README.md](../README.md)**: Main project README with usage instructions and API reference

## Technical Documentation

### Pharmacophore Features
- **[PHARMACOPHORE.md](PHARMACOPHORE.md)**: Comprehensive guide to pharmacophore feature detection
  - What pharmacophore features are
  - Group-level vs atom-level detection
  - Why NOT used for alignment (MCS is better)
  - Usage examples and API reference

### Step-by-Step Capabilities
- **[STEP_BY_STEP_CAPABILITIES.md](STEP_BY_STEP_CAPABILITIES.md)**: Detailed documentation of all 7 pipeline steps
  - Complete parameter reference
  - Examples for each step
  - Low-level API usage

## Design Analysis (Historical)

These documents capture design decisions and analysis during development:

- **[ALIGNMENT_METHODS_ANALYSIS.md](ALIGNMENT_METHODS_ANALYSIS.md)**: Comparison of MCS, Pharmacophore, and Shape-based alignment methods
- **[PHARMACOPHORE_DESIGN.md](PHARMACOPHORE_DESIGN.md)**: Initial design for group-level pharmacophore features
- **[FUNCTIONAL_GROUPS_COMPLETE.md](FUNCTIONAL_GROUPS_COMPLETE.md)**: Implementation details and testing results
- **[DEDUPLICATION_COMPLETE.md](DEDUPLICATION_COMPLETE.md)**: Hierarchical deduplication strategy

## Key Decisions

### 1. Alignment Method: MCS Only

**Decision**: Use Maximum Common Substructure (MCS) for alignment, NOT pharmacophore-based.

**Rationale**:
- MCS provides exact atom-to-atom mapping
- Reference structure guides binding mode
- Works with existing optimization (needs anchors)
- Pharmacophore matching can't provide precise atom mapping

See [ALIGNMENT_METHODS_ANALYSIS.md](ALIGNMENT_METHODS_ANALYSIS.md) for full analysis.

### 2. Pharmacophore: Analysis Tool Only

**Decision**: Pharmacophore feature detection for analysis/visualization, not alignment.

**Use cases**:
- Molecular quality checks
- 3D visualization
- MCS validation
- Future scoring enhancements

See [PHARMACOPHORE.md](PHARMACOPHORE.md) for usage guide.

### 3. Group-Level Features

**Decision**: Detect functional groups (not individual atoms) for cleaner pharmacophore.

**Benefits**:
- Benzene = 1 group (not 6 atoms)
- Carboxylate = 1 group (not 3 separate features)
- 89% feature reduction (2 vs 19 for ibuprofen)

See [FUNCTIONAL_GROUPS_COMPLETE.md](FUNCTIONAL_GROUPS_COMPLETE.md) for implementation.

### 4. No Aliphatic Chains

**Decision**: Exclude aliphatic chains from pharmacophore features.

**Rationale**:
- Pharmacophore = features critical for binding
- Aliphatic chains are structural linkers, not binding features
- Keeps features pharmacologically meaningful

## Quick Links

### For Users
- [Main README](../README.md) - Start here
- [Step-by-Step Guide](STEP_BY_STEP_CAPABILITIES.md) - Detailed API reference

### For Developers
- [Pharmacophore Implementation](PHARMACOPHORE.md) - Feature detection details
- [Design Decisions](ALIGNMENT_METHODS_ANALYSIS.md) - Why MCS not pharmacophore

### For Researchers
- [Methods Analysis](ALIGNMENT_METHODS_ANALYSIS.md) - Comparison of alignment approaches
- [Feature Design](PHARMACOPHORE_DESIGN.md) - Group-level pharmacophore rationale
