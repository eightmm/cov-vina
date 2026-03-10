# CovVina Documentation

User-facing documentation for CovVina.

## Quick Links

- **[Installation & Quick Start](../README.md)**: Get started in 5 minutes
- **[Usage Guide](USAGE.md)**: Detailed usage examples and workflows
- **[API Reference](API_REFERENCE.md)**: Complete Python API documentation
- **[Architecture](ARCHITECTURE.md)**: Technical implementation details

## For Users

### Getting Started

1. **[Installation](../README.md#installation)**: Install CovVina with uv or pip
2. **[Quick Start](../README.md#quick-start)**: Run your first covalent docking in 2 minutes
3. **[Usage Guide](USAGE.md)**: Comprehensive usage examples

### Common Workflows

- **Basic covalent docking**: [USAGE.md](USAGE.md)
- **Custom parameters**: [README.md - Advanced Options](../README.md#advanced-options)
- **Visualization**: [README.md - Visualization](../README.md#visualization)

### API Documentation

- **[Python API Reference](API_REFERENCE.md)**: Complete function signatures and parameters
- **Command-line interface**: See `--help` for each script

## For Developers

### Technical Documentation

- **[Architecture](ARCHITECTURE.md)**: Pipeline stages, scoring, optimization
- **[Results & Benchmarks](../reports/results.md)**: Performance metrics and validation
- **[Development Progress](../reports/PROGRESS.md)**: Change history and design decisions

### Code Organization

```
src/cov_vina/
├── molecular/        # Molecule manipulation (adduct, anchor, conformer)
├── alignment/        # Kinematics for torsion optimization
├── scoring/          # Vina scoring and masks
├── optimization/     # Gradient-based refinement
├── io/               # Input/output utilities
└── pipeline.py       # Main orchestration
```

## Documentation Structure

| File | Purpose | Audience |
|------|---------|----------|
| [USAGE.md](USAGE.md) | Detailed usage guide | Users |
| [API_REFERENCE.md](API_REFERENCE.md) | Python API documentation | Users/Developers |
| [ARCHITECTURE.md](ARCHITECTURE.md) | Technical implementation | Developers |
| [../reports/results.md](../reports/results.md) | Benchmarks and validation | Researchers |
| [../reports/PROGRESS.md](../reports/PROGRESS.md) | Development history | Contributors |

## Reading Order

- **New user**: [../README.md](../README.md) → [USAGE.md](USAGE.md)
- **Developer**: [../README.md](../README.md) → [ARCHITECTURE.md](ARCHITECTURE.md) → [API_REFERENCE.md](API_REFERENCE.md)
- **Researcher**: [../reports/results.md](../reports/results.md) → [ARCHITECTURE.md](ARCHITECTURE.md)
