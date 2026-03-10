# Usage Guide

This document holds practical setup and execution details that are too verbose for the top-level README.

## Installation

```bash
uv venv
source .venv/bin/activate
uv sync
```

If `uv sync` is not suitable for your environment, install the dependencies declared in [pyproject.toml](/home/jaemin/project/protein-ligand/cov-vina/pyproject.toml).

## Covalent Docking CLI

Entry point for covalent docking:

```bash
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "C=CC(=O)Nc1ccccc1" \
  -r CYS145 \
  -o output/
```

Required arguments:

```text
-p, --protein            Protein pocket PDB file
-q, --query_ligand       Query ligand as SMILES or SDF path (must contain a reactive warhead)
```

Covalent-specific arguments:

```text
-r, --reactive_residue   Reactive residue, e.g. "CYS145" or "CYS145:A" (default: auto-detect)
```

Common optional arguments:

```text
-o, --out_dir             Output directory
-n, --num_confs           Number of conformers to generate
--rmsd_threshold          RMSD threshold for clustering
--optimize                Enable gradient optimization
--opt_batch_size          Number of poses optimized together (default: 128)
--optimizer               adam | adamw | lbfgs
--weight_preset           vina | vina_lp | vinardo
--free_anchor             Allow anchor atom to move during optimization
--no_torsion_penalty      Disable the default Vina torsional entropy penalty
```

### Covalent Docking Workflows

#### Basic Covalent Docking

```bash
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "C=CC(=O)Nc1ccccc1" \
  -r CYS145
```

#### Optimized Covalent Docking

```bash
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "C=CC(=O)Nc1ccccc1" \
  -r CYS145 \
  --optimize \
  --optimizer lbfgs \
  --opt_steps 200
```

#### Auto-Detect Reactive Residue

```bash
# Omit -r to auto-detect the first supported reactive residue in the pocket
uv run python scripts/run_covalent_pipeline.py \
  -p pocket.pdb \
  -q "C=CC(=O)Nc1ccccc1"
```

#### Supported Warhead Types

The pipeline auto-detects warheads. Currently supported:

- **Michael acceptors**: acrylamide, acrylate, enone, vinyl sulfonamide, vinyl sulfone, maleimide
- **alpha-Halo carbonyl**: chloroacetamide, bromoacetamide, iodoacetamide, fluoroacetamide, chlorofluoroacetamide
- **Strained rings**: epoxide, aziridine
- **Triple bond electrophiles**: nitrile, propiolamide, propargylamide
- **Reversible**: cyanoacrylamide
- **Others**: disulfide, sulfonyl fluoride, aldehyde

#### Supported Reactive Residues

Currently: **CYS** (SG atom, C-S bond length 1.82 Å). Extensible to SER, LYS, THR.

## MCS-Guided Docking CLI

Entry point for reference-based alignment:

```bash
uv run python scripts/run_pipeline.py \
  -p examples/10gs/10gs_pocket.pdb \
  -r examples/10gs/10gs_ligand.sdf \
  -q "CCO" \
  -o output/
```

Required arguments:

```text
-p, --protein        Protein pocket PDB file
-r, --ref_ligand     Reference ligand SDF file
-q, --query_ligand   Query ligand as SMILES or SDF path
```

Common optional arguments:

```text
-o, --out_dir             Output directory
-n, --num_confs           Number of conformers to generate
--rmsd_threshold          RMSD threshold for clustering
--mcs_mode                auto | single | multi | cross
--optimize                Enable gradient optimization
--opt_batch_size          Number of poses optimized together (default: 128)
--optimizer               adam | adamw | lbfgs
--weight_preset           vina | vina_lp | vinardo
--free_mcs                Allow MCS atoms to move during optimization
--no_torsion_penalty      Disable the default Vina torsional entropy penalty
```

### MCS Docking Workflows

#### Basic Prediction

```bash
uv run python scripts/run_pipeline.py \
  -p pocket.pdb \
  -r ref_ligand.sdf \
  -q "CCO" \
  -o output/
```

#### Recommended Optimized Run

```bash
uv run python scripts/run_pipeline.py \
  -p pocket.pdb \
  -r ref_ligand.sdf \
  -q "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" \
  -n 1000 \
  --rmsd_threshold 1.0 \
  --optimize \
  --optimizer lbfgs
```

Scoring note:

- `Vina` scores include the standard torsional entropy penalty by default
- use `--no_torsion_penalty` only when you explicitly want interaction-only scores
- `opt_batch_size=128` is now the default for same-molecule multi-pose optimization on GPU
- reduce it if your ligand leaves many representative poses or if GPU memory becomes limiting

### MCS Mode Guidance

- `auto`
  - recommended default
  - chooses `multi` for symmetry-equivalent placements
  - chooses `cross` only when multi-fragment matching increases total mapped atoms
  - otherwise uses `single`
- `single`
  - use for the fastest and most conservative contiguous-core match
- `multi`
  - use when the reference is symmetric and you want all equivalent placements enumerated
  - current pipeline still continues with the first candidate after enumeration
- `cross`
  - use when one contiguous MCS is too restrictive and multiple fragments matter
  - current pipeline still continues with the first generated combination

### Query From SDF

```bash
uv run python scripts/run_pipeline.py \
  -p protein.pdb \
  -r ref_ligand.sdf \
  -q query_ligand.sdf \
  -o output_sdf/
```

### Optimize Existing Pose Only

```bash
uv run python scripts/optimize_pose.py \
  -p protein.pdb \
  -l ligand.sdf \
  -o optimized.sdf \
  --steps 200 \
  --optimizer lbfgs
```

## Example Assets

Covalent docking example under `examples/6lu7` (SARS-CoV-2 Mpro, CYS145):

- `examples/6lu7/6lu7.pdb`: full protein structure
- `examples/6lu7/6lu7_pocket.pdb`: 20Å pocket around CYS145:SG

## Relaxation And Score Metadata

Each exported SDF can include run metadata that explains what happened during placement and optimization.

### Covalent Docking Metadata

- `CovVina_Warhead_Type`: detected warhead group (e.g. "acrylamide", "chloroacetamide")
- `CovVina_Reactive_Atom_Idx`: ligand atom index forming the covalent bond
- `CovVina_Anchor_Residue`: protein residue and chain (e.g. "CYS145:A")
- `CovVina_Bond_Length`: covalent bond length used for placement (Å)
- `CovVina_Gradient_Optimized`: whether gradient optimization was applied

### MCS Docking Metadata

- `LigAlign_MCS_Mode`: the mode actually used after `auto` resolution
- `LigAlign_MCS_Mode_Requested`: the mode requested by the user
- `LigAlign_MMFF_Requested`: whether relaxation was requested
- `LigAlign_MMFF_Optimized`: whether relaxation actually ran successfully
- `LigAlign_Relaxation_Summary`: why relaxation was applied, skipped, or fell back

### Common Metadata

- `Vina_Score_Initial`: score before gradient optimization
- `Vina_Score_Final`: score after optimization or final ranking pass
- `Vina_Score_Delta`: final minus initial score

Practical interpretation:

- if `Vina_Score_Delta` is negative, optimization improved the score
- if `Vina_Score_Delta` is near zero, either the pose was already near a local minimum or there were no useful torsional moves available

## Testing

Run all tests with pytest:

```bash
uv run python -m pytest tests/ -v
```
