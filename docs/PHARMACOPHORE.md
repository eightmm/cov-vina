# Pharmacophore Feature Detection

**Purpose**: Analysis and visualization tool for molecular pharmacophore features.

**Note**: This module is for **analysis only**. Alignment uses **MCS (Maximum Common Substructure)** exclusively.

## Overview

This module detects pharmacophore features at the functional group level (not individual atoms) for:
- Molecular quality analysis
- 3D visualization
- MCS validation
- Future scoring improvements

## What are Pharmacophore Features?

Pharmacophore features are functional groups critical for biological binding:

| Feature Type | Importance | Examples | Description |
|-------------|-----------|----------|-------------|
| **Ionizable** | 5 (Highest) | Carboxylate, Ammonium, Sulfonate | Charged/chargeable groups |
| **H-bond Donor** | 5 | Hydroxyl (-OH), Amine (-NH2) | Donate H in H-bonds |
| **H-bond Acceptor** | 5 | Carbonyl (C=O), Ether (-O-) | Accept H in H-bonds |
| **Aromatic** | 3 | Benzene, Naphthalene, Pyridine | Aromatic ring systems |
| **Halogen** | 2 | F, Cl, Br, I | Halogen bonds |

**Excluded**: Aliphatic chains (e.g., methyl, ethyl) - these are structural linkers, NOT pharmacophore features.

## Key Design Decisions

### 1. Group-Level, Not Atom-Level

**Problem with atom-level (RDKit default):**
```python
Ibuprofen pharmacophore (RDKit):
  - 12 hydrophobe features (noise!)
  - 3 separate features for COOH (OH + C=O + O)
  - Total: ~19 noisy features
```

**Our solution (group-level):**
```python
Ibuprofen pharmacophore (ours):
  - 1 Carboxylate (3 atoms as ONE group)
  - 1 Benzene (6 atoms as ONE group)
  - Total: 2 clean features (89% reduction!)
```

### 2. Hierarchical Deduplication

**Rule 1: Ionizable > H-bond**
- If atoms belong to ionizable group (COOH), don't show as separate donor/acceptor
- Example: Carboxylate is ONE feature, not donor + 2 acceptors

**Rule 2: Aromatic Integration**
- Heteroatoms in aromatic rings stay integrated
- Example: Pyridine N is part of aromatic feature, not separate acceptor

**Rule 3: Donor > Acceptor**
- If same atom is both donor and acceptor, show as donor only
- Example: Hydroxyl (-OH) is donor, not both

### 3. No Aliphatic Chains

**Why excluded?**
- Not pharmacophore features (structural linkers only)
- "Pharmacophore" = features critical for binding
- Tert-butyl, isopropyl, etc. are hydrophobic but not distinct pharmacophore points

## Usage

### Basic Usage

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from lig_align.molecular.functional_groups import detect_pharmacophore_features

# Load molecule
mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')  # Ibuprofen
mol.SetProp('_Name', 'Ibuprofen')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)

# Detect pharmacophore features
features = detect_pharmacophore_features(mol)

# Print results
for f in features:
    print(f"{f.type}: {f.subtype} ({len(f.atoms)} atoms, importance={f.importance})")

# Output:
# NegIonizable: Carboxylate (3 atoms, importance=5)
# Aromatic: Benzene (6 atoms, importance=3)
```

### With Visualization

```python
from lig_align.molecular.functional_groups import detect_pharmacophore_features, print_functional_groups

mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
mol.SetProp('_Name', 'Ibuprofen')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)

features = detect_pharmacophore_features(mol)
print_functional_groups(mol, features)
```

### Analysis Use Cases

**1. Molecule Comparison:**
```python
ref_features = detect_pharmacophore_features(ref_mol)
query_features = detect_pharmacophore_features(query_mol)

ref_types = set(f.subtype for f in ref_features)
query_types = set(f.subtype for f in query_features)

common = ref_types & query_types
print(f"Common pharmacophore features: {common}")
# → Validate MCS found these important features
```

**2. Quality Check:**
```python
features = detect_pharmacophore_features(mol)

if len(features) < 2:
    print("Warning: Very few pharmacophore features - may bind poorly")

ionizable = [f for f in features if 'Ionizable' in f.type]
if len(ionizable) == 0:
    print("Note: No ionizable groups - weaker binding expected")
```

**3. Feature Importance Analysis:**
```python
features = detect_pharmacophore_features(mol)

high_importance = [f for f in features if f.importance >= 5]
print(f"Critical binding features: {len(high_importance)}")

for f in high_importance:
    print(f"  {f.subtype} at atoms {f.atoms}")
```

## Examples

### Ibuprofen: `CC(C)Cc1ccc(cc1)C(C)C(=O)O`

```
Pharmacophore features: 2
  NegIonizable: Carboxylate (3 atoms) [importance=5]
  Aromatic: Benzene (6 atoms) [importance=3]

Note: Tert-butyl and isopropyl chains NOT included (not pharmacophore features)
```

### Naproxen: `COc1ccc2cc(ccc2c1)C(C)C(=O)O`

```
Pharmacophore features: 3
  Acceptor: Ether (1 atom) [importance=5]
  NegIonizable: Carboxylate (3 atoms) [importance=5]
  Aromatic: Naphthalene (10 atoms) [importance=3]

Note: Naphthalene fused rings = ONE aromatic feature
```

### Aspirin: `CC(=O)Oc1ccccc1C(=O)O`

```
Pharmacophore features: 4
  Acceptor: Carbonyl (1 atom) [importance=5]  ← Ester
  Acceptor: Ether (1 atom) [importance=5]     ← Ester
  NegIonizable: Carboxylate (3 atoms) [importance=5]
  Aromatic: Benzene (6 atoms) [importance=3]

Note: Ester and carboxylate are separate (different chemical roles)
```

### Pyridine: `c1ccncc1`

```
Pharmacophore features: 1
  Aromatic: Pyridine (6 atoms) [importance=3]

Note: Nitrogen integrated into aromatic feature (not separate acceptor)
```

### Ethanol: `CCO`

```
Pharmacophore features: 1
  Donor: Hydroxyl (1 atom) [importance=5]

Note: Ethyl chain NOT included (not a pharmacophore feature)
```

## Implementation Details

### FunctionalGroup Data Structure

```python
@dataclass
class FunctionalGroup:
    type: str              # "Aromatic", "NegIonizable", "Donor", etc.
    subtype: str           # "Benzene", "Carboxylate", "Hydroxyl", etc.
    atoms: List[int]       # Atom indices [4, 5, 6, 7, 8, 9]
    centroid: np.ndarray   # 3D geometric center
    importance: int        # 1-5 (5 = highest)
```

### Detection Hierarchy

1. **H-bond & Ionizable (importance 5)**: RDKit feature factory
2. **Deduplication**: Remove overlaps (ionizable > donor/acceptor)
3. **Aromatic Systems (importance 3)**: Ring detection + fusion analysis
4. **Aromatic Integration**: Remove heteroatom features inside rings
5. **Halogens (importance 2)**: Simple element matching
6. **Sort by importance**: Descending order

### Aromatic System Detection

```python
# Benzene: 1 ring, 6 atoms → "Benzene"
# Naphthalene: 2 fused rings, 10 atoms → "Naphthalene"
# Biphenyl: 2 non-fused rings → 2 separate "Benzene" groups
# Indole: 2 fused rings (5+6), N in ring → "Indole"
```

Algorithm:
1. Find all aromatic rings
2. Merge fused rings (shared atoms)
3. Classify by ring count and heteroatoms

## Why NOT Used for Alignment

### The Fundamental Problem

**MCS (used):**
- Exact atom-to-atom mapping
- Reference structure guides alignment
- Works with existing optimization

**Pharmacophore (not used):**
- Group-level matching only
- Benzene (6 atoms) ↔ Naphthalene (10 atoms)? Which atoms map?
- Centroid-only alignment loses anchor atoms
- Can't optimize without anchors

### Example Issue

```
Reference: Ibuprofen
  - Benzene at specific position (defines binding mode)

Query: Naproxen
  - Naphthalene (2 rings)

Pharmacophore matching:
  - Match "aromatic" features ✓
  - But which of the 2 naphthalene rings aligns to benzene? ❌
  - Centroid is between rings - wrong! ❌

MCS matching:
  - Find 6-atom benzene substructure in naphthalene ✓
  - Exact atom mapping ✓
  - Uses reference position info ✓
```

**Conclusion**: MCS is superior for template-based alignment.

## Comparison Table

| Aspect | Atom-Level (RDKit) | Our Group-Level | Notes |
|--------|-------------------|----------------|-------|
| Ibuprofen features | 19 noisy | 2 clean | 89% reduction |
| Benzene | 6 hydrophobes | 1 aromatic | Meaningful unit |
| Carboxylate | 3 separate | 1 ionizable | No duplicates |
| Pyridine N | Aromatic + Acceptor | Aromatic only | Integrated |
| Aliphatic | Counted | Excluded | Not pharmacophore |
| Usability | Low (noise) | High (clear) | Analysis-ready |

## API Reference

### Main Functions

```python
def detect_pharmacophore_features(mol: Chem.Mol) -> List[FunctionalGroup]:
    """
    Detect pharmacophore features (recommended).

    Returns:
        List of FunctionalGroup objects sorted by importance (descending)
    """

def detect_functional_groups(mol: Chem.Mol, include_aliphatic: bool = False) -> List[FunctionalGroup]:
    """
    Backward compatibility wrapper.

    Use detect_pharmacophore_features() instead.
    """

def print_functional_groups(mol: Chem.Mol, groups: List[FunctionalGroup]):
    """
    Pretty-print pharmacophore features grouped by type.
    """
```

### Helper Functions

```python
detect_hbond_ionizable(mol)     # RDKit feature factory
detect_aromatic_systems(mol)    # Ring detection + fusion
detect_halogens(mol)            # Element matching
detect_aliphatic_regions(mol)   # Chain detection (rarely used)
```

## Testing

Run comprehensive tests:

```bash
uv run python tests/test_functional_groups_final.py
```

Expected output:
```
======================================================================
PHARMACOPHORE FEATURE DETECTION - COMPREHENSIVE TEST
======================================================================
✓ Ibuprofen test passed (pharmacophore features only)
✓ Naproxen test passed (pharmacophore features only)
✓ Aspirin test passed (pharmacophore features only)
✓ Caffeine test passed (pharmacophore features only)
✓ No duplicates test passed (Ethanol -OH is donor only)
✓ Pyridine integration test passed (N stays in aromatic)

✓ ALL TESTS PASSED!
  - Pharmacophore features only (NO aliphatic chains)
  - No duplicate donor/acceptor features
  - Ionizable groups prioritized over H-bond features
  - Aromatic heteroatoms integrated into aromatic groups
  - Group-level granularity (benzene = 1 group)
  - Clean, pharmacologically meaningful features
======================================================================
```

## Future Directions

### Potential Uses (Not Yet Implemented)

1. **Scoring Enhancement:**
   ```python
   # Bonus if pharmacophore features align well
   if ionizable_near_ionizable_pocket:
       score_bonus -= 2.0  # kcal/mol
   ```

2. **Quality Metrics:**
   ```python
   # Report pharmacophore overlap after alignment
   overlap_score = compute_pharmacophore_overlap(ref, aligned)
   ```

3. **Visualization:**
   ```python
   # Color-code by pharmacophore type in 3D viewer
   plot_molecule_with_pharmacophore(mol, features)
   ```

4. **Documentation:**
   ```python
   # Auto-generate pharmacophore summary
   generate_pharmacophore_report(mol, features)
   ```

## References

- RDKit Feature Factory: https://www.rdkit.org/docs/RDKit_Book.html#pharmacophore-fingerprints
- Pharmacophore concept: https://en.wikipedia.org/wiki/Pharmacophore
- This implementation: Custom group-level approach for clarity

## Summary

✅ **Use for**: Molecular analysis, visualization, quality checks
❌ **Don't use for**: Alignment (use MCS instead)

**Key advantages:**
- Clean, pharmacologically meaningful features
- No noise from aliphatic chains
- No duplicate features
- Group-level granularity
- Easy to interpret

**Limitation:**
- Analysis tool only, not suitable for alignment due to lack of exact atom mapping
