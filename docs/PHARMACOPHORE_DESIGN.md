# Pharmacophore Design: Group-Level vs Atom-Level

## Current Problem: Atom-Level Features

### RDKit Default Output:
```python
Ibuprofen features:
  Donor           atoms=(14,)              # 1 atom
  Acceptor        atoms=(13,)              # 1 atom
  Aromatic        atoms=(4,5,6,7,8,9)      # 6 atoms - GOOD! Already grouped
  Hydrophobe      atoms=(0,)               # 1 atom - BAD!
  Hydrophobe      atoms=(1,)               # 1 atom - BAD!
  Hydrophobe      atoms=(2,)               # 1 atom - BAD!
  ...
  Hydrophobe      atoms=(11,)              # 12 separate hydrophobes - CHAOS!
  LumpedHydrophobe atoms=(4,5,6,7,8,9)     # 6 atoms - GOOD! Ring as one
  LumpedHydrophobe atoms=(1,0,2)           # 3 atoms - GOOD! isopropyl as one
```

### Problems:
1. **12 individual Hydrophobe** features → meaningless noise
2. **Can't distinguish** isobutyl vs benzene (both hydrophobic)
3. **No chemical context** - atom 0 alone means nothing
4. **Matching ambiguity** - which hydrophobe maps to which?

---

## Proposed Solution: Functional Group Level

### Design Philosophy:
**Pharmacophore = Chemically Meaningful Units**

Not: "This carbon is hydrophobic"
But: "This isobutyl group is hydrophobic"

Not: "6 aromatic atoms"
But: "One benzene ring"

### Functional Group Hierarchy:

```
LEVEL 1: H-Bond Groups (atomic precision needed)
  - Donor:        -OH, -NH, -NH2
  - Acceptor:     =O, -O-, =N-, -N<
  - PosIonizable: -NH3+, guanidinium
  - NegIonizable: -COO-, -SO3-, tetrazole

LEVEL 2: Aromatic Systems (ring as unit)
  - Aromatic:     benzene, naphthalene, indole, etc.
                  → Treat entire ring system as ONE feature

LEVEL 3: Aliphatic Groups (carbon chains/branches)
  - Aliphatic:    -CH3, -CH2-, isopropyl, tert-butyl
                  → Treat connected carbons as ONE feature

LEVEL 4: Special Groups (domain-specific)
  - Halogen:      F, Cl, Br, I (each separate)
  - Sulfur:       -SH, -S-, -SO2-
```

---

## Implementation Strategy

### Step 1: Detect Functional Groups (Not Atoms)

```python
class FunctionalGroup:
    """Represents a chemically meaningful group"""
    def __init__(self, type, atoms, centroid, importance):
        self.type = type              # "Donor", "Aromatic", "Aliphatic", etc.
        self.atoms = atoms            # List[int] - all atoms in this group
        self.centroid = centroid      # np.array - 3D position
        self.importance = importance  # Weight for matching

    def size(self):
        return len(self.atoms)


def detect_functional_groups(mol) -> List[FunctionalGroup]:
    """
    Detect functional groups at appropriate granularity.

    Returns list of FunctionalGroup objects.
    """
    groups = []

    # 1. H-Bond Groups (high priority, atomic precision)
    groups.extend(detect_hbond_groups(mol))

    # 2. Aromatic Systems (entire ring as one group)
    groups.extend(detect_aromatic_rings(mol))

    # 3. Aliphatic Regions (connected saturated carbons)
    groups.extend(detect_aliphatic_groups(mol))

    # 4. Special Groups (halogens, sulfur, etc.)
    groups.extend(detect_special_groups(mol))

    return groups
```

### Step 2: Aromatic Systems as Single Units

```python
def detect_aromatic_rings(mol) -> List[FunctionalGroup]:
    """
    Treat each aromatic ring system as ONE functional group.

    Benzene = 1 group (6 atoms)
    Naphthalene = 1 group (10 atoms)
    Indole = 1 group (9 atoms)
    """
    groups = []

    # Get ring info
    ring_info = mol.GetRingInfo()

    # Get aromatic rings
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(set(ring))

    # Merge fused rings (e.g., naphthalene)
    merged_systems = merge_fused_rings(aromatic_rings)

    # Create functional group for each system
    for ring_system in merged_systems:
        atoms = list(ring_system)
        centroid = compute_centroid(mol, atoms)
        groups.append(FunctionalGroup(
            type="Aromatic",
            atoms=atoms,
            centroid=centroid,
            importance=3  # High importance
        ))

    return groups


def merge_fused_rings(rings: List[Set[int]]) -> List[Set[int]]:
    """
    Merge fused aromatic rings into single systems.

    Example: Two benzenes sharing edge → naphthalene (one system)
    """
    systems = []
    for ring in rings:
        merged = False
        for system in systems:
            if ring & system:  # If they share atoms
                system.update(ring)
                merged = True
                break
        if not merged:
            systems.append(ring)
    return systems
```

### Step 3: Aliphatic Groups (Connected Carbons)

```python
def detect_aliphatic_groups(mol) -> List[FunctionalGroup]:
    """
    Group connected aliphatic carbons together.

    CH3-CH2-CH3 → One group (3 carbons)
    C(CH3)3 → One group (4 carbons, branched)
    """
    groups = []

    # Find all aliphatic carbons
    aliphatic_carbons = [
        i for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetSymbol() == 'C'
        and not mol.GetAtomWithIdx(i).GetIsAromatic()
    ]

    # Group by connectivity
    visited = set()
    for start in aliphatic_carbons:
        if start in visited:
            continue

        # BFS to find connected aliphatic region
        group_atoms = bfs_aliphatic_region(mol, start, aliphatic_carbons)
        visited.update(group_atoms)

        if len(group_atoms) >= 1:  # At least 1 carbon
            centroid = compute_centroid(mol, group_atoms)
            groups.append(FunctionalGroup(
                type="Aliphatic",
                atoms=group_atoms,
                centroid=centroid,
                importance=1  # Lower importance than H-bonds/Aromatic
            ))

    return groups


def bfs_aliphatic_region(mol, start, valid_atoms):
    """BFS to find connected aliphatic carbons"""
    visited = {start}
    queue = [start]

    while queue:
        atom_idx = queue.pop(0)
        atom = mol.GetAtomWithIdx(atom_idx)

        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in valid_atoms and n_idx not in visited:
                visited.add(n_idx)
                queue.append(n_idx)

    return list(visited)
```

### Step 4: H-Bond Groups (Atomic, but with context)

```python
def detect_hbond_groups(mol) -> List[FunctionalGroup]:
    """
    Detect H-bond donors/acceptors.

    Keep atomic precision but store context.
    """
    groups = []

    # Use RDKit feature factory
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    import os

    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    feats = factory.GetFeaturesForMol(mol)

    for feat in feats:
        if feat.GetFamily() in ['Donor', 'Acceptor', 'PosIonizable', 'NegIonizable']:
            atoms = list(feat.GetAtomIds())
            pos = feat.GetPos()
            centroid = np.array([pos.x, pos.y, pos.z])

            groups.append(FunctionalGroup(
                type=feat.GetFamily(),
                atoms=atoms,
                centroid=centroid,
                importance=5  # Highest importance!
            ))

    return groups
```

---

## Example: Ibuprofen Analysis

### Old (Atom-Level):
```
19 features:
  - 1 Donor
  - 2 Acceptors
  - 1 Aromatic (6 atoms)
  - 12 Hydrophobes (NOISE!)
  - 2 LumpedHydrophobes
```

### New (Group-Level):
```
5 functional groups:
  1. Carboxylic Acid (Donor + Acceptor combined)
     Type: NegIonizable
     Atoms: [12, 13, 14]
     Importance: 5

  2. Benzene Ring
     Type: Aromatic
     Atoms: [4, 5, 6, 7, 8, 9]
     Importance: 3

  3. Isobutyl Chain
     Type: Aliphatic
     Atoms: [0, 1, 2, 3]
     Importance: 1

  4. Methyl on benzene
     Type: Aliphatic
     Atoms: [10, 11]
     Importance: 1

  5. Carboxylic methyl
     Type: Aliphatic (or ignore - part of main group)
     Atoms: [11]
     Importance: 0.5
```

**Result: 5 meaningful groups instead of 19 ambiguous features!**

---

## Matching Strategy

### Group-to-Group Matching:

```python
def match_functional_groups(ref_groups, query_groups):
    """
    Match functional groups between ref and query.

    Scoring:
      - Type match: Required
      - Size similarity: Preferred (±2 atoms)
      - Importance: Weight matches
    """
    matches = []

    for ref_g in ref_groups:
        for query_g in query_groups:
            # Must be same type
            if ref_g.type != query_g.type:
                continue

            # Size similarity (flexible for aliphatic)
            size_diff = abs(ref_g.size() - query_g.size())
            if ref_g.type == "Aliphatic":
                max_diff = 3  # Allow flexibility
            elif ref_g.type == "Aromatic":
                max_diff = 4  # Benzene vs naphthalene ok
            else:
                max_diff = 1  # H-bonds must match closely

            if size_diff > max_diff:
                continue

            # Spatial distance
            dist = np.linalg.norm(ref_g.centroid - query_g.centroid)

            # Score
            score = ref_g.importance * (1.0 / (1.0 + dist))

            matches.append({
                'ref_group': ref_g,
                'query_group': query_g,
                'score': score,
                'distance': dist
            })

    # Sort by score
    matches.sort(key=lambda x: x['score'], reverse=True)

    return matches
```

### Convert to Atom Mapping:

```python
def groups_to_atom_mapping(group_matches):
    """
    Convert group matches to atom-level mapping.

    Strategy:
      - H-bond groups: Use all atoms (small, precise)
      - Aromatic: Use all atoms (well-defined ring)
      - Aliphatic: Use subset if size differs
    """
    atom_mapping = []

    for match in group_matches:
        ref_g = match['ref_group']
        query_g = match['query_group']

        if ref_g.type in ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable']:
            # H-bonds: Map all atoms
            for r_atom, q_atom in zip(ref_g.atoms, query_g.atoms):
                atom_mapping.append((r_atom, q_atom))

        elif ref_g.type == 'Aromatic':
            # Aromatic: Map all ring atoms
            # If sizes differ (benzene vs naphthalene), use MCS on rings
            if len(ref_g.atoms) == len(query_g.atoms):
                for r_atom, q_atom in zip(ref_g.atoms, query_g.atoms):
                    atom_mapping.append((r_atom, q_atom))
            else:
                # Run mini-MCS on just these atoms
                ring_mcs = find_ring_mcs(ref_g.atoms, query_g.atoms)
                atom_mapping.extend(ring_mcs)

        elif ref_g.type == 'Aliphatic':
            # Aliphatic: Use centroid atoms only (less critical)
            # Or skip entirely (focus on H-bonds and aromatic)
            pass

    return atom_mapping
```

---

## Comparison: Group-Level vs Atom-Level

### Ibuprofen vs Naproxen Example:

**Atom-Level (Current RDKit):**
```
Feature matches: 117 matches (NOISE!)
  - 12 hydrophobe × 11 hydrophobe = 132 combinations
  - Which carbon maps to which? AMBIGUOUS
```

**Group-Level (Proposed):**
```
Functional group matches: 4 matches (CLEAR!)

1. Carboxylic Acid (Ibu) ↔ Carboxylic Acid (Nap)
   ✓ Both NegIonizable, 3 atoms each
   Score: 5.0 (high importance)

2. Benzene (Ibu) ↔ Naphthalene (Nap)
   ~ Aromatic, 6 vs 10 atoms
   Score: 2.5 (size penalty but same type)
   Action: Run mini-MCS on rings → get 6 atom pairs

3. Isobutyl (Ibu) ↔ Methoxy (Nap)
   ~ Both aliphatic-ish, different
   Score: 0.5 (low importance)
   Action: Skip or minimal mapping

4. Methyls (Ibu) ↔ Methyls (Nap)
   ~ Low importance
   Action: Skip
```

**Result:**
- **Atom-level**: 117 ambiguous matches → hard to pick best
- **Group-level**: 4 clear matches → focus on important ones

---

## Advantages of Group-Level Design

### 1. Chemical Interpretability
```
"Map carboxylic acid to carboxylic acid"
vs
"Map atom 13 to atom 15"  ← What does this mean?
```

### 2. Reduced Ambiguity
```
1 aromatic ring vs 1 aromatic ring
vs
6 aromatic atoms vs 6 aromatic atoms in which order?
```

### 3. Importance Weighting
```python
# Automatically prioritize
importance = {
    'NegIonizable': 5,    # Critical for binding
    'Donor': 5,
    'Acceptor': 5,
    'Aromatic': 3,        # Important for shape
    'Aliphatic': 1        # Less critical
}
```

### 4. Flexible Matching
```
Benzene (6) can match Naphthalene (10)
  → Run mini-MCS on aromatic systems
  → Get 6 atom pairs from best substructure

Isobutyl (4 C) can match tert-butyl (5 C)
  → Similar hydrophobic character
  → Minor size difference ok
```

### 5. Better MCS Constraints
```python
# Instead of "use these 23 random atoms"
pharm_atoms = [0,1,2,4,5,6,7,8,9,12,13,14,...]

# Do "use these 3 important groups"
important_groups = [
    carboxylic_acid,  # 3 atoms
    benzene_ring,     # 6 atoms
    donor_oh          # 2 atoms
]
# Total: 11 atoms from meaningful groups
```

---

## Implementation Plan

### Phase 1: Group Detection
```python
# src/lig_align/molecular/functional_groups.py (new)

class FunctionalGroup:
    """Represents a chemically meaningful group"""
    pass

def detect_functional_groups(mol) -> List[FunctionalGroup]:
    """Main entry point"""
    pass

def detect_aromatic_rings(mol):
    """Aromatic systems as single units"""
    pass

def detect_aliphatic_groups(mol):
    """Connected aliphatic regions"""
    pass

def detect_hbond_groups(mol):
    """H-bond donors/acceptors"""
    pass
```

### Phase 2: Group Matching
```python
# src/lig_align/molecular/group_matching.py (new)

def match_functional_groups(ref_groups, query_groups):
    """Match groups between molecules"""
    pass

def groups_to_atom_mapping(group_matches):
    """Convert group matches to atom pairs"""
    pass
```

### Phase 3: Integration
```python
# src/lig_align/aligner.py (modify)

def step2_find_mcs(self, ref_mol, query_mol,
                   method="mcs",
                   use_functional_groups=False):
    """
    method="mcs": Standard MCS
    method="pharmacophore": Old atom-level (deprecate?)
    method="functional_groups": New group-level (recommended!)
    """
    if method == "functional_groups":
        ref_groups = detect_functional_groups(ref_mol)
        query_groups = detect_functional_groups(query_mol)
        matches = match_functional_groups(ref_groups, query_groups)
        mapping = groups_to_atom_mapping(matches)
        return mapping
```

---

## Example Usage

```python
from lig_align import run_pipeline

# Standard MCS (current)
results = run_pipeline(
    protein_pdb="protein.pdb",
    ref_ligand="ref.sdf",
    query_ligand="SMILES",
    alignment_method="mcs"
)

# Functional group based (proposed)
results = run_pipeline(
    protein_pdb="protein.pdb",
    ref_ligand="ref.sdf",
    query_ligand="SMILES",
    alignment_method="functional_groups",
    important_groups=['NegIonizable', 'Donor', 'Acceptor', 'Aromatic'],
    include_aliphatic=False  # Focus on key interactions
)
```

---

## Conclusion

### ✅ Group-Level is Superior:

1. **Chemically meaningful**: "This benzene" not "these 6 carbons"
2. **Less ambiguous**: 5 groups not 19 features
3. **Prioritizable**: Focus on important interactions
4. **Flexible**: Can handle size differences
5. **Interpretable**: Easy to explain why molecules match

### Implementation Complexity:
- **Medium** (~300 lines)
- Group detection: 150 lines
- Group matching: 100 lines
- Integration: 50 lines

### Expected Improvements:
- Better alignment for diverse molecules
- More robust to size differences
- Focus on pharmacologically relevant features
- Easier to explain and debug

**Recommendation: Implement group-level functional groups instead of atom-level pharmacophore features.**
