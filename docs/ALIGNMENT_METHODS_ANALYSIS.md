# Alignment Methods: Deep Analysis

## Current: MCS-Based Alignment

### How it works:
```
1. Find MCS → Get atom pairs [(ref_idx, query_idx), ...]
2. Coordinate Surgery → Teleport MCS atoms to reference positions
3. MMFF Relax → Relax non-MCS atoms (MCS fixed)
4. Score → Vina scoring
5. Optimize → Gradient-based torsion optimization (MCS frozen as anchor)
```

### Why it works:
- ✅ **Anchor atoms**: MCS atoms act as fixed anchors
- ✅ **Torsion optimization**: Rotate around rotatable bonds
- ✅ **Differentiable**: Can backpropagate through torsions
- ✅ **Kinematic tree**: MCS = root, other atoms = children
- ✅ **Forward kinematics**: Rodrigues rotation formula

### Key advantage:
**MCS atoms = rigid scaffold**, rest of molecule can optimize while staying aligned.

---

## Proposed 1: Pharmacophore-Constrained MCS

### How it works:
```
1. Detect pharmacophore features (Donor, Acceptor, Aromatic)
2. Get atoms for these features
3. Run MCS ONLY on pharmacophore atoms
4. Rest is SAME as current MCS pipeline!
```

### Why it works:
- ✅ Still returns atom-level mapping [(ref_idx, query_idx), ...]
- ✅ These atoms become anchors
- ✅ Existing optimization works unchanged
- ✅ Chemically more meaningful (focus on key interactions)

### Implementation:
```python
# Easy - just constrain MCS search
pharm_atoms_ref = get_pharmacophore_atoms(ref_mol, ['Donor', 'Acceptor', 'Aromatic'])
pharm_atoms_query = get_pharmacophore_atoms(query_mol, ['Donor', 'Acceptor', 'Aromatic'])

# Create submolecules with only these atoms
mcs = find_mcs(ref_submol, query_submol)

# Returns: [(ref_idx, query_idx), ...] ← Same as regular MCS!
# Rest of pipeline unchanged!
```

### Verdict: ✅ **FEASIBLE AND PRACTICAL**

---

## Proposed 2: Shape-Based Alignment

### How it would work (naive):
```
1. Overlay molecules by shape (whole molecule alignment)
2. Get transformation matrix (rotation R, translation t)
3. Apply transform to query: coords' = R @ coords + t
4. Score
5. Optimize ??? ← PROBLEM!
```

### ❌ **CRITICAL PROBLEM: No Anchors!**

#### Problem 1: No Fixed Reference
```python
# Current MCS optimization:
def optimize_torsions():
    # MCS atoms = FIXED (anchors)
    # Rotate around bonds
    # MCS stays aligned to reference
    ✓ Works!

# Shape-based - no anchors:
def optimize_torsions():
    # What atoms are fixed?
    # If nothing is fixed, molecule drifts
    # Loses alignment!
    ✗ Doesn't work!
```

#### Problem 2: Whole-Molecule Transform
```
Shape alignment gives: T = [R | t]  (4×4 transform matrix)

This is a GLOBAL transform - applies to entire molecule.

But torsion optimization is LOCAL:
  - Rotate around individual bonds
  - Each rotation changes LOCAL geometry
  - No concept of global transform

These are INCOMPATIBLE!
```

#### Problem 3: Optimization Objective
```python
# Current optimization:
minimize: Vina_score(coords)
subject to: MCS atoms at reference positions  ← Constraint!

# Shape-based - no constraint:
minimize: Vina_score(coords)
subject to: ???

Without constraint, molecule will:
  1. Rotate to random orientation (loses alignment)
  2. Translate away from reference position
  3. Completely lose shape overlay
```

### Could we fix this?

#### Attempt 1: Add alignment loss to objective
```python
# Optimize both score AND alignment
loss = vina_score + λ * alignment_loss

# Where alignment_loss = shape_difference(query, reference)
```

**Problem:**
- Shape comparison is expensive (requires computing moments, overlays)
- Not differentiable through torsions easily
- Optimization becomes much slower
- Needs careful λ tuning

#### Attempt 2: Pick "pseudo-anchors" from shape
```python
# After shape overlay, pick atoms that are close to reference
pseudo_anchors = find_closest_atoms(query, reference, top_k=10)

# Freeze these during optimization
```

**Problem:**
- Arbitrary! No chemical meaning
- Atoms might not even be equivalent
- Essentially reinventing MCS with worse chemistry

#### Attempt 3: No optimization - just scoring
```python
# Shape align → Score → Done (no optimization)
```

**Problem:**
- Much worse results (no refinement)
- Defeats the purpose of gradient-based optimization
- Why use PyTorch at all?

---

## Shape-Based: When Would It Be Useful?

### Use Case 1: Virtual Screening (Pre-filter)
```python
# Screen 1M molecules
for mol in database:
    shape_similarity = compute_shape_similarity(mol, reference)
    if shape_similarity > threshold:
        candidates.append(mol)

# Then use MCS pipeline on candidates
for mol in candidates:
    results = run_pipeline(mol, method="mcs")  # With optimization!
```

**This works!** But shape is just a filter, not the alignment method.

### Use Case 2: Very Dissimilar Molecules (No MCS)
```
Reference: Small rigid molecule
Query:     Large flexible molecule, minimal overlap

Problem: MCS too small or nonexistent
Solution: Shape overlay as last resort
```

**But:** Without optimization, results will be poor.

---

## Comparison Table

| Method | Anchor Atoms | Optimization | Differentiable | Use Case |
|--------|-------------|--------------|----------------|----------|
| **MCS** | ✅ Yes (MCS atoms) | ✅ Yes (torsions) | ✅ Yes | Default, best for similar molecules |
| **Pharmacophore-MCS** | ✅ Yes (pharm atoms) | ✅ Yes (torsions) | ✅ Yes | Focus on key interactions |
| **Shape-Based** | ❌ No | ❌ Very difficult | ⚠️ Complex | Screening only, no optimization |

---

## Deep Dive: Why Torsion Optimization Needs Anchors

### Kinematic Tree Structure

```
Current MCS-based approach:

         [MCS atoms] ← FIXED (root)
              |
         Rotatable Bond 1
              |
         +---------+
         |         |
      Atom A    Atom B
         |         |
    Bond 2    Bond 3
         |         |
        ...       ...

- MCS = root of tree
- Each rotatable bond = edge
- Rotate bonds → atoms move
- But MCS stays fixed → maintains alignment
```

### Without Anchors (Shape-based):

```
         [All atoms floating]
              |
         Rotatable Bond 1
              |
         +---------+
         |         |
      Atom A    Atom B

Problem: Where is the root?
- No fixed reference
- Rotating Bond 1 → entire left subtree rotates
- But what's it rotating relative to?
- No anchor → coordinate system drifts
```

### Mathematical View

**MCS-based (constrained optimization):**
```
minimize: E(θ)                    # Energy as function of torsions
subject to: x_mcs = x_ref         # MCS positions fixed

Forward kinematics:
  x_i = FK(θ, x_mcs)              # Compute positions from torsions + anchors

Gradient:
  ∂E/∂θ = ∂E/∂x · ∂x/∂θ          # Chain rule, x_mcs constant
```

**Shape-based (unconstrained):**
```
minimize: E(θ)                    # Energy as function of torsions
subject to: ???                   # No constraints!

Problem:
  x_i = FK(θ, ???)                # No anchor positions!

Even if we optimize, molecule can:
  - Translate anywhere
  - Rotate arbitrarily
  - Lose alignment completely
```

---

## Hybrid Approach: Shape-Aligned + Soft Anchors?

### Could we do this?

```python
# 1. Shape overlay → get initial alignment
transform = shape_align(query, reference)
query_coords = apply_transform(query_coords, transform)

# 2. Pick "soft anchors" = atoms closest to reference surface
soft_anchors = find_closest_atoms(query_coords, ref_coords, k=5)

# 3. Optimize with soft constraint
loss = vina_score + λ * distance_to_soft_anchors(soft_anchors)
```

**Analysis:**
- ⚠️ **Soft anchors are arbitrary** - no chemical meaning
- ⚠️ **λ tuning is hard** - too small = drift, too large = can't optimize
- ⚠️ **Distance loss is not differentiable through torsions easily**
- ⚠️ **Computationally expensive**
- ⚠️ **Results likely worse than MCS anyway**

**Verdict:** Not worth the complexity.

---

## Real-World Analogy

### MCS-based = Building with Foundation
```
You have a foundation (MCS atoms) bolted to the ground.
You can rotate and adjust the walls (non-MCS atoms).
But the foundation stays fixed.
✓ Structure remains aligned while you optimize.
```

### Shape-based = Building on Water
```
You have a floating platform (shape overlay).
You try to adjust the walls.
But the entire structure drifts and rotates.
✗ You lose your reference frame.
```

---

## Pharmacophore Feature Analysis

### Which features make good anchors?

| Feature | Atoms per Feature | Good Anchor? | Why? |
|---------|------------------|--------------|------|
| **Donor** | 1 | ✅ Excellent | Single atom, clear chemistry |
| **Acceptor** | 1 | ✅ Excellent | Single atom, clear chemistry |
| **NegIonizable** | 2-3 | ✅ Good | Small group, important interaction |
| **PosIonizable** | 1-2 | ✅ Good | Small group, important interaction |
| **Aromatic** | 6-10 | ⚠️ Medium | Multiple atoms, but well-defined |
| **Hydrophobe** | 1-12 | ❌ Poor | Too many, too ambiguous |
| **LumpedHydrophobe** | 2-3 | ⚠️ Medium | Better than individual, still vague |

### Recommended Feature Set for Pharmacophore-MCS:
```python
key_features = [
    "Donor",          # H-bond donors (OH, NH)
    "Acceptor",       # H-bond acceptors (O, N)
    "Aromatic",       # Aromatic rings
    "NegIonizable",   # Carboxylates, etc.
    "PosIonizable"    # Amines, etc.
]

# Skip Hydrophobe - too many, too ambiguous
```

---

## Implementation Recommendation

### Phase 1: Pharmacophore-Constrained MCS ✅

**Implementation Plan:**

```python
# src/lig_align/molecular/pharmacophore.py (new file)

def get_pharmacophore_atoms(mol, features=['Donor', 'Acceptor', 'Aromatic']):
    """
    Extract atom IDs for specified pharmacophore features.

    Returns: set of atom indices
    """
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    import os

    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    mol_feats = factory.GetFeaturesForMol(mol)

    atom_set = set()
    for feat in mol_feats:
        if feat.GetFamily() in features:
            atom_set.update(feat.GetAtomIds())

    return atom_set


def find_mcs_pharmacophore_constrained(ref_mol, query_mol, features):
    """
    Find MCS constrained to pharmacophore atoms.

    Returns: [(ref_idx, query_idx), ...] - same as regular MCS!
    """
    # 1. Get pharmacophore atoms
    ref_pharm_atoms = get_pharmacophore_atoms(ref_mol, features)
    query_pharm_atoms = get_pharmacophore_atoms(query_mol, features)

    # 2. Create atom map to original indices
    ref_idx_map = {new_idx: orig_idx for new_idx, orig_idx in enumerate(ref_pharm_atoms)}
    query_idx_map = {new_idx: orig_idx for new_idx, orig_idx in enumerate(query_pharm_atoms)}

    # 3. Create submolecules with only these atoms
    ref_sub = create_submol(ref_mol, ref_pharm_atoms)
    query_sub = create_submol(query_mol, query_pharm_atoms)

    # 4. Find MCS on submolecules
    mcs_sub = find_mcs(ref_sub, query_sub)

    # 5. Map back to original indices
    mcs_orig = [(ref_idx_map[r], query_idx_map[q]) for r, q in mcs_sub]

    return mcs_orig


# src/lig_align/aligner.py (modify)

def step2_find_mcs(self, ref_mol, query_mol,
                   mcs_mode="single",
                   use_pharmacophore=False,
                   pharm_features=None):
    """
    Find MCS with optional pharmacophore constraint.
    """
    if use_pharmacophore:
        from .molecular.pharmacophore import find_mcs_pharmacophore_constrained

        if pharm_features is None:
            pharm_features = ['Donor', 'Acceptor', 'Aromatic']

        mapping = find_mcs_pharmacophore_constrained(ref_mol, query_mol, pharm_features)
        return mapping
    else:
        # Existing MCS logic
        ...
```

**Usage:**
```python
# Regular MCS
mapping = aligner.step2_find_mcs(ref_mol, query_mol)

# Pharmacophore-constrained MCS
mapping = aligner.step2_find_mcs(
    ref_mol,
    query_mol,
    use_pharmacophore=True,
    pharm_features=['Donor', 'Acceptor', 'Aromatic']
)

# Rest of pipeline UNCHANGED!
aligned = aligner.step3_batched_kabsch_alignment(...)
scores = aligner.step4_vina_scoring(...)
optimized = aligner.step6_refine_pose(...)  # Works exactly the same!
```

---

## Final Verdict

### ✅ IMPLEMENT: Pharmacophore-Constrained MCS
- Easy to implement (~100 lines of code)
- Compatible with existing pipeline
- Provides chemical insight
- Enables optimization
- **Best of both worlds: MCS structure + Pharmacophore chemistry**

### ❌ DON'T IMPLEMENT: Shape-Based Alignment
- Cannot support gradient optimization effectively
- No anchor atoms = drifting problem
- Would need entirely different pipeline
- Use case is limited (screening filter only)
- Not worth the complexity

### 🎯 Recommendation
Focus on Pharmacophore-Constrained MCS. It gives you:
1. Chemical relevance (focus on key interactions)
2. Optimization capability (anchors for torsions)
3. Easy implementation (small addition to existing code)
4. Better results for diverse molecules (when MCS is small)

Shape-based methods are better suited for separate tools (screening, filtering), not gradient-based pose optimization.
