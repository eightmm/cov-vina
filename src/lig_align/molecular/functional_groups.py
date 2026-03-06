"""
Pharmacophore Feature Detection for Molecular Analysis and Visualization.

This module detects pharmacophore features at the functional group level
(not individual atoms) for molecular analysis, visualization, and quality checks.

**Purpose**: Analysis and visualization tool, NOT used for alignment.
Alignment uses MCS (Maximum Common Substructure) only.

Standard Pharmacophore Features:
- H-Bond Donors: -OH, -NH, -NH2, -SH (importance: 5)
- H-Bond Acceptors: =O, -O-, =N-, -N< (importance: 5)
- Positive Ionizable: Amines, Guanidines, Imidazoles (importance: 5)
- Negative Ionizable: Carboxylates, Sulfonates, Phosphates (importance: 5)
- Aromatic: Benzene, Naphthalene, Pyridine, etc. (importance: 3)
- Halogen: F, Cl, Br, I (importance: 2)

Note: Aliphatic chains are NOT pharmacophore features and are excluded.
"""

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, Descriptors
from rdkit import RDConfig
import numpy as np
import os
from typing import List, Set, Tuple
from dataclasses import dataclass


@dataclass
class FunctionalGroup:
    """
    Represents a chemically meaningful functional group.

    Attributes:
        type: Group type (e.g., "Donor", "Aromatic", "Aliphatic")
        subtype: More specific type (e.g., "Carboxylic", "Benzene", "Alkyl")
        atoms: List of atom indices in this group
        centroid: 3D center position
        importance: Weight for matching (1-5, higher = more important)
    """
    type: str
    subtype: str
    atoms: List[int]
    centroid: np.ndarray
    importance: int

    def size(self) -> int:
        return len(self.atoms)

    def __repr__(self):
        return f"{self.subtype}({self.type}, {len(self.atoms)} atoms, importance={self.importance})"


def detect_pharmacophore_features(mol: Chem.Mol) -> List[FunctionalGroup]:
    """
    Detect pharmacophore features in a molecule.

    Pharmacophore features are functional groups critical for binding:
    - Ionizable groups (importance 5)
    - H-bond donors/acceptors (importance 5)
    - Aromatic systems (importance 3)
    - Halogens (importance 2)

    Aliphatic chains are NOT included (not pharmacophore features).

    Args:
        mol: RDKit molecule (must have 3D coordinates)

    Returns:
        List of FunctionalGroup objects, sorted by importance (descending)

    Example:
        >>> mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        >>> mol = Chem.AddHs(mol)
        >>> AllChem.EmbedMolecule(mol)
        >>> features = detect_pharmacophore_features(mol)
        >>> print(f"Found {len(features)} pharmacophore features")
        Found 2 pharmacophore features
        >>> for f in features:
        ...     print(f"{f.type}: {f.subtype}")
        NegIonizable: Carboxylate
        Aromatic: Benzene
    """
    groups = []

    # Priority 1: H-Bond Donors/Acceptors and Ionizable (importance 5)
    hbond_ionizable_groups = detect_hbond_ionizable(mol)

    # Track atoms in ionizable groups to avoid duplicates
    ionizable_atoms = set()
    for g in hbond_ionizable_groups:
        if 'Ionizable' in g.type:
            ionizable_atoms.update(g.atoms)

    # Filter out donor/acceptor features that overlap with ionizable groups
    # Also merge donor/acceptor that share atoms (e.g., -OH is both)
    donor_acceptor_atoms = set()
    for g in hbond_ionizable_groups:
        if g.type in ['Donor', 'Acceptor']:
            # Skip if any atoms overlap with ionizable groups
            if any(atom in ionizable_atoms for atom in g.atoms):
                continue
            # Skip acceptor if already have donor for same atoms (donor takes priority)
            if g.type == 'Acceptor' and any(atom in donor_acceptor_atoms for atom in g.atoms):
                continue
            # Add donor atoms to track
            if g.type == 'Donor':
                donor_acceptor_atoms.update(g.atoms)
        groups.append(g)

    # Priority 2: Aromatic Systems (importance 3)
    aromatic_groups = detect_aromatic_systems(mol)

    # Track aromatic atoms
    aromatic_atoms = set()
    for g in aromatic_groups:
        aromatic_atoms.update(g.atoms)

    # Filter out donor/acceptor in aromatic rings (like pyridine N)
    groups_filtered = []
    for g in groups:
        if g.type in ['Donor', 'Acceptor']:
            # Skip if in aromatic ring
            if any(atom in aromatic_atoms for atom in g.atoms):
                continue
        groups_filtered.append(g)

    groups = groups_filtered
    groups.extend(aromatic_groups)

    # Priority 3: Halogen Bonds (importance 2)
    groups.extend(detect_halogens(mol))

    # NOTE: Aliphatic regions are NOT pharmacophore features - excluded

    # Sort by importance (descending)
    groups.sort(key=lambda g: g.importance, reverse=True)

    return groups


# Alias for backward compatibility
def detect_functional_groups(mol: Chem.Mol,
                            include_aliphatic: bool = False,
                            min_aliphatic_size: int = 2) -> List[FunctionalGroup]:
    """
    Detect functional groups in a molecule.

    DEPRECATED: Use detect_pharmacophore_features() instead.
    This function is kept for backward compatibility only.

    Args:
        mol: RDKit molecule (must have 3D coordinates)
        include_aliphatic: Whether to detect aliphatic groups (default: False)
        min_aliphatic_size: Minimum carbon count for aliphatic groups

    Returns:
        List of FunctionalGroup objects, sorted by importance (descending)
    """
    if include_aliphatic:
        # Full implementation with aliphatic
        groups = list(detect_pharmacophore_features(mol))
        groups.extend(detect_aliphatic_regions(mol, min_size=min_aliphatic_size))
        groups.sort(key=lambda g: g.importance, reverse=True)
        return groups
    else:
        return detect_pharmacophore_features(mol)


def detect_hbond_ionizable(mol: Chem.Mol) -> List[FunctionalGroup]:
    """
    Detect H-bond donors, acceptors, and ionizable groups.

    These are the most important for binding - atomic precision maintained.

    Features detected:
    - Donor: -OH, -NH, -NH2, -SH
    - Acceptor: =O, -O-, =N-, -N<
    - NegIonizable: -COO-, -SO3-, -PO3-
    - PosIonizable: -NH3+, guanidinium, imidazole
    """
    groups = []

    # Load RDKit feature factory
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    feats = factory.GetFeaturesForMol(mol)

    for feat in feats:
        family = feat.GetFamily()

        # Only process H-bond and ionizable features
        if family not in ['Donor', 'Acceptor', 'NegIonizable', 'PosIonizable']:
            continue

        atoms = list(feat.GetAtomIds())
        pos = feat.GetPos()
        centroid = np.array([pos.x, pos.y, pos.z])

        # Determine subtype based on atoms
        subtype = _classify_hbond_group(mol, atoms, family)

        groups.append(FunctionalGroup(
            type=family,
            subtype=subtype,
            atoms=atoms,
            centroid=centroid,
            importance=5  # Highest priority
        ))

    return groups


def _classify_hbond_group(mol: Chem.Mol, atoms: List[int], family: str) -> str:
    """Classify H-bond/ionizable groups into subtypes"""

    if family == "NegIonizable":
        # Check for carboxylate, sulfonate, phosphate
        for atom_idx in atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                return "Carboxylate"
            elif atom.GetSymbol() == 'S':
                return "Sulfonate"
            elif atom.GetSymbol() == 'P':
                return "Phosphate"
        return "NegIonizable"

    elif family == "PosIonizable":
        # Check for amine, guanidine, imidazole
        smarts_patterns = {
            "Guanidine": Chem.MolFromSmarts("NC(=N)N"),
            "Imidazole": Chem.MolFromSmarts("c1ncnc1"),
            "Amine": Chem.MolFromSmarts("[N;+0,+1]")
        }
        for name, pattern in smarts_patterns.items():
            if pattern and mol.HasSubstructMatch(pattern):
                match = mol.GetSubstructMatch(pattern)
                if any(a in match for a in atoms):
                    return name
        return "PosIonizable"

    elif family == "Donor":
        for atom_idx in atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                return "Hydroxyl"
            elif atom.GetSymbol() == 'N':
                return "Amine"
            elif atom.GetSymbol() == 'S':
                return "Thiol"
        return "Donor"

    elif family == "Acceptor":
        for atom_idx in atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                # Check if carbonyl or ether
                if any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds()):
                    return "Carbonyl"
                else:
                    return "Ether"
            elif atom.GetSymbol() == 'N':
                return "Nitrogen"
        return "Acceptor"

    return family


def detect_aromatic_systems(mol: Chem.Mol) -> List[FunctionalGroup]:
    """
    Detect aromatic ring systems as single functional groups.

    Each aromatic ring system (including fused rings like naphthalene) is
    treated as ONE group, not individual atoms.

    Examples:
    - Benzene: 1 group (6 atoms)
    - Naphthalene: 1 group (10 atoms, fused)
    - Biphenyl: 2 groups (6 atoms each, not fused)
    """
    groups = []

    ring_info = mol.GetRingInfo()

    # Find all aromatic rings
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(set(ring))

    if not aromatic_rings:
        return groups

    # Merge fused aromatic systems
    merged_systems = _merge_fused_rings(aromatic_rings)

    # Create functional group for each system
    for ring_system in merged_systems:
        atoms = sorted(list(ring_system))
        centroid = _compute_centroid(mol, atoms)

        # Classify aromatic system
        subtype = _classify_aromatic_system(mol, atoms)

        groups.append(FunctionalGroup(
            type="Aromatic",
            subtype=subtype,
            atoms=atoms,
            centroid=centroid,
            importance=3
        ))

    return groups


def _merge_fused_rings(rings: List[Set[int]]) -> List[Set[int]]:
    """Merge fused aromatic rings into single systems"""
    systems = []

    for ring in rings:
        merged = False
        for system in systems:
            if ring & system:  # If they share atoms (fused)
                system.update(ring)
                merged = True
                break

        if not merged:
            systems.append(ring.copy())

    # Continue merging until no more changes
    changed = True
    while changed:
        changed = False
        new_systems = []
        used = set()

        for i, sys1 in enumerate(systems):
            if i in used:
                continue

            for j, sys2 in enumerate(systems[i+1:], start=i+1):
                if j in used:
                    continue

                if sys1 & sys2:
                    sys1.update(sys2)
                    used.add(j)
                    changed = True

            new_systems.append(sys1)

        systems = new_systems

    return systems


def _classify_aromatic_system(mol: Chem.Mol, atoms: List[int]) -> str:
    """Classify aromatic system type"""
    size = len(atoms)

    # Check heteroatoms
    has_nitrogen = any(mol.GetAtomWithIdx(i).GetSymbol() == 'N' for i in atoms)
    has_oxygen = any(mol.GetAtomWithIdx(i).GetSymbol() == 'O' for i in atoms)
    has_sulfur = any(mol.GetAtomWithIdx(i).GetSymbol() == 'S' for i in atoms)

    if size == 5:
        if has_nitrogen:
            return "Pyrrole/Imidazole"
        elif has_oxygen:
            return "Furan"
        elif has_sulfur:
            return "Thiophene"
        return "5-Ring"

    elif size == 6:
        if has_nitrogen:
            return "Pyridine"
        return "Benzene"

    elif size == 9:
        if has_nitrogen:
            return "Indole"
        return "Indene"

    elif size == 10:
        return "Naphthalene"

    elif size == 14:
        return "Anthracene"

    else:
        return f"Aromatic-{size}"


def detect_halogens(mol: Chem.Mol) -> List[FunctionalGroup]:
    """
    Detect halogen atoms (F, Cl, Br, I).

    Halogens can form halogen bonds and affect binding.
    Each halogen is a separate group (single atom).
    """
    groups = []

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
            idx = atom.GetIdx()
            pos = mol.GetConformer().GetAtomPosition(idx)
            centroid = np.array([pos.x, pos.y, pos.z])

            groups.append(FunctionalGroup(
                type="Halogen",
                subtype=atom.GetSymbol(),
                atoms=[idx],
                centroid=centroid,
                importance=2
            ))

    return groups


def detect_aliphatic_regions(mol: Chem.Mol, min_size: int = 2) -> List[FunctionalGroup]:
    """
    Detect aliphatic (non-aromatic, hydrophobic) regions.

    Connected aliphatic carbons are grouped together as single units.

    Examples:
    - CH3-CH2-CH3: 1 group (3 carbons)
    - C(CH3)3: 1 group (4 carbons)
    """
    groups = []

    # Find all aliphatic carbons
    aliphatic_carbons = []
    for atom in mol.GetAtoms():
        if (atom.GetSymbol() == 'C' and
            not atom.GetIsAromatic() and
            atom.GetAtomicNum() == 6):
            aliphatic_carbons.append(atom.GetIdx())

    if not aliphatic_carbons:
        return groups

    # Group by connectivity
    visited = set()

    for start in aliphatic_carbons:
        if start in visited:
            continue

        # BFS to find connected aliphatic region
        region = _bfs_aliphatic_region(mol, start, aliphatic_carbons)
        visited.update(region)

        if len(region) >= min_size:
            atoms = sorted(list(region))
            centroid = _compute_centroid(mol, atoms)

            # Classify aliphatic type
            subtype = _classify_aliphatic(mol, atoms)

            groups.append(FunctionalGroup(
                type="Aliphatic",
                subtype=subtype,
                atoms=atoms,
                centroid=centroid,
                importance=1
            ))

    return groups


def _bfs_aliphatic_region(mol: Chem.Mol, start: int, valid_atoms: List[int]) -> Set[int]:
    """BFS to find connected aliphatic carbons"""
    valid_set = set(valid_atoms)
    visited = {start}
    queue = [start]

    while queue:
        atom_idx = queue.pop(0)
        atom = mol.GetAtomWithIdx(atom_idx)

        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in valid_set and n_idx not in visited:
                visited.add(n_idx)
                queue.append(n_idx)

    return visited


def _classify_aliphatic(mol: Chem.Mol, atoms: List[int]) -> str:
    """Classify aliphatic group type"""
    size = len(atoms)

    # Check branching
    max_degree = max(mol.GetAtomWithIdx(i).GetDegree() for i in atoms)

    if size == 1:
        return "Methyl"
    elif size == 2:
        return "Ethyl"
    elif size == 3:
        if max_degree >= 3:
            return "Isopropyl"
        return "Propyl"
    elif size == 4:
        if max_degree >= 4:
            return "Tert-butyl"
        elif max_degree >= 3:
            return "Isobutyl"
        return "Butyl"
    else:
        if max_degree >= 3:
            return f"Branched-C{size}"
        return f"Alkyl-C{size}"


def _compute_centroid(mol: Chem.Mol, atoms: List[int]) -> np.ndarray:
    """Compute geometric center of atoms"""
    conf = mol.GetConformer()
    coords = []
    for i in atoms:
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    return np.array(coords).mean(axis=0)


def print_functional_groups(mol: Chem.Mol, groups: List[FunctionalGroup]):
    """Pretty print functional groups for debugging"""
    print(f"\n{'='*70}")
    print(f"Functional Groups Analysis: {mol.GetProp('_Name') if mol.HasProp('_Name') else 'Molecule'}")
    print(f"{'='*70}")
    print(f"Total: {len(groups)} groups\n")

    # Group by type
    from collections import defaultdict
    by_type = defaultdict(list)
    for g in groups:
        by_type[g.type].append(g)

    for type_name, type_groups in sorted(by_type.items(), key=lambda x: -x[1][0].importance):
        print(f"\n{type_name} ({len(type_groups)} groups):")
        for g in type_groups:
            atom_str = f"atoms={g.atoms}" if len(g.atoms) <= 6 else f"atoms=[...{len(g.atoms)} atoms...]"
            print(f"  {g.subtype:20s} importance={g.importance}, {atom_str}")

    print(f"\n{'='*70}\n")


if __name__ == "__main__":
    # Test with example molecules
    from rdkit.Chem import AllChem

    test_molecules = {
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    }

    for name, smiles in test_molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp('_Name', name)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        groups = detect_functional_groups(mol)
        print_functional_groups(mol, groups)
