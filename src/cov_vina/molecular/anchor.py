"""
Covalent warhead detection and anchor point generation.

Replaces MCS-based alignment for covalent docking:
- Detects reactive warhead groups on ligands via SMARTS patterns
- Extracts reactive residue atom coordinates from protein PDB
- Creates coordinate constraints (coordMap) for conformer generation
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D


# ---------------------------------------------------------------------------
# Warhead SMARTS registry
# ---------------------------------------------------------------------------
# Each entry: (smarts, reactive_atom_map_num, description)
# The atom tagged with map number `reactive_atom_map_num` is the one that
# forms the new covalent bond with the protein nucleophile.

_WARHEAD_REGISTRY: list[tuple[str, int, str]] = [
    # --- Michael acceptors (1,4-conjugate addition) ---
    ("[CH2:1]=[CH]C(=O)[N,n]",        1, "acrylamide"),
    ("[CH2:1]=[CH]C(=O)[OH]",         1, "acrylic_acid"),
    ("[CH2:1]=[CH]C(=O)O[#6]",        1, "acrylate"),
    ("[CH2:1]=[CH]C(=O)[#6]",         1, "enone"),
    ("[CH2:1]=[CH]S(=O)(=O)[N,n]",    1, "vinyl_sulfonamide"),
    ("[CH2:1]=[CH]S(=O)(=O)[#6]",     1, "vinyl_sulfone"),
    ("O=C1[CH:1]=[CH]C(=O)[NH,N]1",   1, "maleimide"),

    # --- alpha-Halo carbonyl (SN2 displacement) ---
    ("Cl[CH2:1]C(=O)[N,n]",           1, "chloroacetamide"),
    ("Br[CH2:1]C(=O)[N,n]",           1, "bromoacetamide"),
    ("I[CH2:1]C(=O)[N,n]",            1, "iodoacetamide"),
    ("F[CH2:1]C(=O)[N,n]",            1, "fluoroacetamide"),
    ("Cl[C:1](F)C(=O)[N,n]",          1, "chlorofluoroacetamide"),

    # --- Strained ring opening (SN2) ---
    ("[CH:1]1OC1",                     1, "epoxide"),
    ("[CH:1]1NC1",                     1, "aziridine"),
    ("[CH:1]1SC1",                     1, "thiirane"),

    # --- Triple bond electrophiles ---
    ("N#[C:1]c",                       1, "aryl_nitrile"),        # Aromatic nitrile (more specific)
    ("N#[C:1]C([#6])",                 1, "alkyl_nitrile"),       # Aliphatic nitrile
    ("[CH:1]#CC(=O)[N,n]",            1, "propiolamide"),
    ("[C:1]#CC(=O)[N,n]",             1, "propargylamide"),

    # --- Reversible Michael (cyanoacrylamide) ---
    ("N#CC=[C:1]C(=O)[N,n]",          1, "cyanoacrylamide"),

    # --- Disulfide exchange ---
    ("[S:1]S[#6]",                     1, "disulfide"),

    # --- Sulfonyl fluoride (SuFEx) ---
    ("F[S:1](=O)(=O)[c,C]",           1, "sulfonyl_fluoride"),

    # --- Aldehyde (hemithioacetal) ---
    ("[CH1:1]=O",                      1, "aldehyde"),           # Aldehyde carbon
    ("O=[C:1]C(=O)[N,n]",              1, "alpha_ketoamide"),    # α-ketoamide (glyoxal amide)

    # --- Isothiocyanate ---
    ("[C:1](=N)=S",                    1, "isothiocyanate"),     # Isothiocyanate carbon

    # --- Activated esters ---
    ("[C:1](=O)On1c(=O)cccc1",        1, "nhs_ester"),          # N-hydroxysuccinimide
    ("[C:1](=O)OC(F)(F)F",            1, "tfe_ester"),          # Trifluoroethyl ester

    # --- Acyl fluoride ---
    ("[C:1](=O)F",                     1, "acyl_fluoride"),

    # --- Boronic acid (for SER) ---
    ("[B:1]([OH])[OH]",                1, "boronic_acid"),

    # --- Phosphonate (for SER) ---
    ("[P:1](=O)([OH])[OH]",           1, "phosphonate"),
]


# ---------------------------------------------------------------------------
# Warhead-Residue compatibility matrix
# ---------------------------------------------------------------------------
# Maps (warhead_type, residue_name) → compatibility level
# "GOOD": Well-established, fast reaction
# "SLOW": Possible but slow/rare
# "NO": Not compatible
#
# If not in this dict, defaults to "SLOW" (cautious approach)

WARHEAD_RESIDUE_COMPATIBILITY = {
    # Michael acceptors - CYS only (strong nucleophile needed)
    ("acrylamide", "CYS"): "GOOD",
    ("acrylic_acid", "CYS"): "GOOD",
    ("acrylate", "CYS"): "GOOD",
    ("enone", "CYS"): "GOOD",
    ("vinyl_sulfonamide", "CYS"): "GOOD",
    ("vinyl_sulfone", "CYS"): "GOOD",
    ("maleimide", "CYS"): "GOOD",
    ("cyanoacrylamide", "CYS"): "GOOD",

    # Michael acceptors - other residues (not compatible)
    ("acrylamide", "SER"): "NO",
    ("acrylamide", "THR"): "NO",
    ("acrylamide", "LYS"): "NO",
    ("vinyl_sulfonamide", "SER"): "NO",

    # Alpha-halo carbonyls - CYS best, others possible but slow
    ("chloroacetamide", "CYS"): "GOOD",
    ("bromoacetamide", "CYS"): "GOOD",
    ("iodoacetamide", "CYS"): "GOOD",
    ("chloroacetamide", "SER"): "SLOW",
    ("chloroacetamide", "HIS"): "SLOW",
    ("chloroacetamide", "LYS"): "SLOW",

    # Epoxides - multiple residues possible
    ("epoxide", "CYS"): "GOOD",
    ("epoxide", "LYS"): "GOOD",
    ("epoxide", "HIS"): "GOOD",
    ("epoxide", "SER"): "SLOW",
    ("epoxide", "THR"): "SLOW",

    # Nitriles - CYS and LYS (reversible)
    ("aryl_nitrile", "CYS"): "GOOD",
    ("aryl_nitrile", "LYS"): "GOOD",
    ("alkyl_nitrile", "CYS"): "GOOD",
    ("alkyl_nitrile", "LYS"): "GOOD",
    ("aryl_nitrile", "SER"): "SLOW",

    # Boronic acid - SER/THR/TYR only (serine proteases)
    ("boronic_acid", "SER"): "GOOD",
    ("boronic_acid", "THR"): "GOOD",
    ("boronic_acid", "TYR"): "GOOD",
    ("boronic_acid", "CYS"): "NO",
    ("boronic_acid", "LYS"): "NO",

    # Phosphonate - SER/THR only
    ("phosphonate", "SER"): "GOOD",
    ("phosphonate", "THR"): "GOOD",
    ("phosphonate", "CYS"): "NO",

    # Sulfonyl fluoride - very reactive, multiple residues
    ("sulfonyl_fluoride", "CYS"): "GOOD",
    ("sulfonyl_fluoride", "SER"): "GOOD",
    ("sulfonyl_fluoride", "TYR"): "GOOD",
    ("sulfonyl_fluoride", "LYS"): "GOOD",
    ("sulfonyl_fluoride", "HIS"): "GOOD",

    # Aldehyde - reversible, CYS/SER/LYS
    ("aldehyde", "CYS"): "GOOD",
    ("aldehyde", "SER"): "GOOD",
    ("aldehyde", "LYS"): "GOOD",

    # α-ketoamide - similar to aldehyde, reversible hemithioacetal
    ("alpha_ketoamide", "CYS"): "GOOD",
    ("alpha_ketoamide", "SER"): "GOOD",
    ("alpha_ketoamide", "LYS"): "GOOD",

    # Isothiocyanate - CYS/LYS/HIS
    ("isothiocyanate", "CYS"): "GOOD",
    ("isothiocyanate", "LYS"): "GOOD",
    ("isothiocyanate", "HIS"): "GOOD",
    ("isothiocyanate", "SER"): "SLOW",

    # Acyl fluoride - very reactive
    ("acyl_fluoride", "CYS"): "GOOD",
    ("acyl_fluoride", "SER"): "GOOD",
    ("acyl_fluoride", "THR"): "GOOD",
    ("acyl_fluoride", "TYR"): "GOOD",
    ("acyl_fluoride", "LYS"): "GOOD",
    ("acyl_fluoride", "HIS"): "GOOD",
}


def check_warhead_residue_compatibility(
    warhead_type: str,
    residue_name: str,
    strict: bool = False
) -> tuple[bool, str]:
    """
    Check if warhead and residue are chemically compatible.

    Args:
        warhead_type: Type of warhead (e.g., "acrylamide")
        residue_name: Type of residue (e.g., "CYS")
        strict: If True, reject SLOW reactions. If False, only reject NO.

    Returns:
        (is_compatible, message)
    """
    key = (warhead_type, residue_name)
    compat = WARHEAD_RESIDUE_COMPATIBILITY.get(key, "SLOW")

    if compat == "NO":
        return False, f"{warhead_type} is NOT compatible with {residue_name} (no reaction expected)"
    elif compat == "SLOW" and strict:
        return False, f"{warhead_type} with {residue_name} may be slow or rare (use strict=False to allow)"
    elif compat == "SLOW":
        return True, f"Warning: {warhead_type} with {residue_name} may be slow (reaction possible but uncommon)"
    else:  # GOOD
        return True, f"{warhead_type} with {residue_name}: well-established combination"


# ---------------------------------------------------------------------------
# Reactive residue definitions (extensible per-residue config)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ResidueConfig:
    """Configuration for a reactive residue type."""
    atom_name: str       # PDB atom name of the nucleophilic atom (e.g. "SG")
    bond_length: float   # Typical covalent bond length to warhead (Å)


# Reactive residues and their nucleophilic atoms
# Ordered by frequency of use in covalent drug discovery
REACTIVE_RESIDUES: dict[str, ResidueConfig] = {
    # === Core residues (most common) ===

    # Thiol (strongest nucleophile, >90% of covalent drugs)
    "CYS": ResidueConfig(atom_name="SG", bond_length=1.82),  # C-S bond

    # Hydroxyl (serine proteases: trypsin, elastase, thrombin)
    "SER": ResidueConfig(atom_name="OG", bond_length=1.43),  # C-O bond

    # Amine (HDAC inhibitors, lysine-targeting drugs)
    "LYS": ResidueConfig(atom_name="NZ", bond_length=1.47),  # C-N bond

    # === Optional residues (less common but useful) ===

    # Hydroxyl (similar to SER, some protease variants)
    "THR": ResidueConfig(atom_name="OG1", bond_length=1.43), # C-O bond

    # Phenolic hydroxyl (protein tyrosine phosphatase inhibitors)
    "TYR": ResidueConfig(atom_name="OH", bond_length=1.43),  # C-O bond (phenolic)

    # Imidazole (moderate nucleophile, some metalloproteases)
    "HIS": ResidueConfig(atom_name="NE2", bond_length=1.47), # C-N bond (epsilon nitrogen)

    # Note: ASP/GLU removed - carboxylates are poor nucleophiles
    # (negative charge repels electrophiles, extremely rare in literature)
}


# ---------------------------------------------------------------------------
# Warhead detection
# ---------------------------------------------------------------------------

@dataclass
class WarheadHit:
    """Result of warhead detection on a ligand."""
    warhead_type: str
    reactive_atom_idx: int
    matched_atoms: tuple[int, ...]


def detect_warheads(mol: Chem.Mol) -> list[WarheadHit]:
    """
    Detect all warhead groups present in a molecule.

    Returns a list of WarheadHit sorted by specificity (most atoms matched first).
    """
    hits: list[WarheadHit] = []
    seen_reactive: set[int] = set()

    for smarts, map_num, name in _WARHEAD_REGISTRY:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue

        # Find which pattern atom carries the map number
        reactive_pattern_idx = None
        for atom in pattern.GetAtoms():
            if atom.GetAtomMapNum() == map_num:
                reactive_pattern_idx = atom.GetIdx()
                break
        if reactive_pattern_idx is None:
            continue

        for match in matches:
            reactive_idx = match[reactive_pattern_idx]
            if reactive_idx in seen_reactive:
                continue
            seen_reactive.add(reactive_idx)
            hits.append(WarheadHit(
                warhead_type=name,
                reactive_atom_idx=reactive_idx,
                matched_atoms=match,
            ))

    # Prefer more specific matches (larger SMARTS → more atoms)
    hits.sort(key=lambda h: len(h.matched_atoms), reverse=True)
    return hits


# ---------------------------------------------------------------------------
# Reactive residue extraction from protein PDB
# ---------------------------------------------------------------------------

@dataclass
class AnchorPoint:
    """Protein-side anchor for covalent bond placement."""
    residue_name: str
    residue_num: int
    chain_id: str
    atom_name: str
    coord: np.ndarray           # [3] xyz of nucleophilic atom (S for CYS)
    bond_vector: np.ndarray     # [3] unit vector along Cβ→Sγ (approach direction)
    bond_length: float          # expected covalent bond length (Å)
    cb_coord: Optional[np.ndarray] = None  # [3] xyz of Cβ atom (for chain rotation)


def find_reactive_residues(pocket_mol: Chem.Mol,
                           residue_spec: Optional[str] = None) -> list[AnchorPoint]:
    """
    Find reactive residues in a protein pocket.

    Args:
        pocket_mol: RDKit Mol loaded from PDB
        residue_spec: Optional residue specifier, e.g. "CYS145", "CYS145:A".
                      If None, auto-detect all supported reactive residues.

    Returns:
        List of AnchorPoint objects.
    """
    # Parse residue_spec if provided
    target_resname = None
    target_resnum = None
    target_chain = None

    if residue_spec is not None:
        parts = residue_spec.split(":")
        res_part = parts[0]
        if len(parts) > 1:
            target_chain = parts[1]

        # Extract residue name and number: "CYS145" -> ("CYS", 145)
        alpha_part = ""
        num_part = ""
        for ch in res_part:
            if ch.isdigit():
                num_part += ch
            else:
                alpha_part += ch
        target_resname = alpha_part.upper()
        target_resnum = int(num_part) if num_part else None

    conf = pocket_mol.GetConformer()
    anchors: list[AnchorPoint] = []

    # Group atoms by residue
    residue_atoms: dict[tuple[str, int, str], dict[str, int]] = {}
    for atom in pocket_mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        resname = info.GetResidueName().strip()
        resnum = info.GetResidueNumber()
        chain = info.GetChainId().strip()
        atom_name = info.GetName().strip()

        key = (resname, resnum, chain)
        if key not in residue_atoms:
            residue_atoms[key] = {}
        residue_atoms[key][atom_name] = atom.GetIdx()

    for (resname, resnum, chain), atom_map in residue_atoms.items():
        # Filter by user-specified residue
        if target_resname is not None and resname != target_resname:
            continue
        if target_resnum is not None and resnum != target_resnum:
            continue
        if target_chain is not None and chain != target_chain:
            continue

        # Check if this residue type is supported
        res_cfg = REACTIVE_RESIDUES.get(resname)
        if res_cfg is None:
            continue

        nuc_atom_name = res_cfg.atom_name
        if nuc_atom_name not in atom_map:
            continue

        nuc_idx = atom_map[nuc_atom_name]
        nuc_pos = np.array(conf.GetAtomPosition(nuc_idx))

        # Compute bond approach vector (CB → nucleophile) and get CB coordinates
        # All standard amino acids have CB (except GLY)
        bond_vec = np.array([0.0, 0.0, 1.0])  # fallback
        cb_coord = None
        if "CB" in atom_map:
            cb_pos = np.array(conf.GetAtomPosition(atom_map["CB"]))
            cb_coord = cb_pos  # Store CB position
            vec = nuc_pos - cb_pos
            norm = np.linalg.norm(vec)
            if norm > 1e-6:
                bond_vec = vec / norm

        anchors.append(AnchorPoint(
            residue_name=resname,
            residue_num=resnum,
            chain_id=chain,
            atom_name=nuc_atom_name,
            coord=nuc_pos,
            bond_vector=bond_vec,
            bond_length=res_cfg.bond_length,
            cb_coord=cb_coord,
        ))

    return anchors


# ---------------------------------------------------------------------------
# Coordinate constraint generation
# ---------------------------------------------------------------------------

def create_covalent_coordmap(
    cb_atom_idx: Optional[int],
    s_atom_idx: int,
    anchor: AnchorPoint
) -> dict[int, Point3D]:
    """
    Create a coordMap that fixes CB and S atoms at protein anchor positions.

    This allows ligand conformers to explore diverse orientations around
    the fixed CB-S anchor point, maximizing conformational diversity.

    Args:
        cb_atom_idx: Index of CB atom in adduct molecule (or None if not present)
        s_atom_idx: Index of S atom in adduct molecule
        anchor: Protein anchor point (contains CB and S coordinates)

    Returns:
        coordMap: {atom_idx: Point3D} for use with EmbedMultipleConfs.
                  Contains CB and S positions only (NOT reactive atom).
    """
    coord_map = {}

    # Fix CB position (if present)
    if cb_atom_idx is not None and anchor.cb_coord is not None:
        coord_map[cb_atom_idx] = Point3D(
            float(anchor.cb_coord[0]),
            float(anchor.cb_coord[1]),
            float(anchor.cb_coord[2]),
        )

    # Fix S position
    coord_map[s_atom_idx] = Point3D(
        float(anchor.coord[0]),
        float(anchor.coord[1]),
        float(anchor.coord[2]),
    )

    return coord_map
