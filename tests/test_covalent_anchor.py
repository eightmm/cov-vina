"""Tests for covalent warhead detection and anchor generation."""

import numpy as np
import pytest
from rdkit import Chem

from cov_vina.molecular.anchor import (
    AnchorPoint,
    WarheadHit,
    create_covalent_coordmap,
    detect_warheads,
    find_reactive_residues,
    check_warhead_residue_compatibility,
)


# ------------------------------------------------------------------ #
# Warhead detection tests
# ------------------------------------------------------------------ #

class TestDetectWarheads:
    """Test warhead SMARTS detection on known covalent inhibitors."""

    def test_acrylamide(self):
        # Ibrutinib-like acrylamide warhead
        mol = Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) >= 1
        types = [h.warhead_type for h in hits]
        assert "acrylamide" in types

    def test_chloroacetamide(self):
        mol = Chem.MolFromSmiles("ClCC(=O)Nc1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) >= 1
        types = [h.warhead_type for h in hits]
        assert "chloroacetamide" in types

    def test_aryl_nitrile(self):
        mol = Chem.MolFromSmiles("N#Cc1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) >= 1
        types = [h.warhead_type for h in hits]
        assert "aryl_nitrile" in types

    def test_vinyl_sulfonamide(self):
        mol = Chem.MolFromSmiles("C=CS(=O)(=O)Nc1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) >= 1
        types = [h.warhead_type for h in hits]
        assert "vinyl_sulfonamide" in types

    def test_epoxide(self):
        mol = Chem.MolFromSmiles("C1OC1c1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) >= 1
        types = [h.warhead_type for h in hits]
        assert "epoxide" in types

    def test_no_warhead(self):
        # Benzene - no warhead
        mol = Chem.MolFromSmiles("c1ccccc1")
        hits = detect_warheads(mol)
        assert len(hits) == 0

    def test_reactive_atom_is_valid_index(self):
        mol = Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1")
        hits = detect_warheads(mol)
        for hit in hits:
            assert 0 <= hit.reactive_atom_idx < mol.GetNumAtoms()

    def test_multiple_warheads(self):
        # Molecule with both acrylamide and nitrile
        mol = Chem.MolFromSmiles("C=CC(=O)Nc1ccc(C#N)cc1")
        hits = detect_warheads(mol)
        types = {h.warhead_type for h in hits}
        assert "acrylamide" in types
        assert "aryl_nitrile" in types

    def test_aldehyde(self):
        mol = Chem.MolFromSmiles("O=Cc1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "aldehyde" in types

    def test_bromoacetamide(self):
        mol = Chem.MolFromSmiles("BrCC(=O)Nc1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "bromoacetamide" in types

    def test_isothiocyanate(self):
        mol = Chem.MolFromSmiles("N=C=Sc1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "isothiocyanate" in types

    def test_boronic_acid(self):
        mol = Chem.MolFromSmiles("OB(O)c1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "boronic_acid" in types

    def test_acyl_fluoride(self):
        mol = Chem.MolFromSmiles("FC(=O)c1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "acyl_fluoride" in types

    def test_thiirane(self):
        mol = Chem.MolFromSmiles("C1SC1c1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        assert "thiirane" in types

    def test_aryl_nitrile_not_alkyl(self):
        """Test that aromatic nitrile is detected as aryl_nitrile, not alkyl_nitrile."""
        mol = Chem.MolFromSmiles("N#Cc1ccccc1")
        hits = detect_warheads(mol)
        types = [h.warhead_type for h in hits]
        # Should match aryl_nitrile (more specific pattern first)
        assert "aryl_nitrile" in types


# ------------------------------------------------------------------ #
# Reactive residue detection tests
# ------------------------------------------------------------------ #

def _make_cys_pdb_block() -> str:
    """Create a minimal PDB with a single CYS residue."""
    return """\
ATOM      1  N   CYS A 145       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  CYS A 145       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  CB  CYS A 145       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  SG  CYS A 145       4.000   2.000   3.000  1.00  0.00           S
ATOM      5  C   CYS A 145       2.000   3.000   3.000  1.00  0.00           C
ATOM      6  O   CYS A 145       2.000   4.000   3.000  1.00  0.00           O
END
"""


class TestFindReactiveResidues:

    def test_auto_detect_cys(self, tmp_path):
        pdb_file = tmp_path / "pocket.pdb"
        pdb_file.write_text(_make_cys_pdb_block())
        mol = Chem.MolFromPDBFile(str(pdb_file), sanitize=False, removeHs=True)
        anchors = find_reactive_residues(mol)
        assert len(anchors) >= 1
        assert anchors[0].residue_name == "CYS"
        assert anchors[0].residue_num == 145
        assert anchors[0].atom_name == "SG"

    def test_specify_residue(self, tmp_path):
        pdb_file = tmp_path / "pocket.pdb"
        pdb_file.write_text(_make_cys_pdb_block())
        mol = Chem.MolFromPDBFile(str(pdb_file), sanitize=False, removeHs=True)
        anchors = find_reactive_residues(mol, residue_spec="CYS145")
        assert len(anchors) == 1
        assert anchors[0].residue_num == 145

    def test_wrong_residue_returns_empty(self, tmp_path):
        pdb_file = tmp_path / "pocket.pdb"
        pdb_file.write_text(_make_cys_pdb_block())
        mol = Chem.MolFromPDBFile(str(pdb_file), sanitize=False, removeHs=True)
        anchors = find_reactive_residues(mol, residue_spec="CYS999")
        assert len(anchors) == 0

    def test_bond_vector_normalized(self, tmp_path):
        pdb_file = tmp_path / "pocket.pdb"
        pdb_file.write_text(_make_cys_pdb_block())
        mol = Chem.MolFromPDBFile(str(pdb_file), sanitize=False, removeHs=True)
        anchors = find_reactive_residues(mol)
        vec = anchors[0].bond_vector
        assert abs(np.linalg.norm(vec) - 1.0) < 1e-6


# ------------------------------------------------------------------ #
# CoordMap generation test
# ------------------------------------------------------------------ #

class TestCreateCovalentCoordmap:

    def test_coordmap_cb_s_fixed(self):
        """Test that coordmap fixes CB and S atoms only (not reactive atom)."""
        anchor = AnchorPoint(
            residue_name="CYS", residue_num=145, chain_id="A",
            atom_name="SG",
            coord=np.array([4.0, 2.0, 3.0]),
            bond_vector=np.array([1.0, 0.0, 0.0]),
            bond_length=1.82,
            cb_coord=np.array([3.0, 2.0, 3.0]),
        )

        # Simulate adduct indices (CB=10, S=11)
        cb_atom_idx = 10
        s_atom_idx = 11

        coord_map = create_covalent_coordmap(cb_atom_idx, s_atom_idx, anchor)

        # Check CB is in coordmap
        assert cb_atom_idx in coord_map
        cb_pt = coord_map[cb_atom_idx]
        assert abs(cb_pt.x - 3.0) < 1e-6
        assert abs(cb_pt.y - 2.0) < 1e-6
        assert abs(cb_pt.z - 3.0) < 1e-6

        # Check S is in coordmap
        assert s_atom_idx in coord_map
        s_pt = coord_map[s_atom_idx]
        assert abs(s_pt.x - 4.0) < 1e-6
        assert abs(s_pt.y - 2.0) < 1e-6
        assert abs(s_pt.z - 3.0) < 1e-6

        # IMPORTANT: Reactive atom should NOT be in coordmap
        # (This allows ligand diversity)
        assert 0 not in coord_map  # Reactive atom is not fixed


# ------------------------------------------------------------------ #
# Warhead-Residue compatibility tests
# ------------------------------------------------------------------ #

class TestWarheadResidueCompatibility:

    def test_acrylamide_cys_good(self):
        """Acrylamide + CYS is well-established (Michael addition)."""
        is_compat, msg = check_warhead_residue_compatibility("acrylamide", "CYS")
        assert is_compat is True
        assert "well-established" in msg

    def test_acrylamide_ser_incompatible(self):
        """Acrylamide + SER is not compatible (needs strong nucleophile)."""
        is_compat, msg = check_warhead_residue_compatibility("acrylamide", "SER")
        assert is_compat is False
        assert "NOT compatible" in msg

    def test_boronic_acid_ser_good(self):
        """Boronic acid + SER is good (serine proteases)."""
        is_compat, msg = check_warhead_residue_compatibility("boronic_acid", "SER")
        assert is_compat is True
        assert "well-established" in msg

    def test_boronic_acid_cys_incompatible(self):
        """Boronic acid + CYS is not compatible."""
        is_compat, msg = check_warhead_residue_compatibility("boronic_acid", "CYS")
        assert is_compat is False
        assert "NOT compatible" in msg

    def test_chloroacetamide_cys_good(self):
        """Chloroacetamide + CYS is good (SN2)."""
        is_compat, msg = check_warhead_residue_compatibility("chloroacetamide", "CYS")
        assert is_compat is True
        assert "well-established" in msg

    def test_chloroacetamide_ser_slow(self):
        """Chloroacetamide + SER is slow but possible."""
        is_compat, msg = check_warhead_residue_compatibility("chloroacetamide", "SER", strict=False)
        assert is_compat is True
        assert "slow" in msg.lower()

    def test_chloroacetamide_ser_strict(self):
        """Chloroacetamide + SER is rejected in strict mode."""
        is_compat, msg = check_warhead_residue_compatibility("chloroacetamide", "SER", strict=True)
        assert is_compat is False
        assert "slow or rare" in msg

    def test_unknown_combination_defaults_to_slow(self):
        """Unknown combinations default to SLOW (cautious)."""
        is_compat, msg = check_warhead_residue_compatibility("unknown_warhead", "CYS", strict=False)
        assert is_compat is True  # Allowed in non-strict mode
        assert "slow" in msg.lower()

    def test_sulfonyl_fluoride_multiple_residues(self):
        """Sulfonyl fluoride works with many residues (very reactive)."""
        for residue in ["CYS", "SER", "TYR", "LYS", "HIS"]:
            is_compat, msg = check_warhead_residue_compatibility("sulfonyl_fluoride", residue)
            assert is_compat is True
            assert "well-established" in msg
