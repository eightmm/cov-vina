"""Integration test for the covalent docking pipeline."""

import os
import pytest
from rdkit import Chem

from cov_vina import run_covalent_pipeline


# A minimal PDB with CYS145 + a few surrounding residues for pocket
_MINI_POCKET_PDB = """\
ATOM      1  N   GLY A 144       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A 144       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A 144       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   GLY A 144       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  N   CYS A 145       3.326   1.500   0.000  1.00  0.00           N
ATOM      6  CA  CYS A 145       3.950   2.820   0.000  1.00  0.00           C
ATOM      7  CB  CYS A 145       5.470   2.730   0.000  1.00  0.00           C
ATOM      8  SG  CYS A 145       6.120   1.070   0.000  1.00  0.00           S
ATOM      9  C   CYS A 145       3.500   3.700   1.160  1.00  0.00           C
ATOM     10  O   CYS A 145       3.950   4.850   1.200  1.00  0.00           O
ATOM     11  N   HIS A 163       5.000   5.000   2.000  1.00  0.00           N
ATOM     12  CA  HIS A 163       5.500   6.300   2.000  1.00  0.00           C
ATOM     13  C   HIS A 163       7.000   6.400   2.000  1.00  0.00           C
ATOM     14  O   HIS A 163       7.600   7.400   2.000  1.00  0.00           O
ATOM     15  CB  HIS A 163       5.000   7.200   3.100  1.00  0.00           C
END
"""

# A simple acrylamide-bearing ligand
_ACRYLAMIDE_SMILES = "C=CC(=O)Nc1ccccc1"


@pytest.fixture
def pocket_pdb(tmp_path):
    pdb_file = tmp_path / "pocket.pdb"
    pdb_file.write_text(_MINI_POCKET_PDB)
    return str(pdb_file)


class TestCovalentPipeline:

    def test_basic_run(self, pocket_pdb, tmp_path):
        """Smoke test: pipeline runs end-to-end and produces an SDF."""
        output_dir = str(tmp_path / "output")
        results = run_covalent_pipeline(
            protein_pdb=pocket_pdb,
            query_ligand=_ACRYLAMIDE_SMILES,
            reactive_residue="CYS145",
            output_dir=output_dir,
            num_confs=50,
            device="cpu",
            verbose=False,
        )
        assert os.path.exists(results["output_file"])
        assert results["num_poses"] > 0
        assert results["warhead_type"] == "acrylamide"
        assert results["anchor_residue"] == "CYS145"
        assert isinstance(results["best_score"], float)

    def test_auto_detect_residue(self, pocket_pdb, tmp_path):
        """Without specifying residue, should auto-detect CYS145."""
        output_dir = str(tmp_path / "output")
        results = run_covalent_pipeline(
            protein_pdb=pocket_pdb,
            query_ligand=_ACRYLAMIDE_SMILES,
            reactive_residue=None,
            output_dir=output_dir,
            num_confs=50,
            device="cpu",
            verbose=False,
        )
        assert results["anchor_residue"] == "CYS145"

    def test_no_warhead_raises(self, pocket_pdb, tmp_path):
        """A ligand with no warhead should raise ValueError."""
        with pytest.raises(ValueError, match="No reactive warhead"):
            run_covalent_pipeline(
                protein_pdb=pocket_pdb,
                query_ligand="c1ccccc1",  # benzene
                output_dir=str(tmp_path / "output"),
                num_confs=10,
                device="cpu",
                verbose=False,
            )

    def test_wrong_residue_raises(self, pocket_pdb, tmp_path):
        """Specifying a non-existent residue should raise ValueError."""
        with pytest.raises(ValueError, match="No reactive residue"):
            run_covalent_pipeline(
                protein_pdb=pocket_pdb,
                query_ligand=_ACRYLAMIDE_SMILES,
                reactive_residue="CYS999",
                output_dir=str(tmp_path / "output"),
                num_confs=10,
                device="cpu",
                verbose=False,
            )

    def test_output_sdf_has_metadata(self, pocket_pdb, tmp_path):
        """Output SDF should contain covalent docking metadata."""
        output_dir = str(tmp_path / "output")
        results = run_covalent_pipeline(
            protein_pdb=pocket_pdb,
            query_ligand=_ACRYLAMIDE_SMILES,
            reactive_residue="CYS145",
            output_dir=output_dir,
            num_confs=50,
            device="cpu",
            verbose=False,
        )
        suppl = Chem.SDMolSupplier(results["output_file"])
        mol = suppl[0]
        assert mol is not None
        assert mol.HasProp("CovVina_Warhead_Type")
        assert mol.GetProp("CovVina_Warhead_Type") == "acrylamide"

    def test_with_optimization(self, pocket_pdb, tmp_path):
        """Pipeline with gradient optimization enabled."""
        output_dir = str(tmp_path / "output")
        results = run_covalent_pipeline(
            protein_pdb=pocket_pdb,
            query_ligand=_ACRYLAMIDE_SMILES,
            reactive_residue="CYS145",
            output_dir=output_dir,
            num_confs=20,
            optimize=True,
            opt_steps=10,
            device="cpu",
            verbose=False,
        )
        assert results["num_poses"] > 0
        suppl = Chem.SDMolSupplier(results["output_file"])
        mol = suppl[0]
        assert mol.HasProp("CovVina_Gradient_Optimized")
