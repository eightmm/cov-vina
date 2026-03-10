"""Shared test fixtures: minimal protein pocket + reference ligand."""

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem


# Minimal pocket PDB with a few residues
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
ATOM     16  N   GLU A 166       8.000   3.000   1.000  1.00  0.00           N
ATOM     17  CA  GLU A 166       9.000   4.000   1.000  1.00  0.00           C
ATOM     18  C   GLU A 166      10.000   3.500   2.000  1.00  0.00           C
ATOM     19  O   GLU A 166      10.500   2.400   2.000  1.00  0.00           O
ATOM     20  CB  GLU A 166       9.500   4.500   -0.300  1.00  0.00          C
END
"""


@pytest.fixture
def pocket_pdb(tmp_path):
    """Create a minimal pocket PDB file."""
    pdb_file = tmp_path / "pocket.pdb"
    pdb_file.write_text(_MINI_POCKET_PDB)
    return str(pdb_file)


@pytest.fixture
def ref_ligand_sdf(tmp_path):
    """Create a small reference ligand SDF file (phenylacetic acid)."""
    mol = Chem.MolFromSmiles("OC(=O)Cc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)

    sdf_file = tmp_path / "ref_ligand.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol)
    writer.close()
    return str(sdf_file)
