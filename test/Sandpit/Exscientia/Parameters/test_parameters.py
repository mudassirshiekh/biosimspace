import BioSimSpace.Sandpit.Exscientia as BSS

import os
import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()

# Make sure required AMBER executables are present.
if BSS._amber_home is not None:
    tleap = "%s/bin/tleap" % BSS._amber_home
    if os.path.isfile(tleap):
        has_tleap = True
    else:
        has_tleap = False
else:
    has_tleap = False


@pytest.fixture(scope="session")
def molecule0():
    return BSS.IO.readMolecules(f"{url}/4LYT_Fixed.pdb.bz2")[0]


@pytest.fixture(scope="session")
def molecule1():
    return BSS.IO.readMolecules(f"{url}/3G8K_Fixed.pdb.bz2")[0]


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
@pytest.mark.parametrize("ff", BSS.Parameters.amberProteinForceFields())
def test_disulphide(molecule0, ff):
    """Test parameterisation in the presence of disulphide bridges."""

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule0).getMolecule()

    # Check that we actually generate records for four disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule0._sire_object
    )
    assert len(bonds) == 4

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 4


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
def test_disulphide_renumber(molecule1, ff="ff14SB"):
    """
    Test parameterisation in the presence of disulphide bridges using a
    multi-chain PDB with duplicate residue numbering
    """

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule1).getMolecule()

    # Check that we actually generate records for eight disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule1._sire_object
    )
    assert len(bonds) == 8

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 8
