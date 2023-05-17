import BioSimSpace.Sandpit.Exscientia as BSS

import glob
import os
import pytest
import tempfile


def test_file_cache():
    """
    Simple test to see if cached files are used when repeatedly writing
    the same system to the same file format.
    """

    # Clear the file cache.
    BSS.IO._file_cache._cache = {}

    # Load the molecular system.
    s = BSS.IO.readMolecules(["tests/input/ala.crd", "tests/input/ala.top"])

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name

    # Write the system to PDB and PRM7 format.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp", s, ["pdb", "prm7"])

    # Check that the file cache has two entries.
    assert len(BSS.IO._file_cache._cache) == 2

    # Write to PDB and GroTop format. The PDB from the cache should be reused.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp2", s, ["pdb", "grotop"])

    # Check that the file cache has three entries, i.e. only the GroTop was added.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 4 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 4

    # Now delete one of the files on disk and re-write. This should create a new
    # entry in the cache.
    os.remove(f"{tmp_path}/tmp.pdb")
    BSS.IO.saveMolecules(f"{tmp_path}/tmp3", s, ["pdb", "grotop"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 5 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 5

    # Now "corrupt" a file on disk so that its MD5 checksum is no longer
    # valid.
    with open(f"{tmp_path}/tmp2.pdb", "w") as f:
        pass

    # Write back to PDB and Gro87 format. The PDB file is now invalid, so
    # a new one will be written and added to the cache.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp4", s, ["pdb", "gro87"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 4

    # The directory should now contain 7 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 7


@pytest.mark.parametrize(
    "ext,is_cached",
    [
        ("prm7", True),
        ("rst7", False),
        ("grotop", True),
        ("gro87", False),
        ("pdb", False),
    ],
)
@pytest.mark.parametrize("cache_topology", [True, False])
def test_topology_cache(tmp_path_factory, ext, is_cached, cache_topology):
    # Load the molecular system.
    s = BSS.IO.readMolecules(["tests/input/ala.crd", "tests/input/ala.top"])
    assert s._cache_topology is True
    s._cache_topology = cache_topology
    assert len(s._topology_cache) == 0

    # Save.
    workdir = tmp_path_factory.mktemp("test_topology_cache")
    BSS.IO.saveMolecules(f"{workdir}/system", s, ext)

    # Check if we have cached the topology.
    if is_cached and cache_topology:
        assert len(s._topology_cache) == 1
    else:
        assert len(s._topology_cache) == 0
        return

    # Now we replace the cache to see if we are using it when we write again.
    (key,) = s._topology_cache
    s._topology_cache[key] = ["dummy\n", "value\n"]

    # We save again.
    (out_file,) = BSS.IO.saveMolecules(f"{workdir}/system", s, ext)
    assert open(out_file).readlines() == ["dummy\n", "value\n"]

    # Now we modify the system - this should reset the cache.
    s.removeMolecules(s[0])
    assert len(s._topology_cache) == 0
