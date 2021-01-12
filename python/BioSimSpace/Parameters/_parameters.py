######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2020
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Functionality for parameterising molecules.
"""

__author__ = "Lester Hedges"
__email_ = "lester.hedges@gmail.com"

__all__ = ["parameterise",
           "ff99",
           "ff99SB",
           "ff99SBildn",
           "ff14SB",
           "gaff",
           "gaff2",
           "forceFields",
           "openForceFields"]

from BioSimSpace import _amber_home, _gmx_exe, _gromacs_path

from BioSimSpace._Exceptions import MissingSoftwareError as _MissingSoftwareError
from BioSimSpace._SireWrappers import Molecule as _Molecule
from BioSimSpace.Types import Charge as _Charge

from ._process import Process as _Process
from . import Protocol as _Protocol

def parameterise(molecule, forcefield, work_dir=None, property_map={}):
    """Parameterise a molecule using a specified force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       forcefield : str
           The force field. Run BioSimSpace.Parameters.forceFields() to get a
           list of the supported force fields.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if type(forcefield) is not str:
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        forcefield = forcefield.replace(" ", "").lower()

        if forcefield not in _forcefields_lower:
            raise ValueError("Supported force fields are: %s" % forceFields())

    return _forcefield_dict[forcefield](molecule, work_dir=work_dir, property_map=property_map)

def ff99(molecule, work_dir=None, property_map={}):
    """Parameterise using the ff99 force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99' is not supported. "
                                    "Please install AMBER (http://ambermd.org) or "
                                    "GROMACS (http://www.gromacs.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99(property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff99SB(molecule, work_dir=None, property_map={}):
    """Parameterise using the ff99SB force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99SB' is not supported. "
                                    "Please install AMBER (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99SB(property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff99SBildn(molecule, work_dir=None, property_map={}):
    """Parameterise using the ff99SBildn force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff99SBildn' is not supported. "
                                    "Please install AMBER (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF99SBILDN(property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff03(molecule, work_dir=None, property_map={}):
    """Parameterise using the ff03 force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None and (_gmx_exe is None or _gromacs_path is None):
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff03' is not supported. "
                                    "Please install AMBER (http://ambermd.org) "
                                    "or GROMACS (http://www.gromacs.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF03(property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def ff14SB(molecule, work_dir=None, property_map={}):
    """Parameterise using the ff14SB force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.ff14SB' is not supported. "
                                    "Please install AMBER (http://ambermd.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.FF14SB(property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def gaff(molecule, work_dir=None, net_charge=None, property_map={}):
    """Parameterise using the gaff force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
           The net charge on the molecule.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.gaff' is not supported. "
                                    "Please install AMBER (http://ambermd.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if net_charge is not None:
        # Get the magnitude of the charge.
        if type(net_charge) is _Charge:
            net_charge = net_charge.magnitude()

        if type(net_charge) is float:
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF(net_charge=net_charge, property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def gaff2(molecule, work_dir=None, net_charge=None, property_map={}):
    """Parameterise using the gaff force field.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
           The net charge on the molecule.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError("'BioSimSpace.Parameters.gaff2' is not supported. "
                                    "Please install AMBER (http://ambermd.org).")

    # Validate arguments.

    if type(molecule) is not _Molecule:
        raise TypeError("'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if net_charge is not None:
        # Get the magnitude of the charge.
        if type(net_charge) is _Charge:
            net_charge = net_charge.magnitude()

        if type(net_charge) is float:
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError("'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'")

        if net_charge % 1 != 0:
            raise ValueError("'net_charge' must be integer valued.")

    if type(property_map) is not dict:
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.GAFF2(net_charge=net_charge, property_map=property_map)

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)

def _parameterise_openff(molecule, forcefield, work_dir=None, property_map={}):
    """Parameterise a molecule using a force field from the Open Force Field
       Initiative.

       Parameters
       ----------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to parameterise.

       forcefield : str
           The force field. Run BioSimSpace.Parameters.openForceFields() to get a
           list of the supported force fields.

       work_dir : str
           The working directory for the process.

       property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The parameterised molecule.
    """

    print(forcefield)

# Create a list of the force field names.
# This needs to come after all of the force field functions.
_forcefields = []           # List of force fields (actual names).
_forcefields_lower = []     # List of lower case names.
_forcefield_dict = {}       # Mapping between lower case names and functions.
import sys as _sys
_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var[0].upper() != "P":
        _forcefields.append(_var)
        _forcefields_lower.append(_var.lower())
        _forcefield_dict[_var.lower()] = getattr(_namespace, _var)
del _sys
del _var

# Wrapper function to dynamically generate functions with a given name.
# Here "name" refers to the name of a supported force field from the Open
# Force Field Initiative. The force field name has been "tidied" so that
# it conforms to sensible function naming standards, i.e. "-" and "."
# characters replaced by underscores.
def _make_function(name):
    def _function(molecule, work_dir=None, property_map={}):
        _parameterise_openff(molecule, name, work_dir, property_map)
    return _function

# Dynamically create functions for all available force fields from the Open
# Force Field Initiative.
from glob import glob as _glob
import openforcefields as _openforcefields
import os as _os
_openff_dirs = _openforcefields.get_forcefield_dirs_paths()
_open_forcefields = []
# Loop over all force field directories.
for _dir in _openff_dirs:
    # Glob all offxml files in the directory.
    _ff_list = _glob(f"{_dir}" + "/*.offxml")
    for _ff in _ff_list:
        # Get the force field name (base name minus extension).
        _base = _os.path.basename(_ff)
        _ff = _os.path.splitext(_base)[0]
        # Append to the list of available force fields.
        _forcefields.append(_ff)
        _open_forcefields.append(_ff)

        # Create a sane function name, i.e. replace "-" and "."
        # characters with "_".
        _func_name = _ff.replace("-", "_")
        _func_name = _func_name.replace(".", "_")

        # Generate the function and bind it to the namespace.
        _function = _make_function(_func_name)
        setattr(_namespace, _func_name, _function)

        # Expose the function to the user.
        __all__.append(_func_name)

        # Convert force field name to lower case and map to its function.
        _forcefields_lower.append(_ff.lower())
        _forcefield_dict[_ff.lower()] = getattr(_namespace, _func_name)
del _glob
del _openforcefields
del _os
del _namespace

def forceFields():
    """Return a list of the supported force fields.

       Returns
       -------

       force_fields : [str]
           A list of the supported force fields.
    """
    return _forcefields

def openForceFields():
    """Return a list of the supported force fields from the Open Force Field
       Initiative.

       Returns
       -------

       force_fields : [str]
           A list of the supported force fields.
    """
    return _open_forcefields
