#!/usr/bin/env python
##############################################################################
#
# isotools            by Frandsen Group
#                     Benjamin A. Frandsen benfrandsen@byu.edu
#                     (c) 2022 Benjamin Allen Frandsen
#                     All rights reserved
#
# File coded by:    Parker Hamilton and Benjamin Frandsen
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Tools for adapting ISODISTORT outputs to DiffPy.
"""

import re
import itertools
from collections import OrderedDict
import diffpy.structure as dps

_pf = r'[+-]?\d[^; \n\t]*'

def get_xyz(atomName,iso):
    """Extract the xyz position of a given atom from mode amplitudes.

    Args:
        atomName (string): label of the atom for which the coordinates are
            desired.
        iso: instance of IsoInfo class.

    Returns: tuple with the x, y, and z fractional coordinates.
    """
    return tuple(iso.positions[atomName+'_'+c][1] for c in 'xyz')

def get_shifts(iso):
    """Get fractional shifts of symmetry-unique atoms from mode amplitudes.

    Args:
        iso: instance of IsoInfo class.

    Returns:
        shiftDict: OrderedDict with fractional shifts of each atom.
    """
    shiftDict = OrderedDict()
    for key in iso.deltas.keys():
        equation = iso.deltas[key][0]
        pieces = equation.split()
        signs = pieces[::2]
        terms = pieces[1::2]
        coefficients = []
        modenums = []
        for term in terms:
            coefficient, modenum = term.split('*')
            coefficients.append(coefficient)
            modenums.append(modenum)
        modenums2 = []
        for modenum in modenums:
            modenum2 = "iso.modepars['"+modenum+"']"
            modenums2.append(modenum2)
        newequation = ''
        for idx, sign in enumerate(signs):
            newequation += sign+coefficients[idx]+'*'+modenums2[idx]
        shiftDict[key] = eval(newequation)
    return shiftDict

def get_positions(iso):
    """Get fractional positions of symmetry-unique atoms from mode amplitudes.

    Args:
        iso: instance of IsoInfo class.

    Returns:
        posDict: OrderedDict with positions of each atom.
    """
    posDict = OrderedDict()
    shiftDict = get_shifts(iso)
    for key in iso.positions.keys():
        equation = iso.positions[key][0]
        k = key[:-1]+'d'+key[-1]
        pParams = re.findall(k,equation)
        for pParam in pParams:
            equation = equation.replace(pParam,"shiftDict['"+pParam+"']")
        posDict[key] = eval(equation)
    return posDict

def update_positions(strucObj, iso):
    """Update fractional positions of diffpy structure object from mode amps.

    Args:
        strucObj: the diffpy structure object to be updated.
        iso: instance of IsoInfo class.

    Returns: Nothing. Just updates the strucObj.
    """
    posDict = get_positions(iso)
    for idx, key in enumerate(posDict.keys()):
        k = key[:-2]
        XYZ = key[-1]
        string = 'strucObj['+str(idx//3)+'].'+ XYZ +' = posDict[key]'
        exec(string)

def update_positions_pyobjcryst(crystObj, iso):
    """Update fractional positions of pyobjcryst Crystal from mode amps.

    Args:
        crystObj: the pyobjcryst Crystal instance to be updated.
        iso: instance of IsoInfo class.

    Returns: Nothing. Just updates the strucObj.
    """
    posDict = get_positions(iso)
    for key in posDict.keys():
        k = key[:-2]
        XYZ = key[-1].capitalize()
        eval('crystObj.GetScatterer(k).Set'+ XYZ +'(posDict[key])')

def set_amps(iso,ampArray):
    """Set the symmetry mode amplitudes to desired values.

    Args:
        iso: instance of IsoInfo class that will be updated with desired amps
        ampArray: array or list of amplitudes, in the same order as
            iso.modepars.

    Returns: Nothing. Just updates iso.
    """
    for idx, key in enumerate(iso.modepars.keys()):
        iso.modepars[key] = ampArray[idx]

def get_amps(iso):
    """Get the symmetry mode amplitudes from current iso instance.

    Args:
        iso: instance of IsoInfo class from which amplitudes will be extracted.
        ampArray: array or list of amplitudes, in the same order as
            iso.modepars.

    Returns: list containing the mode amplitudes.
    """
    amps = []
    for idx, key in enumerate(iso.modepars.keys()):
        amps.append(iso.modepars[key])
    return amps

def build_struc(iso):
    """Create a PDFfitStructure instance from the IsoInfo instance.
    
    Args:
        iso: instance of IsoInfo class to be converted to a diffpy structure.

    Returns:
        struc: diffpy structure object corresponding to info in IsoInfo.
    """
    # create lattice
    a, b, c, alpha, beta, gamma = iso.abcdegabg 
    lat = dps.Lattice(a, b, c, alpha, beta, gamma)
    
    # create the atoms
    atoms = []
    for key in list(iso.positions.keys())[::3]: # grab every third entry because x, y, z are listed for each atom
        letters = re.findall(r"[a-zA-Z]",key)
        elem = ''.join(letters[:-1])
        atoms.append(dps.Atom(atype=elem, xyz=get_xyz(key[:-2],iso), anisotropy=True))

    struc = dps.Structure(atoms=atoms)
    struc.lattice = lat
    return struc

class IsoInfo:
    """Container for structural information contained in ISODISTORT output file.

    This class takes a *.str TOPAS-format output file from ISODISTORT,
    extracts the structural/symmetry information, and prepares it for use by
    other DiffPy-based functions.
    
    Args:
        filename (string): name of the .str ISODISTORT file.
    """
    def __init__(self, filename):
        from math import radians
        with open(filename) as fp:
            lines = fp.readlines()
        self.abcdegabg = self._parse_abcabg(lines)
        radabg = tuple(radians(x) for x in self.abcdegabg[3:])
        self.abcradabg = self.abcdegabg[:3] + radabg
        self.modepars = self._parse_modepars(lines)
        self.deltas = self._parse_deltas(lines)
        self.positions = self._parse_positions(lines)
        # find spacegroup
        sgline = next((i for i, line in enumerate(lines)
                       if 'space_group' in line))
        self.spacegroup = lines[sgline - 1].strip().strip("'\"")
        return

    @staticmethod
    def _parse_abcabg(lines):
        """Extract lattice constants.
        """
        mx = re.compile(fr'^\s*([abc]|al|be|ga)\s+({_pf})')
        words = (mx.search(s).group(2) for s in lines if mx.search(s))
        rv = tuple(float(x) for x in itertools.islice(words, 6))
        return rv


    @staticmethod
    def _parse_modepars(lines):
        """Parse the 'mode definitions' block.
        """
        mx = re.compile(fr'^\s*prm *!?(a\d+) +({_pf})')
        words = (mx.search(s).groups() for s in lines if mx.search(s))
        rv = OrderedDict((n, float(v)) for n, v in words)
        return rv


    @staticmethod
    def _parse_deltas(lines):
        """
        Parse the "mode-amplitude to delta transformation" block.
        """
        mx = re.compile(fr'^\s*prm *(\w+_d[xyz]) += +(\S[^;]*);: *({_pf})')
        words = (mx.search(s).groups() for s in lines if mx.search(s))
        rv = OrderedDict((w[0], (w[1], float(w[2]))) for w in words)
        return rv


    @staticmethod
    def _parse_positions(lines):
        """Parse the 'distorted parameters' block.
        """
        mx = re.compile(fr'^\s*prm *!?(\w+_[xyz]) += +(\S[^;]*);: *({_pf})')
        words = (mx.search(s).groups() for s in lines if mx.search(s))
        rv = OrderedDict((w[0], (w[1], float(w[2]))) for w in words)
        return rv

# end of class IsoInfo
