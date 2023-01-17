#!/usr/bin/env python
##############################################################################
#
# isotools            by Frandsen Group
#                     Benjamin A. Frandsen benfrandsen@byu.edu
#                     (c) 2022 Benjamin Allen Frandsen
#                     All rights reserved
#
# File coded by:    Benjamin Frandsen
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Tools for adapting ISODISTORT outputs to DiffPy fits.
"""

from diffpy.srfit.fitbase.parameterset import ParameterSet

def getvars(fitnode):
    """Return dictionary of variable names and values in srfit fit object.
    """
    state = {p.name: p.value for p in fitnode}
    return state


def setvars(fitnode, state):
    """Set fitnode variables according to their state dictionary.

    Skip any state items that are constrained.
    """
    for n, v in dict(state).items():
        p = fitnode.get(n)
        p.setValue(v)
    return


def loadmodes(fit,iso,phase):
    """Turns symmetry mode amplitudes into fit variables.
    
    This function creates a ParameterSet consisting of the symmetry mode
    amplitudes, constrains the xyz positions of the atoms in relation to
    the mode amplitudes, and sets the amplitudes as fitting variables. The
    mode amplitude fitting variables are given the tag 'smd'.
    
    Args:
        fit: diffpy.srfit FitRecipe instance
        iso: IsoInfo instance from iso2diffpy
        phase: diffpy.structure parameters generated from the structure
            that was added to the PDF contribution.
    """
    smd = ParameterSet('smd')
    phase.addParameterSet(smd)
    for n, v in iso.modepars.items():
        smd.newParameter(n, v)
    for n, ev in iso.deltas.items():
        smd.newParameter(n, ev[-1])
        smd.constrain(n, ev[0])

    atoms = phase.getScatterers()
    for idx, (v, ev) in enumerate(iso.positions.items()):
        atomIdx, mod = divmod(idx, 3)
        if mod == 0:
            smd.constrain(atoms[atomIdx].x,ev[0])
        if mod == 1:
            smd.constrain(atoms[atomIdx].y,ev[0])
        if mod == 2:
            smd.constrain(atoms[atomIdx].z,ev[0])

    for pa in smd.iterPars('a\d+'):
        fit.addVar(pa,tag='smd')
