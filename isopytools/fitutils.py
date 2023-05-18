#!/usr/bin/env python
##############################################################################
#
# isopytools            by Frandsen Group
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
from collections import OrderedDict
import diffpy.structure.symmetryutilities as su
import diffpy.structure.spacegroups as sg
from isopytools.iso2diffpy import get_positions
import numpy as np
import sympy as sp

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

def generate_constraints(iso):
    """Generate general expressions for the position of each atom.
    
    Args:
        iso: IsoInfo instance from iso2diffpy
       
    Returns:
        Dictionary with coordinates of each atom expressed in terms of the
        deltas (which can in turn be expressed in terms of the mode
        amplitudes).
    """
    names, coords = get_positions(iso, returnDict=False)
    sgobj = sg.GetSpaceGroup(iso.spacegroupnum)
    alltransformations = []
    uniqueidxs = []
    # go through list of symmetry operations once to sort out the duplicates
    for idx, name in enumerate(names):
        eqpos, symoplist, mult = su.expandPosition(sgobj, coords[idx])
        allpositions = []
        for sublist in symoplist:
            for symop in sublist:
                newposN = (np.dot(symop.R, coords[idx]) + symop.t) % 1
                allpositions.append(newposN)
        allpositions = np.round(np.array(allpositions), decimals=6)
        unique, idxs = np.unique(np.round(allpositions, decimals=4), axis=0, return_index=True)
        uniqueidxs.append(idxs)
    # now go through it again to apply the symmetry transformations to the mode shifts
    counter = 0
    for idx, name in enumerate(names):
        pos = sp.Matrix([iso.positions[name+'_x'][0],iso.positions[name+'_y'][0],iso.positions[name+'_z'][0]])
        eqpos, symoplist, mult = su.expandPosition(sgobj, coords[idx])
        partialtransformations = []
        for sublist in symoplist:
            for symop in sublist:
                Rmat = sp.Matrix(symop.R)
                tmat = sp.Matrix(symop.t)
                newpos = Rmat*pos + tmat
                partialtransformations.append(newpos)
        subset = [partialtransformations[i] for i in uniqueidxs[counter]]
        for vector in subset:
            for idx, element in enumerate(vector):
                try:
                    float(element)
                    vector[idx] = element % 1
                except TypeError:
                    pass

        alltransformations.append(subset)
        counter += 1
    # create a dictionary with the transformation constriants
    counter = 0
    transformed_positions = OrderedDict()
    for tr in alltransformations:
        for transformed_pos in tr:
            counter += 1
            transformed_positions['atom'+str(counter)+'_x'] = str(transformed_pos[0])
            transformed_positions['atom'+str(counter)+'_y'] = str(transformed_pos[1])
            transformed_positions['atom'+str(counter)+'_z'] = str(transformed_pos[2])
    return transformed_positions

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
    transform_dict = generate_constraints(iso)
    for idx, (v, ev) in enumerate(transform_dict.items()):
        atomIdx, mod = divmod(idx, 3)
        if mod == 0:
            smd.constrain(atoms[atomIdx].x,ev)
        if mod == 1:
            smd.constrain(atoms[atomIdx].y,ev)
        if mod == 2:
            smd.constrain(atoms[atomIdx].z,ev)

    for pa in smd.iterPars('a\d+'):
        fit.addVar(pa,tag='smd')
