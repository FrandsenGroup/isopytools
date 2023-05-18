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

"""Unit tests basic pydistort functionalities.
"""


import unittest
import isotools
import os
import sys

##############################################################################
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

class testIso(unittest.TestCase):
    def test_iso(self):
        path = os.path.dirname(os.path.abspath(__file__))
        fname = find('MnTe_iso.txt', path)
        iso = isotools.iso2diffpy.IsoInfo(fname)
        teststr = iso.deltas['Te_1_dx'][0]
        self.assertEqual(teststr, '+  0.09736*a2 -  0.16864*a3 +  0.09736*a5 -  0.16864*a6')
        
    def test_struc(self):
        path = os.path.dirname(os.path.abspath(__file__))
        fname = find('MnTe_iso.txt', path)
        iso = isotools.iso2diffpy.IsoInfo(fname)
        struc = isotools.iso2diffpy.build_struc(iso)
        self.assertEqual(struc[0].element, 'Te')

# End of class

if __name__ == '__main__':
    unittest.main()

# End of file


