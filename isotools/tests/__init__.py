#!/usr/bin/env python
##############################################################################
#
# isotools       by Billinge Group
#                     Simon J. L. Billinge sb2896@columbia.edu
#                     (c) 2016 trustees of Columbia University in the City of
#                           New York.
#                      All rights reserved
#
# File coded by:    Benjamin Frandsen
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Unit tests for isotools.
"""


# create logger instance for the tests subpackage
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)
del logging


def testsuite():
    '''Build a unit tests suite for the isotools package.
    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        isotools.tests.testtest
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    mobj = None
    for mname in modulenames:
        exec ('import %s as mobj' % mname)
        suite.addTests(loader.loadTestsFromModule(mobj))
    return suite


def test():
    '''Execute all unit tests for the isotools package.
    Return a unittest TestResult object.
    '''
    import unittest
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file
