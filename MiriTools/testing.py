#-*- coding: utf-8 -*-
#!/usr/bin/env python

"""

Package testing provides utilities for running tests

== OBSOLETE MODULE (replaced by nosetests or pytest collection) ==

:History:

29 Nov 2010: Created with function discover_tests
12 Mar 2012: Corrected typos.
26 Apr 2013: Test list is sorted to run modules in alphabetical order.
12 Feb 2014: Removed the hard-wired 'miri.tools.tests'. Close the
             open files from which the test modules have been imported.

@author: Julien Morin (DIAS)

"""

from __future__ import division

import glob
import imp
import inspect
import os
import unittest


def discover_tests(test_package=None, top_level='miri.tests.'):
    """
    
    Discovers the TestCase classes contained in the package test_package and
    returns a TestSuite.
    
    For an example of how to use this utility, see the test() functions
    contained in a variety of MIRI lib/__init__.py modules.
    
    NOTE: In Python >=2.7 unittest provides a test discovery utility, and this
    functionality has been backported in the unittest2 package available for
    Python >=2.4. The discover_tests function provides a lighter/simpler
    implementation and avoids unnecessary dependency.
    
    :Parameters:
    
    test_package: Python package
        A Python package to be searched for test modules.
        If no package is provided, no tests will be run.
    top_level: str (optional)
        The name of the top level tests package. This package must exist.
        The name defaults to the miri.tools.tests package.

    """
#     print("Discover_tests: " + str(test_package))
    if test_package is None:
        # No tests will be run. This allows nosetests to skip over this module.
        return
    # Locate all the files named with the pattern 'test_*.py' found
    # within the test package.
    tests_modules = glob.glob(\
        os.path.join(test_package.__path__[0], 'test_*.py'))
    # Sort the tests into alphabetical order
    tests_modules.sort()
#     print("Test modules found: " + str(tests_modules))

    # Step through each of the test modules in turn and use the Python
    # importer (imp) to find and load each module.
    # NOTE: imp.load_module cares more about the file object (mod_descr[0])
    # rather than the name of the module (top_level + mod_name).
    tests_list = []
    tests_names = []
    for mod in tests_modules: 
        mod_name = (os.path.basename(mod).split('.')[0])
#         print("Analysing", mod_name)
        mod_descr = imp.find_module(mod_name, test_package.__path__)
        try:
            mod_obj = imp.load_module(top_level + mod_name, mod_descr[0],
                                      mod_descr[1], mod_descr[2])
            # Inspect the test module, find all the test case instances and
            # append them to the tests_list.
            for name, obj in inspect.getmembers(mod_obj):
                if isinstance(obj, type) and issubclass(obj, unittest.TestCase):
                    tests_names.append(name)
                    tests_list.append(\
                        unittest.TestLoader().loadTestsFromTestCase(obj))
        except:
            # Close each module when finished.
            if mod_descr[0] is not None:
                mod_descr[0].close()
    # Finally, execute the list of test cases and return the results.
#     print("Executing: " + str(tests_names))
    return unittest.TestSuite(tests_list)
