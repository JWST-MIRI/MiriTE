#!/usr/bin/env python

"""

Module test_filesearching - Contains the unit tests for the
ParameterFileManager class.

:History:

25 Aug 2010: Created
15 Nov 2010: Moved to miri.tools.
16 Jan 2012: Correction: There isn't necessarily a test_filesearching.py
             file in the PYTHONPATH, but there is always a __init__.py
04 Apr 2012: Improvements suggested by pylint.
08 Sep 2015: Made compatible with Python 3
23 Mar 2016: Verbosity parameter replaced by Python logger.
30 Mar 2016: Added a search path parameter to the ParameterFileManager,
             so the search result can be made predictable.
16 Mar 2017: Removed unreliable '.' directory from test_find_files_matching
             function.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
#import numpy as np

# Python logging facility
import logging
logging.basicConfig(level=logging.ERROR) # Turn off logging 
LOGGER = logging.getLogger("miri.tools") # Get a default parent logger

import miri.datamodels
import miri.tools
from miri.tools.filesearching import make_searchpath, find_files_matching, \
    find_file_in_tree, find_file_in_path, find_file_in_pythonpath, find_file, \
    find_file_prefix, ParameterFileManager


class TestParameterFileManager(unittest.TestCase):
    
    def setUp(self):
        # Create a parameters object from the example parameters file.
        dir_list = ['.', os.path.dirname(__file__),
                    miri.datamodels.__path__[0]]
        self.search_path = make_searchpath(dir_list)     
        self.pfmobject = ParameterFileManager('example_properties.py',
                                              search_path=self.search_path,
                                              logger=LOGGER)
        
    def tearDown(self):
        # Tidy up
        del self.pfmobject
        
    def test_creation(self):
        # Check for exceptions when creating bad ParameterFileManager objects.
        # The file name must be non-null and the file much exist.
        self.assertRaises(NameError, ParameterFileManager, '',
                          search_path=self.search_path,
                          logger=LOGGER)
        self.assertRaises(NameError, ParameterFileManager, 'nosuchfile.py',
                          search_path=self.search_path,
                          logger=LOGGER)
        # Check what happens when the file cannot be interpreted
        # by executing a non-Python text file.
        self.assertRaises(Exception, ParameterFileManager,
                          'example_measurement.txt',
                          search_path=self.search_path,
                          logger=LOGGER)

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error
        # and generate non-null strings.
        descr1 = str(self.pfmobject)
        self.assertIsNotNone(descr1)

    def test_lookup(self):
        # Look up a few parameters from the ParameterFileManager object.
        # All the keywords returned by keys() must exist and contain
        # valid values.
        kwlist = self.pfmobject.keys()
        for kw in kwlist:
            self.assertTrue(kw in self.pfmobject)
#            self.assertTrue(self.pfmobject.has_key(kw))
            value = self.pfmobject[kw]
            self.assertIsNotNone(value)

        # Check that valid parameters can be read successfully.
        forty_two = self.pfmobject['INTEGER_PARAMETER']
        self.assertEqual(forty_two, 42)
        pi = self.pfmobject['FLOAT_PARAMETER']
        self.assertAlmostEqual(pi, 3.14159265)
        
        level2 = self.pfmobject.get('DICT_PARAM', 'one')
        self.assertEqual(level2, 1)

        level2 = self.pfmobject.get('TUPLE_PARAMETER', 4)
        self.assertEqual(level2, 5)
        
        level3 = self.pfmobject.get('DETECTORS_DICT', 'Holmes', 'NAME')
        self.assertEqual(level3, 'Holmes')
        
        # Check that an attempt to read a non-existent parameter
        # or a null parameter results in an exception.
        self.assertRaises(KeyError, self.pfmobject.get, 'NO_SUCH_KEYWORD')
        self.assertRaises(KeyError, self.pfmobject.get, '')
        
        # Check that attempting to access a scalar parameter as it if were
        # nested, or a nested parameter with the wrong type of keyword,
        # also results in an exception.
        self.assertRaises(TypeError, self.pfmobject.get, 'INTEGER_PARAMETER',
                          'LOOKTHISUP')
        self.assertRaises(TypeError, self.pfmobject.get, 'INTEGER_PARAMETER',
                          42)
        self.assertRaises((TypeError,KeyError), self.pfmobject.get,
                          'DICT_PARAM', 42)
        self.assertRaises((ValueError,KeyError), self.pfmobject.get,
                          'TUPLE_PARAMETER', 'LOOKTHISUP')
        self.assertRaises(TypeError, self.pfmobject.get, 'DICT_PARAM', 'one',
                          4)
                
    def test_readonly(self):
        # Check that parameters cannot be modified or deleted once
        # read from the file.
        self.assertRaises(TypeError, self.pfmobject.__setitem__,
                          'INTEGER_PARAMETER', 50)
        self.assertRaises(TypeError, self.pfmobject.__delitem__,
                          'INTEGER_PARAMETER')


class TestFileSearchFunctions(unittest.TestCase):
    
    def setUp(self):
        self.test_dir = miri.tools.__path__[0]  # Was '.'
        self.test_file_name = 'this_is_a_test_file.tmp'
        self.test_file_prefix = 'this_is_a_test_'
        self.test_patterns = '*.py;*.tmp'
        
        # Make sure there is a test file in the current working directory.
        fp = open(self.test_file_name, 'w')
        fp.write("Please delete this file.\n")
        fp.close()
        
    def tearDown(self):
        # Tidy up
        os.remove(self.test_file_name)
    
    def test_find_files_matching(self):
        # Check that find_files_matching finds files that exist.
        for mfile in find_files_matching(self.test_dir,
                                         patterns=self.test_patterns):
            self.assertTrue(os.path.isfile(mfile))

        # Check that nothing is returned when no files are found
        count = 0
        for mfile in find_files_matching(self.test_dir,
                                         patterns='nosuchfile.txt'):
            count += 1
        self.assertEqual(count, 0)

    def test_find_file_in_tree(self):
        # Check that find_file_in_tree finds a file that exists.
        found = find_file_in_tree(self.test_file_name, '.')
        self.assertTrue(os.path.isfile(found))
        # It returns an empty string if nothing is found.
        found = find_file_in_tree('nosuchfile.txt', '.')
        self.assertEqual(found, '')
        
    def test_find_file_in_path(self):
        # Check that find_file_in_tree finds a file that exists.
        found = find_file_in_path(self.test_file_name)
        self.assertTrue(os.path.isfile(found))
        # It returns an empty string if nothing is found.
        found = find_file_in_path('nosuchfile.txt')
        self.assertEqual(found, '')
        
    def test_find_file_in_pythonpath(self):
        # Check that find_file_in_pythonpath finds a file that exists.
        # There must always be an '__init__.py' somewhere in PYTHONPATH.
        found = find_file_in_pythonpath('__init__.py')
        self.assertTrue(os.path.isfile(found))
        # It returns an empty string if nothing is found.
        found = find_file_in_pythonpath('nosuchfile.txt')
        self.assertEqual(found, '')
        
    def test_find_file(self):
        # Check that find_file can find various important files.
        found = find_file(self.test_file_name)
        self.assertTrue(os.path.isfile(found))
        found = find_file('test_filesearching.py')
        self.assertTrue(os.path.isfile(found))
        found = find_file('example_properties.py')
        self.assertTrue(os.path.isfile(found))
        
        # The function can either return an empty string or raise an
        # exception if the file is not found.
        found = find_file('nosuchfile.txt', canraise=False)
        self.assertEqual(found, '')
        self.assertRaises(NameError, find_file, 'nosuchfile.txt',
                          canraise=True)
        
    def test_find_file_prefix(self):
        # Check that find_file_prefix can find various files.
        # This time the string returned is not a valid file name,
        # but it should not be the same as the string given.
        found = find_file_prefix(self.test_file_prefix)
        self.assertNotEqual(self.test_file_prefix,found)
        found = find_file_prefix('test_')
        self.assertNotEqual('test_',found)
        found = find_file_prefix('example_')
        self.assertNotEqual('example_',found)
        
        # The function can either return the prefix or raise an
        # exception if the file is not found.
        found = find_file_prefix('nosuchfile', canraise=False)
        self.assertEqual(found, 'nosuchfile')
        self.assertRaises(NameError, find_file_prefix, 'nosuchfile',
                          canraise=True)
        
# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
