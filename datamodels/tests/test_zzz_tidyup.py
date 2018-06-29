#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_zzz_tidyup - Dummy test module designed to run last and
tidy up temporary test files left over by previous tests. Files might
be left over if they can't be removed because of a locking problem.

NOTE: The work-around only works properly when this tidyup test is
run last, which depends on the nature of the test discovery module.

ONLY NEEDED AS LONG AS THE DATA MODEL/PYFITS FILE LOCKING PROBLEM REMAINS.

:History:

24 Jan 2013: Created.
07 Oct 2015: Made exception catching Python 3 compatible.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
# import time
import unittest
import warnings

from miri.tools.filesearching import find_files_matching

class TestTidyUp(unittest.TestCase):
    
    # A dummy test suite designed just to remove the temporary files
    # left over from earlier tests.
    
    def setUp(self):        
        # Search for all the left over test files in the current
        # directory matching the known pattern.
        self.tempfiles = find_files_matching( '.', patterns='Miri*_test.fits',
                                              single_level=True,
                                              sortfiles=False,
                                              yield_folders=False,
                                              yield_path_only=False)
        
    def test_dummy(self):
        pass
    
    def tearDown(self):
        # Remove temporary files, if they exist and if able to.
        count = 0
#         time.sleep(1) # Wait a second
        for tempfile in self.tempfiles: 
            if os.path.isfile(tempfile):
                try:
                    os.remove(tempfile)
                    count += 1
                except Exception as e:
                    strg = "Could not remove temporary file, " + tempfile + \
                        "\n   " + str(e)
                    warnings.warn(strg)
        del self.tempfiles
        if count > 0:
            print("Successfully tidied up %d temporary files." % count)


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
