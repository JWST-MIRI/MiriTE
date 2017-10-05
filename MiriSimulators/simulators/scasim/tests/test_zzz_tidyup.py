#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_zzz_tidyup - Dummy test module designed to run last and
tidy up temporary CDP files left over by previous tests.

NOTE: Only works properly when this tidyup test is run last, which
depends on the nature of the test discovery module.

:History:

04 May 2016: Created.
07 Jun 2016: Only remove files contained within "tests" directory.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

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
        # directory matching the known pattern. To protect against
        # accidental removal of needed CDP files, only remove files
        # from a current directory called "tests".
        self.tempdir = './CDP'
        if os.path.isdir(self.tempdir) and \
            ('tests' in os.path.abspath(self.tempdir)):
            self.tidyup = True
            self.tempfiles = find_files_matching( './CDP',
                                patterns='MIRI_*.fits',
                                single_level=True,
                                sortfiles=False,
                                yield_folders=False,
                                yield_path_only=False)
        else:
            self.tidyup = False
            self.tempfiles = []
        
    def test_dummy(self):
        pass
    
    def tearDown(self):
        # Remove temporary files, if they exist and if able to.
        count = 0
#         time.sleep(1) # Wait a second
        if self.tidyup:
            for tempfile in self.tempfiles:
                if os.path.isfile(tempfile):
                    try:
                        os.remove(tempfile)
                        count += 1
                    except Exception as e:
                        strg = "Could not remove CDP file, " + tempfile + \
                            "\n   " + str(e)
                        warnings.warn(strg)
            os.rmdir(self.tempdir)
        if count > 0:
            print("Successfully tidied up %d CDP files." % count)


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
