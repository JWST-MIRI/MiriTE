#!/usr/bin/env python

"""

miri_run_tests: helper script to run all the unit tests contained within
the miri package.

:History:

12 Mar 2012: Created using Julian Morin's old test discovery system.
19 May 2017: Completely new version based on nosetests.

@author: Steven Beard (UKATC)

"""
from __future__ import absolute_import, unicode_literals, division, print_function

import os, sys
import warnings
import nose

# If being run as a main program, run the tests.
if __name__ == '__main__':
    # Run the unit tests for each MIRI package.
#    print("Testing MiriTools/datamodels...")
#    os.chdir("MiriTools/datamodels")

    params = [ \
        '--exe',
        '--with-xcoverage',
        '--with-xunit',
        '--cover-package=.',
        '--cover-erase']

    testfolders = [ \
        'MiriTools/tools/tests/',
        'MiriTools/datamodels/tests/',
        'MiriSimulators/simulators/tests/',
        'MiriSimulators/simulators/scasim/tests/']
    
    argv = params + testfolders

    # Works
    result = nose.run(  )
    # FIXME: Fails with "OSError: No such file" errors and "Module miri was never imported."
    # FIXME: Removing the "where" statements from setup.cfg changes error to "'module' object has no attribute 'datamodels'"
    # result = nose.run( argv=argv )

    if result:
       print("Tests completed successfully.")
    else:
       print("Tests failed!")
