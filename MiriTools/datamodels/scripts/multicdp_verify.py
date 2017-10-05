#! /usr/bin/env python
#
# :History:
#
# 09 Oct 2013: Created.
# 21 Oct 2013: Count the number of files passing and failing the test.
# 06 Jul 2015: Check for a MemoryError (which caused a blank failure
#              string to be reported.
# 04 Dec 2015: Make the "pass" messages optional, so the failures are
#              easier to spot. Added --nopass option.
# 07 Dec 2015: Added a basic file structure test, to detect really
#              strangely formatted CDP files.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Ruyman Azzollini (DIAS), Steven Beard (UKATC)
#
# -*- coding: utf-8 -*-
"""

This script verifies that FITS files in path, AND containing MIRI CDPs 
are compatible with the current JWST and MIRI software release. 

To pass the test, the file must be capable of being read and written 
without error, and must be unchanged when written to a file and read 
back again.

NOTE: This script creates a temporary copies of the inspected files
with the string "_copy.fits" appended to the name. The script will
fail if a file with that name already exists, unless --overwrite is
specified as a parameter. The temporary files will be removed unless
--keepfile is specified as a parameter.


The command takes the following options::

    --path <type-string>
        The directory in which to look for MIRI CDP files. If not given
        the script will inspect the current directory and all directories
        below it.

    --overwrite or -o
        Overwrite any existing FITS file when making a temporary copy.
        
    --keepfile or -k
        Keep the temporary copy file instead of removing it when finished.
        
    --nopass or -n
        Pass messages are suppressed so that only failures are shown on
        the screen. This helps keep the output tidy if there are a lot
        of files but only a few fail.

"""

from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, time
from pdb import set_trace as stop

# Import the MIRI CDP utilities.
from miri.datamodels.util import verify_fits_file, verify_cdp_file
from miri.tools.filesearching import find_files_matching

def main(path, overwrite=False, keepfile=False, nopass=False):
    # Search for FITS files in the given path.
    # single_level=False means all subfolders will be searched.
    filenames = []
    for filename in find_files_matching( path, patterns='*.fits',
                                         single_level=False):
        filenames.append(filename)
    
    print("%i files to be tested...\n" % len(filenames))
    success_count = failure_count = untested_count = 0
    
    for filename in filenames:
        bfilename = os.path.basename(filename)
        try:
            verify_fits_file(filename, cdp_checks=True)
            time.sleep(0.1) # Allow the file to close.
            datatype = verify_cdp_file( filename, overwrite=overwrite, 
                                        keepfile=keepfile)
            if not nopass:
                print("File \'%s\' of data type \'%s\' has passed the test." % \
                      (bfilename, datatype))
            success_count += 1
        except MemoryError:
            print("*** Insufficient memory to verify file \'%s\'!" % bfilename)
            untested_count += 1
        except Exception as e:
            print("File \'%s\' has failed the test." % bfilename)
            print("  " + str(e))
            failure_count += 1
        
    print("\n%d files passed the test." % success_count)
    if untested_count > 0:
        print("*** WARNING: %d files could not be tested! Scroll up for details. ***" % untested_count)
    if failure_count > 0:
        print("*** WARNING: %d files failed the test! Scroll up for details. ***" % failure_count)
 
if __name__ == "__main__":
    
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt]\n"
    usage += "Verify the contents of a set of FITS files containing "
    usage += "MIRI calibration data products (hanging from some path).\n\n"
    usage += "If no path is provided the current working directory is\n" 
    usage += "searched for fits files to be inspected."
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--path", dest="path", type="string",
                     default=None, help="path in which to search for the CDPs")
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy of copy file if it already exists"
                     )
    parser.add_option("-k", "--keepfile", dest="keepfile", action="store_true",
                      help="Do not remove the copy file when finished"
                     )
    parser.add_option("-n", "--nopass", dest="nopass", action="store_true",
                      help="Suppress pass messages. Only show failures."
                     )
                     
    (options, args) = parser.parse_args()
    
    if options.path != None:
        path = options.path
    else:
        print("Inspecting current directory... hope that was your intention!\n")
        path = os.getcwd()
    
    overwrite = options.overwrite
    keepfile = options.keepfile   
    nopass = options.nopass 
    
    main( path, overwrite=overwrite, keepfile=keepfile, nopass=nopass )
