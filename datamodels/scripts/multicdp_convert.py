#! /usr/bin/env python
#
# :History:
#
# 01 Oct 2014: Created (based on multicdp_verify)
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Ruyman Azzollini (DIAS), Steven Beard (UKATC)
#
# -*- coding: utf-8 -*-
"""

This script converts all files in the current path from CDP-2 format
to CDP-3. 

NOTE: This script writes new versions of the files with new names
of the form <old-name>_new.fits. The script will fail if a file with
that name already exists, unless --overwrite is specified as a parameter.

NOTE: The script applies the given detector settings to all files.
If the path contains files with different detector settings, it is
better to use the individual "cdp_convert" script.

The command takes the following options:

    --path <type-string>
        The directory in which to look for MIRI CDP files. If not given
        the script will inspect the current directory and all directories
        below it.
        
    --settings <detector-settings>
        The detector settings used to create the CDP file.
        Allowed values are 'RAL1', 'JPL1' or 'ANY'        
        
    --overwrite or -o:
        Overwrite any existing FITS file when making a temporary copy.

"""

import optparse
import os
from pdb import set_trace as stop

# Import the MIRI CDP utilities.
from miri.datamodels.util import convert_cdp_2to3
from miri.tools.filesearching import find_files_matching

def main(path, settings=None, overwrite=False):
    # Search for FITS files in the given path.
    # single_level=False means all subfolders will be searched.
    filenames = []
    for filename in find_files_matching( path, patterns='*.fits',
                                         single_level=False):
        filenames.append(filename)
    
    print("%i files to be converted...\n" % len(filenames))
    success_count = failure_count = 0
    
    for filename in filenames:
        bfilename = os.path.basename(filename)
        outputfile = filename + "_new.fits"
        try:
            # TODO: When different conversions are possible, insert if-then-else block here.
            datatype = convert_cdp_2to3(filename, outputfile, settings=settings,
                                        overwrite=overwrite)
            print("File \'%s\' of data type \'%s\'\n  has been converted." % \
                  (bfilename, datatype))
            success_count += 1
        except Exception as e:
            print("File \'%s\' failed to convert." % bfilename)
            print("  " + str(e))
            failure_count += 1
        
    print("\n%d files were successfully converted." % success_count)
    if failure_count > 0:
        print("*** WARNING: %d files failed to convert! Scroll up for details. ***" % failure_count)
 
if __name__ == "__main__":
    
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt]\n"
    usage += "Convert the contents of a set of FITS files containing "
    usage += "MIRI calibration data products (hanging from some path)"
    usage += "from CDP-2 format to CDP-3 format.\n\n"
    usage += "If no path is provided the current working directory is\n" 
    usage += "searched for fits files to be converted."
    parser = optparse.OptionParser(usage)
    parser.add_option("-p", "--path", dest="path", type="string",
                     default=None, help="path in which to search for the CDPs")
    parser.add_option("-s", "--settings", dest="settings", type="string",
                     default=None, help="Detector settings (RAL or JPL)"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy of copy file if it already exists"
                     )
                     
    (options, args) = parser.parse_args()
    
    if options.path != None:
        path = options.path
    else:
        print("Converting all files in current directory... hope that was your intention!\n")
        path = os.getcwd()
    
    settings = options.settings
    if settings == 'RAL':
        settings = 'RAL1'
    elif settings == 'JPL':
        settings = 'JPL1'
    overwrite = options.overwrite
    
    main( path, settings=settings, overwrite=overwrite )
