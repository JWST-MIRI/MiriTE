#!/usr/bin/env python
#
# :History:
#    
# 09 Oct 2013: Created
# 11 Oct 2013: Corrected: missing "keepfile" option. R. Azzollini (DIAS)
# 02 Oct 2014: CDP metadata is now verified.
# 07 Dec 2015: Added a basic file structure test, to detect really
#              strangely formatted CDP files.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Steven Beard (UKATC), R. Azzollini (DIAS)
#

"""

This script verifies that a FITS file containing a MIRI CDP is compatible 
with the current JWST and MIRI software release. To pass the test, the 
file must be capable of being read and written without error, and must be 
unchanged when written to a file and read back again.

NOTE: This script creates a temporary copy of the specified file
with the string "_copy.fits" appended to the name. The script will
fail if a file with that name already exists, unless --overwrite is
specified as a parameter. The temporary file will be removed unless
--keepfile is specified as a parameter.

The following command arguments are defined by position::

    inputfile[0]
    
        The path+name of the file to be read. Compulsory.

The command also takes the following options::

    --datatype <type-string>
        The name of the data type to be used to read the product.
        If specified, this option overrides the TYPE keyword
        contained in the input file.

    --overwrite or -o
        Overwrite any existing FITS file when making a temporary copy.
        
    --keepfile or -k
        Keep the temporary copy file instead of removing it when finished.

"""



import optparse
import sys, time

# Import the MIRI CDP utilities.
from miri.datamodels.cdp import CDP_DICT
from miri.datamodels.util import verify_fits_file, verify_cdp_file
     
if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile\n"
    usage += "Verify the contents of a FITS file containing "
    usage += "any MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--datatype", dest="datatype", type="string",
                     default=None, help="Data type to use (overriding TYPE)"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy file if it already exists"
                     )
    parser.add_option("-k", "--keepfile", dest="keepfile", action="store_true",
                      help="Do not remove the copy file when finished"
                     )

    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    overwrite = options.overwrite
    keepfile = options.keepfile
    
    if options.datatype:
        # Use the data type specified
        datatype = str(options.datatype)
    else:
        datatype = None
    
    # Call the utility function which verifies a CDP file.
    try:
        verify_fits_file(inputfile, cdp_checks=True)
        time.sleep(0.1) # Allow the file to close.
        datatype = verify_cdp_file( inputfile, datatype, overwrite=overwrite, 
                                    keepfile=keepfile )
        print("File \'%s\' of data type \'%s\' has passed the test." % \
              (inputfile, datatype))
    except Exception as e:
        print("File \'%s\' has failed the test." % inputfile)
        print("  " + str(e))
