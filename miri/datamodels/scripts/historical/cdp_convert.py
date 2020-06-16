#!/usr/bin/env python
#
# NOTE: Although the conversion carried out by this script is obsolete
# (CDP-2 to CDP-3) it could be reused for different conversions in the future.
#
# :History:
# 
# 30 Sep 2014: Created.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_convert` converts a calibration data product from one format
to another.

CDP-2 FORMAT DATA FILES ARE CONVERTED TO CDP-3 FORMAT.
NOTE: VERIFY EACH OUTPUT FILE AFTER USE.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read. Compulsory.
        
    outputfile[1]
        The path+name of the converted file to be written.
        Optional. If not given, defaults to <oldfile>_new.fits

The command also takes the following options::

    --datatype <type-string>
        The name of the data type to be used to read the product.
        If specified, this option overrides the TYPE keyword
        contained in the input file.

    --settings: <type-string>
        The detector settings used to create the CDP file.
        Allowed values are 'RAL1', 'JPL1' or 'ANY'
        
    --verbose or -v
        Print the converted model before saving it.

    --overwrite or -o
        Overwrite any existing FITS file.

"""

import optparse
import sys, time
import warnings

# Import the MIRI CDP utilities.
from miri.datamodels.util import convert_cdp_2to3

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile [outputfile]\n"
    usage += "Lists the contents of a FITS file compatible with "
    usage += "any MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--datatype", dest="datatype", type="string",
                     default=None, help="Data type to use (overriding TYPE)"
                     )
    parser.add_option("", "--conversion", dest="conversion", type="string",
                     default=None, help="Conversion to be made"
                     )
    parser.add_option("", "--settings", dest="settings", type="string",
                     default=None, help="Detector settings (RAL1, JPL1 or ANY)"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the copy of the file if it already exists"
                     )

    (options, args) = parser.parse_args()

    try:
        inputfile = args[0]
        if len(args) > 1:
            outputfile = args[1]
        else:
            outputfile = inputfile + "_new.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    # Get the detector settings. RAL and JPL are shorthand for RAL1 or JPL1.
    settings = options.settings
    if settings == 'RAL':
        settings = 'RAL1'
    elif settings == 'JPL':
        settings = 'JPL1'
    verb = options.verb
    overwrite = options.overwrite

    if options.datatype:
        # Use the data type specified
        datatype = str(options.datatype)
        print("Forcing the data model to be opened with type \'%s\'" % datatype)
    else:
        datatype = ''

    # Call the utility function which converts a CDP file.
    try:
        # TODO: When different conversions are possible, insert if-then-else block here.
        datatype = convert_cdp_2to3(inputfile, outputfile, datatype=datatype,
                                    settings=settings, printmodel=verb,
                                    overwrite=overwrite)
    except Exception as e:
        print("File \'%s\' failed to convert." % inputfile)
        print("  " + str(e))
    else:
        print("File \'%s\' of data type \'%s\'\n  has been converted to \'%s\'." % \
              (inputfile, datatype, outputfile))
        
