#!/usr/bin/env python
#
# :History:
# 
# 06 Nov 2015: Created.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_add_subarray` checks the metadata within a MIRI calibration
data product and ensures that all necessary subarray metadata is present. 
The CDP standard changed from making SUBARRAY keywords optional to
making them compulsory, which meant some FULL frame CDPs without SUBARRAY
keywords started to violate the standard. This script can be used to
correct FULL frame CDPs whose only problem is missing SUBARRAY keywords.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options::

    --subarray <subarray-name-string>
        The name of the subarray to be applied. It is only necessary to
        specify a new subarray if 'FULL' is not appropriate, or if the
        file specifies an incorrect subarray mode.
        By default, the subarray mode is inferred from the file content.
    --verbose or -v
        Generate more verbose output.
    --overwrite or -o
        Overwrite any existing FITS file.

"""

from __future__ import absolute_import, unicode_literals, division, print_function

import warnings
import optparse
import sys, time

from miri.datamodels.util import add_subarray_metadata

# Import all the data models that might be contained in the file.
#from miri.datamodels.miri_measured_model import MiriMeasuredModel
import miri.datamodels        

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Adds SUBARRAY records to "
    usage += "any MIRI calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--subarray", dest="subarray", type="string",
                     default="", help="Name of subarray mode to be applied"
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
            outputfile = inputfile + "_out.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    subname = options.subarray
    verb = options.verb
    overwrite = options.overwrite

    # Open the data model using the class derived from the data type.
    with miri.datamodels.open( init=inputfile ) as datamodel:
        subarray_modified = add_subarray_metadata( datamodel, subname=subname )
                    
        if verb:
            print(datamodel)
            print(datamodel.get_history_str())

        if subarray_modified:
            datamodel.save( outputfile, overwrite=overwrite)
            print("Data saved to %s\n" % outputfile)
        else:
            print("Data not changed. No output file written.\n")             
           
        del datamodel
