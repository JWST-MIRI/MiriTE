#!/usr/bin/env python
#
# :History:
# 
# 06 Sep 2018: Created.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_add_filter_band` checks the metadata within a MIRI calibration
data product and ensures that all necessary filter, channel and band metadata
are present. 

The JWST pipeline requires that CDP which have different variations for both
detector and wavelength range, and exist in both imager and MRS variations,
must include both the FILTER keyword relevant to the imager and the CHANNEL
and BAND keywords relevant to the MRS. Keywords which are not relevant can
be set to 'N/A'. This allows the CRDS to check these keywords without errors.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options::

    --verbose or -v
        Generate more verbose output.
    --overwrite or -o
        Overwrite any existing FITS file.

"""

import warnings
import optparse
import sys, time

# Import all the data models that might be contained in the file.
import miri.datamodels

def add_missing_filter_band( datamodel ):
    """
    
    Ensure the given data model contains FILTER, CHANNEL and BAND keywords.
        
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be updated.
        
    :Returns:
    
    nreplaced: int
        The number of metadata keywords replaced (0, 1, 2 or 3).
    
    """
    nreplaced = 0
    if hasattr(datamodel, 'meta') and hasattr(datamodel.meta, 'instrument'):
        
        # Replace any missing filter, channel or band attribute with 'N/A'
        if datamodel.meta.instrument.filter is None:
            datamodel.meta.instrument.filter = 'N/A'
            nreplaced += 1
        if datamodel.meta.instrument.channel is None:
            datamodel.meta.instrument.channel = 'N/A'
            nreplaced += 1
        if datamodel.meta.instrument.band is None:
            datamodel.meta.instrument.band = 'N/A'
            nreplaced += 1
    else:
        strg = "MIRI instrument metadata attributes missing from data model %s" % \
            datamodel.__class__.__name_
        raise TypeError(strg)
    return nreplaced

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Adds missing FILTER/CHANNEL/BAND records to any MIRI "
    usage += "calibration data product which requires all three."
    parser = optparse.OptionParser(usage)
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

    verb = options.verb
    overwrite = options.overwrite

    # Open the data model using the class derived from the data type.
    with miri.datamodels.open( init=inputfile ) as datamodel:
        nreplaced = add_missing_filter_band( datamodel )
        if nreplaced > 2:
            print("WARNING: All 3 metadata keywords have been set to \'N/A\'!")
        elif nreplaced > 0:
            print("%d metadata keywords set to \'N/A\'" % nreplaced)
        else:
            print("No changes made.")
                    
        if verb:
            print(datamodel)
            print(datamodel.get_history_str())

        if nreplaced > 0:
            datamodel.save( outputfile, overwrite=overwrite)
            print("Data saved to %s\n" % outputfile)
        else:
            print("Data not changed. No output file written.\n")             
           
        del datamodel
