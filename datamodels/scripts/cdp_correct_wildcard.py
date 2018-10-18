#!/usr/bin/env python
#
# :History:
# 
# 10 Aug 2018: Created.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_correct_wildcard` corrects wildcard references within CDP metadata.

Prior to the CDP-7 release, MIRI CDPs would set a metadata keywords to 'ANY'
to indicate that the CDP was valid for any variant of that CDP (e.g. FILTER='ANY').
From CDP-7 onwards, the naming convention is changed so the string 'N/A' is
used instead, which is more compatible with the JWST CRDS searching mechanism.

This script checks the keywords contained in a CDP file and changes all
occurrences of 'ANY' to 'N/A'.

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
        Overwrite wildcard existing FITS file.

"""

# Python logging facility.
import logging
# Set the default logging level.
logging.basicConfig(level=logging.INFO)
# Get a default parent logger
logger = logging.getLogger("cdp_correct_wildcard") 

import optparse
import sys, time
 
import miri.datamodels        

def correct_wildcard_metadata( datamodel ):
    """
    
    Correct the wild card used to .
        
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be updated.
        
    :Returns:
    
    nchanges: int
        Returns the number of changes made to the metadata.
    
    """
    # Check MRS wildcard information.
    nchanges = 0
    if hasattr(datamodel, 'meta') and hasattr(datamodel.meta, 'instrument') and \
       hasattr(datamodel.meta, 'exposure') and hasattr(datamodel.meta, 'subarray'):
        
        if datamodel.meta.instrument.model is not None:
            if str(datamodel.meta.instrument.model).strip() == 'ANY':
                datamodel.meta.instrument.model = 'N/A'
                nchanges += 1
        if datamodel.meta.instrument.detector is not None:
            if str(datamodel.meta.instrument.detector).strip() == 'ANY':
                datamodel.meta.instrument.detector = 'N/A'
                nchanges += 1
        if datamodel.meta.instrument.detector_settings is not None:
            if str(datamodel.meta.instrument.detector_settings).strip() == 'ANY':
                datamodel.meta.instrument.detector_settings = 'N/A'
                nchanges += 1
        if datamodel.meta.instrument.filter is not None:
            if str(datamodel.meta.instrument.filter).strip() == 'ANY':
                datamodel.meta.instrument.filter = 'N/A'
                nchanges += 1
        if datamodel.meta.instrument.channel is not None:
            if str(datamodel.meta.instrument.channel).strip() == 'ANY':
                datamodel.meta.instrument.channel = 'N/A'
                nchanges += 1
        if datamodel.meta.instrument.band is not None:
            if str(datamodel.meta.instrument.band).strip() == 'ANY':
                datamodel.meta.instrument.band = 'N/A'
                nchanges += 1
        if datamodel.meta.exposure.readpatt is not None:
            if str(datamodel.meta.exposure.readpatt).strip() == 'ANY':
                datamodel.meta.exposure.readpatt = 'N/A'
                nchanges += 1
        if datamodel.meta.subarray.name is not None:
            if str(datamodel.meta.subarray.name).strip() == 'ANY':
                datamodel.meta.subarray.name = 'N/A'
                nchanges += 1
    else:
        strg = "MIRI instrument, exposure and subarray metadata attributes missing from data model %s" % \
            datamodel.__class__.__name__
        raise TypeError(strg)
    return nchanges

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Corrects the wildcard usage (\'ANY\'-->\'N/A\') within a "
    usage += "MIRI calibration data product."
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
        # Attempt to correct the wildcards in the metadata keyword
        logger.info("Analysing %s..." % inputfile)
        nchanges = correct_wildcard_metadata( datamodel )
                    
        if verb:
            print(datamodel)
            print(datamodel.get_history_str())

        if nchanges > 0:
            datamodel.save( outputfile, overwrite=overwrite)
            logger.info("%d changes made. Data saved to new file, %s\n" % (nchanges, outputfile))
        else:
            logger.info("Data not changed. No output file written.\n")             
           
        del datamodel
