#!/usr/bin/env python
#
# :History:
# 
# 10 Aug 2018: Created.
#
# @author: Steven Beard (UKATC)
#
"""

Script `cdp_correct_band` checks the MRS BAND keyword within a MIRI calibration
data product and ensures the value is in the correct format.

The script uses the following algorithm:

    Open the MIRI MRS CDP file.
    If the BAND keyword contains a valid value then
        Nothing needs to be done to the data model.
    Else if a --band parameter has been provided
        Use the value contained in the --band parameter to correct
        the BAND keyword.
    Else if the data model contains valid DGAA and DGAB keywords
        Use those keywords to correct the BAND keyword.
    Else
        Attempt to infer the BAND value from the input file name.
    Save a copy of the data model to a new file.

The following command arguments are defined by position::

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options::

    --band <band-name-string>
        The name of the MRS band to be applied. It is only necessary to
        specify a band name if it cannot be inferred from the DGAA and DGAB
        keywords.
    --verbose or -v
        Generate more verbose output.
    --overwrite or -o
        Overwrite any existing FITS file.

"""

# Python logging facility.
import logging
# Set the default logging level.
logging.basicConfig(level=logging.INFO)
# Get a default parent logger
logger = logging.getLogger("cdp_correct_band") 

import optparse
import sys, time
 
from miri.parameters import MIRI_BANDS_SINGLE, MIRI_BANDS_CROSS, MIRI_BANDS
import miri.datamodels        

def correct_band_metadata(datamodel, band='', filename=''):
    """
    
    Correct the MRS band in the metadata of a data model.
        
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be updated.
    band: str (optional)
        The name of the MIRI MRS band to be applied. It is only necessary to
        specify a new band if the data model does not contain valid
        DGAA and DGAB keywords. By default, the band is inferred from the
        DGAA and DGAB keywords.
        
    :Returns:
    
    band_modified: bool
        Returns True if the BAND metadata has been modified.
    
    """
    # Check MRS band information.
    band_modified = False
    if hasattr(datamodel, 'meta') and hasattr(datamodel.meta, 'instrument'):
        
        if hasattr(datamodel.meta.instrument, 'detector'):
            if str(datamodel.meta.instrument.detector) not in ['MIRIFUSHORT', 'MIRIFULONG']:
                logger.error("Not a MIRI MRS data model. Nothing to correct.")
                return False   
        else:
            strg = "MIRI detector metadata attributes missing from data model %s" % \
                datamodel.__class__.__name_
            raise TypeError(strg)

        if hasattr(datamodel.meta.instrument, 'band'):
        
            band_acceptable_values = MIRI_BANDS + ['ANY', 'N/A']
            if datamodel.meta.instrument.band is not None and \
               str(datamodel.meta.instrument.band) in band_acceptable_values:
                # The BAND keyword is already correct. Nothing needs to be done.
                logger.info("BAND keyword already contains \'%s\'. No correction needed." % \
                      str(datamodel.meta.instrument.band))
            else:
                # First attempt to use the BAND parameter provided
                if band is not None and band in MIRI_BANDS:
                    datamodel.meta.instrument.band = band
                    logger.info("BAND keyword corrected to \'%s\' from --band parameter." % \
                          band)
                    band_modified = True
                else:
                    # Invalid or no parameter given. Look for DGAA and DGAB metadata.
                    if hasattr(datamodel.meta.instrument, "dichroic_a") and \
                       hasattr(datamodel.meta.instrument, "dichroic_b") and \
                       datamodel.meta.instrument.dichroic_a is not None and \
                       datamodel.meta.instrument.dichroic_b is not None and \
                       str(datamodel.meta.instrument.dichroic_a) in MIRI_BANDS_SINGLE and \
                       str(datamodel.meta.instrument.dichroic_b) in MIRI_BANDS_SINGLE:
                        # Valid DGAA and DGAB parameters. Construct a BAND keyword from them.
                         dgaa = str(datamodel.meta.instrument.dichroic_a)
                         dgab = str(datamodel.meta.instrument.dichroic_b)
                         band = "%s-%s" % (dgaa, dgab)
                         datamodel.meta.instrument.band = band
                         logger.info("BAND keyword corrected to \'%s\' from DGAA and DGAB." % \
                               band)
                         band_modified = True
                    else:
                        # No DGAA/DGABB keywords and no band suggested.
                        # Attempt to derive a correct BAND name from the filename.
                        # The band name comes immediately after a channel name
                        # string ('12' or '34').
                        if filename:
                            matched_band = ''
                            for test_band in MIRI_BANDS_SINGLE:
                                test_string = '12' + test_band
                                if test_string in filename:
                                    matched_band = test_band
                                test_string = '34' + test_band
                                if test_string in filename:
                                    matched_band = test_band
                            # Test for cross-dichroic strings last so they
                            # override a false matche with a single mode.
                            for test_band in MIRI_BANDS_CROSS:
                                test_string = '12' + test_band
                                if test_string in filename:
                                    matched_band = test_band
                                test_string = '34' + test_band
                                if test_string in filename:
                                    matched_band = test_band
                            if matched_band:
                                 datamodel.meta.instrument.band = matched_band
                                 logger.warning("BAND keyword corrected to \'%s\' from the filename only." % \
                                       matched_band)
                                 band_modified = True                           
                            else:
                                # Run out of ideas
                                strg = "Insufficient information to correct the BAND keyword."
                                strg = " Please provide a --band parameter and try again."
                                logger.error(strg)
                        else:
                            # Run out of ideas
                            strg = "Insufficient information to correct the BAND keyword."
                            strg = " Please provide a --band parameter and try again."
                            logger.error(strg)
        else:
            strg = "MIRI MRS band metadata attributes missing from data model %s" % \
                datamodel.__class__.__name__
            raise TypeError(strg)
    else:
        strg = "MIRI instrument metadata attributes missing from data model %s" % \
            datamodel.__class__.__name__
        raise TypeError(strg)
    return band_modified

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Corrects the BAND keyword within a "
    usage += "MIRI MRS calibration data product."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--band", dest="band", type="string",
                     default="", help="Name of sub-band to be applied"
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

    band = options.band
    verb = options.verb
    overwrite = options.overwrite

    # Open the data model using the class derived from the data type.
    with miri.datamodels.open( init=inputfile ) as datamodel:
        # Attempt to correct the BAND keyword
        logger.info("Analysing %s..." % inputfile)
        band_modified = correct_band_metadata( datamodel, band=band, filename=inputfile )
                    
        if verb:
            print(datamodel)
            print(datamodel.get_history_str())

        if band_modified:
            datamodel.save( outputfile, overwrite=overwrite)
            logger.info("Data saved to new file, %s\n" % outputfile)
        else:
            logger.info("Data not changed. No output file written.\n")             
           
        del datamodel
