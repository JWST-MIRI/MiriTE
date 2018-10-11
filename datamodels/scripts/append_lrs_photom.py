#!/usr/bin/env python
#
# :History:
# 
# 11 Oct 2018: Created
#
# @author: Steven Beard (UKATC)
#

"""

Script `append_lrs_photom` creates a merged photometric CDP for the MIRI
imager from 3 individual files:

   * An imager-only photometric CDP, such as MIRI_FM_MIRIMAGE_PHOTOM_05.02.00.fits
   * An LRS flux conversion CDP, such as MIRI_FM_MIRIMAGE_P750L_PHOTOM_05.02.00
   * An LRS flux conversion CDP for SLITLESSPRISM, such as MIRI_FM_MIRIMAGE_P750L_SLITLESSPRISM_PHOTOM_05.02.00

The script relies on the MIRI imager photometric data models, which are
based on the STScI JWST data models.

The following command arguments are defined by position:

    imagerfile[0]
        The path+name of the file containing imager photometric data.
        Compulsory.
    lrsfullfile[1]
        The path+name of the file containing LRS FULL flux conversion data.
        Compulsory.
    lrsslitnessfile[2]
        The path+name of the file containing LRS SLITLESSPRISM flux conversion
        data. Compulsory.
    outputfile[3]
        The path+name of the file to be written.
        Optional. Defaults to the same name as imagerfile with "_out" appended.

The command also takes the following options:

    --version
        Override the version number in the file and provide a new one.
    --verbose or -v:
        Generate more output.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""

import optparse
import os, sys, time
import warnings

import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_photometric_models import MiriPhotometricModel
from miri.datamodels.miri_fluxconversion_models import MiriLrsFluxconversionModel

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] imagerfile lrsfullfile lrsslitnessfile [outputfile]\n"
    usage += "Merge an imager photometic CDP and 2 LRS flux conversion CDPs "
    usage += "into one combined CDP."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--version", dest="version", type="string",
                      default=None, help="New version number (overriding file contents)"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )

    (options, args) = parser.parse_args()

    try:
        imager_file = args[0]
        lrsfull_file = args[1]
        lrsslitless_file = args[2]
        if len(args) > 3:
            outputfile = args[3]
        else:
            outputfile = imager_file + "_out.fits"
    except IndexError:
        print(help_text)
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    version = options.version
    verb = options.verb
    overwrite = options.overwrite

    # Read the header and data from given file.
    if verb:
        strg = "Appending %s and %s\n\tto %s." % (lrsfull_file, lrsslitless_file, imager_file)
        print(strg)
        
    # Read the 3 input data models.
    imager_model = MiriPhotometricModel( imager_file )
    lrsfull_model = MiriLrsFluxconversionModel( lrsfull_file )
    lrsslitless_model = MiriLrsFluxconversionModel( lrsslitless_file )
    
    # Append the two LRS models to the imager model.
    imager_model.append_lrs( [lrsfull_model, lrsslitless_model])
    
    # Update the history metadata
    merge_history = "Merged with %s" % lrsfull_file
    imager_model.add_history(merge_history)
    merge_history = "Merged with %s" % lrsslitless_file
    imager_model.add_history(merge_history)
    if lrsfull_model.meta.author:
        merge_history = "LRS data has AUTHOR=%s" % lrsfull_model.meta.author
        imager_model.add_history(merge_history)
    if lrsfull_model.meta.description:
        merge_history = "LRS data has DESCRIPTION=%s" % lrsfull_model.meta.description
        imager_model.add_history(merge_history)
    
    # If necessary, update the version metadata
    if version:
        imager_model.meta.version = version
        
    if verb:
        print(imager_model)
         
    # Save the updated data model to the output file.
    imager_model.save( outputfile, overwrite=overwrite)
    if verb:
        print("Combined data model saved to %s\n" % outputfile)
    
