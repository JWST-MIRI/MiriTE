#!/usr/bin/env python
#
# :History:
# 
# 11 Oct 2018: Created
# 18 Oct 2018: Modified to merge CDP-5 format CDPs together to make a new
#              CDP-7 format CDP. Added --document and --description parameters.
#
# @author: Steven Beard (UKATC)
#

"""

Script `append_lrs_photom` creates a merged photometric CDP for the MIRI
imager from 3 individual files:

   * An imager-only photometric CDP (in CDP-5 format), such as
     MIRI_FM_MIRIMAGE_PHOTOM_05.02.00.fits
   * An LRS flux conversion CDP, such as
     MIRI_FM_MIRIMAGE_P750L_PHOTOM_05.02.00.fits
   * An LRS flux conversion CDP for SLITLESSPRISM, such as
     MIRI_FM_MIRIMAGE_P750L_SLITLESSPRISM_PHOTOM_05.02.00.fits

The merged CDP is written in CDP-7 format.

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

    --document
        Override the version number in the file and provide a new one.
    --description
        Override the version number in the file and provide a new one.
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

from miri.datamodels.miri_photometric_models import MiriPhotometricModel, \
    MiriPhotometricModel_CDP5, MiriImagingPhotometricModel
from miri.datamodels.miri_fluxconversion_models import MiriLrsFluxconversionModel

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] imagerfile lrsfullfile lrsslitnessfile [outputfile]\n"
    usage += "Merge an imager photometic CDP and 2 LRS flux conversion CDPs "
    usage += "into one combined CDP."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--document", dest="document", type="string",
                      default=None, help="New document name"
                     )
    parser.add_option("", "--description", dest="description", type="string",
                      default=None, help="New description"
                     )
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

    document = options.document
    description = options.description
    version = options.version
    verb = options.verb
    overwrite = options.overwrite

    # Read the header and data from given file.
    if verb:
        strg = "Appending %s and %s\n\tto %s." % (lrsfull_file, lrsslitless_file, imager_file)
        print(strg)
        
    # Read the 3 input data models.
    imager_model = MiriPhotometricModel_CDP5( imager_file )
    lrsfull_model = MiriLrsFluxconversionModel( lrsfull_file )
    lrsslitless_model = MiriLrsFluxconversionModel( lrsslitless_file )
    
    # Create a new imager data model based on the first model.
    new_model = MiriImagingPhotometricModel(phot_table=imager_model.phot_table,
                                     pixar_sr=imager_model.meta.photometry.pixelarea_steradians,
                                     pixar_a2=imager_model.meta.photometry.pixelarea_arcsecsq)
    # Ensure all the metadata is copied (apart from obsolete comments)
    new_model.copy_metadata( imager_model, ignore=['COMMENT'])
    new_model.meta.instrument.detector_settings = 'N/A' # Change 'ANY' to 'N/A'
    new_model.meta.exposure.readpatt = 'N/A' # Change 'ANY' to 'N/A'
    new_model.meta.subarray.name = imager_model.meta.subarray.name
    new_model.meta.subarray.xstart = imager_model.meta.subarray.xstart
    new_model.meta.subarray.ystart = imager_model.meta.subarray.ystart
    new_model.meta.subarray.xsize = imager_model.meta.subarray.xsize
    new_model.meta.subarray.ysize = imager_model.meta.subarray.ysize   
    
    # Append the two LRS models to the new imager model.
    new_model.append_lrs( lrsfull_model, lrsslitless_model )
    
    # Update the description and document references, if requested.
    # Also add some items missing from the history.
    if description is not None and description:
        new_model.meta.description = description
    new_model.add_referencefile_history(document=document,
                                software='IDL scripts',
                                dataused='See table within document',
                                differences='Merged with LRS data')
    
    # Add extra information to the history metadata
    merge_history = "Merged with %s" % lrsfull_file
    new_model.add_history(merge_history)
    merge_history = "Merged with %s" % lrsslitless_file
    new_model.add_history(merge_history)
    if lrsfull_model.meta.author:
        merge_history = "LRS data has AUTHOR=%s" % lrsfull_model.meta.author
        new_model.add_history(merge_history)
    if lrsfull_model.meta.description:
        merge_history = "LRS data has DESCRIPTION=%s" % lrsfull_model.meta.description
        new_model.add_history(merge_history)
    
    # If necessary, update the version metadata
    if version:
        new_model.meta.version = version
    # Display the contents of the data model in verbose mode.
    if verb:
        print(new_model)
         
    # Save the updated data model to the output file.
    new_model.save( outputfile, overwrite=overwrite)
    if verb:
        print("Combined data model saved to %s\n" % outputfile)
