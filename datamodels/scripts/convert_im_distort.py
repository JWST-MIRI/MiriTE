#!/usr/bin/env python
#
# Script 'convert_im_distort'
#
# :History:
# 
# 19 Feb 2013: Created
# 20 Feb 2013: Filename changed to convert_im_distort.py.
# 21 Feb 2013: Changed to make a copy of the data product (initialized
#              with the inputfile) before modifying it.
# 16 Sep 2013: Removed the ORDER parameter, since it is derivable from
#              the size of the matrices.
# 31 Oct 2013: Corrected import.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS)
#

"""

Script `convert_im_distort` creates a FITS file in standard MIRI format 
from input ASCII table. The script relies on the MiriDistortionModel data 
product, which is based on the STScI data model.

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options:

    --verbose or -v:
        Generate more output.
    --plot or -p:
        Plot the bad pixel mask.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""

import optparse
import os, sys, time

# import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_distortion_models import \
    MiriImagingDistortionModel

def load_distortion_fits( filename ):
    """
    
    Reads distortion parameters from FITS file.
    
    """
    # Read the FITS header and primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        fitsheader = hdulist[0].header
        fitsdata = hdulist[0].data
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return fitsheader, fitsdata

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a distortion FITS file into "
    usage += "standard MIRI CDP format."
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot bad pixel mask"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
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
    makeplot = options.makeplot
    overwrite = options.overwrite

    # Read the given file.
    if verb:
        print("Reading %s" % inputfile)
    header, inputdata,  = load_distortion_fits( inputfile )
    
    # Convert the distortion data into a MIRI Distortion CDP
    with MiriImagingDistortionModel( init=inputfile ) as oldproduct:
        
        # Make a copy of the data product just created.
        distortproduct = oldproduct.copy()
        
        # Copy over distortion data
        distortproduct.cmatrix = inputdata[0]
        distortproduct.rmatrix = inputdata[1]
        
        # default modifications for CDP-1 specifications
        distortproduct.meta.filename_original = os.path.basename(inputfile)
        distortproduct.meta.filename = os.path.basename(outputfile)

        if verb:
            print(distortproduct)
        if makeplot:
            distortproduct.plot()

        distortproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del distortproduct, oldproduct
