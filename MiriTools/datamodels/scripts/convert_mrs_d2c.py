#!/usr/bin/env python
#
# :History:
# 
# 13 Feb 2013: Created. Make a copy of the data product before modifying it.
# 15 Feb 2013: Added modifications to keywords in header of output file
#              to conform with Calibration Data Product delivery 1 
#              specifications. Added import of os. (Vincent Geers, DIAS)
# 07 Jun 2013: wavelengthdata renamed to wavelength, alphadata renamed to
#              alpha, and slicedata renamed to slicenumber. beta added but
#              commented out for now. jwst_lib problem prevents wavelength
#              array being copied properly.
# 18 Jun 2013: Problem with wavelength array now solved in jwst_lib.
# 31 Oct 2013: BETA array removed from MRS D2C model.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Steven Beard (UKATC)
#

"""

Script `convert_mrs_d2c` creates a FITS file in standard MIRI format
from an MRS D2C file. The script relies on the MIRI MRS D2C data product,
which is based on the STScI data model.

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

from __future__ import division, print_function

import optparse
import os, sys, time

# import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_distortion_models import MiriMrsD2CModel

def load_fm_d2c( filename ):
    """
    
    Reads a header and mask from a MIRI FM format D2C FITS file.
    
    """
    # Read the FITS header and primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        fitsheader = hdulist[0].header
        wavelength = hdulist[0].data[0]
        alpha = hdulist[0].data[1]
        slicenumber = hdulist[0].data[2]
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return fitsheader, wavelength, alpha, slicenumber

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a FM format D2C file into "
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

    # Read the header and data from given file.
    if verb:
        print("Reading %s" % inputfile)
    (header, wavelength, alpha, slicenumber) = load_fm_d2c( inputfile )

    # Convert the DHAS mask into a MIRI CDP mask
    with MiriMrsD2CModel( init=inputfile ) as oldproduct:
        # Make a copy of the data product just created.
        d2cproduct = oldproduct.copy()
        d2cproduct.wavelength = wavelength
        d2cproduct.alpha = alpha
        d2cproduct.slicenumber = slicenumber

        
        # Apply modifications required to meet CDP-1 specifications
        if 'READOUT' in header and not d2cproduct.meta.exposure.readpatt:
            d2cproduct.meta.exposure.readpatt = header['READOUT']
            
        d2cproduct.meta.filename_original = os.path.basename(inputfile)
        d2cproduct.meta.filename = os.path.basename(outputfile)
        d2cproduct.meta.origin = 'MIRI European Consortium'
        
        if verb:
            print(d2cproduct)
        if makeplot:
            d2cproduct.plot()

        d2cproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del oldproduct, d2cproduct
