#!/usr/bin/env python
#
# :History:
# 
# 22 Feb 2013: Created.
# 27 Feb 2013: Added inputtype as an input parameter, to handle 
#              PSF of IM, LRS, MRS differently.
# 28 Feb 2013: IM PSF input data now assumed to have empty primary HDU,
#              SCI data in first extension, and a PSF Look-up Table 
#              (psf_lut) in the second extension. IM PSF files are saved
#              with this psf_lut, while LRS and MRS PSF are not.
# 22 Aug 2013: columnnames renamed to fieldnames.
# 03 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
#              - a solution to the "Types in column 0 do not match" problem
#              suggested by Michael Droettboom at STScI.
# 31 Oct 2013: Corrected imports.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS)
#

"""

Script `convert_psf` creates a FITS file in standard MIRI format. The 
script relies on the MIRI PointSpreadFunction model product, which is
based on the STScI data model.

The following command arguments are defined by position:

    inputtype[0]
        The type of PSF file. Valid values are IMPSF, LRSPSF, MRSPSF.
    inputfile[1]
        The path+name of the file to be read. Compulsory.
    outputfile[2]
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

import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_psf_models import MiriPointSpreadFunctionModel
from miri.datamodels.miri_psf_models import MiriImagingPointSpreadFunctionModel

def load_psf( filename, inputtype ):
    """
    
    Reads a header and mask from a PSF FITS file.
    
    """
    # Read the FITS header and primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        if inputtype == 'IMPSF':
            fitsheader = hdulist[0].header
            fitsdata = hdulist[1].data
            fitslut = hdulist[2].data
            if len(fitslut[0][0]) == 1:
                lut_im = [(fitslut[0][0][0],fitslut[0][1],fitslut[0][2],fitslut[0][3])]
            else:
                lut_im = list(zip(fitslut[0][0],fitslut[0][1],fitslut[0][2],fitslut[0][3]))
            psf_lut = lut_im
#             psf_lut = np.rec.fromrecords(lut_im, dtype=None,
#                                  names=MiriImagingPointSpreadFunctionModel.fieldnames,
#                                  titles=None)
        if inputtype == 'LRSPSF':
            fitsheader = hdulist[0].header
            fitsdata = hdulist[0].data        
            psf_lut = None      
        if inputtype == 'MRSPSF':
            fitsheader = hdulist[0].header
            fitsdata = hdulist[0].data  
            psf_lut = None      
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return fitsheader, fitsdata, psf_lut

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputtype inputfile outputfile\n"
    usage += "Converts a PSF file into "
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
        inputtype = args[0]
        inputfile = args[1]
        if len(args) > 2:
            outputfile = args[2]
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
    header, fitsdata, psf_lut = load_psf( inputfile, inputtype )

    # Convert the PSF into a MIRI CDP
    if inputtype == 'IMPSF':
        oldproduct = MiriImagingPointSpreadFunctionModel( init=inputfile, psf_lut=psf_lut )
    else:
        oldproduct = MiriPointSpreadFunctionModel( init=inputfile )

    # Make a copy of the data product just created.
    psfproduct = oldproduct.copy()
    
    # Copy over the PSF data
    psfproduct.data = fitsdata
        
    # Apply modifications required to meet CDP-1 specifications
    if 'READOUT' in header and not psfproduct.meta.exposure.readpatt:
        psfproduct.meta.exposure.readpatt = header['READOUT']
    psfproduct.meta.filename_original = os.path.basename(inputfile)
    psfproduct.meta.filename = os.path.basename(outputfile)
    psfproduct.meta.origin = 'MIRI European Consortium'

    # Copy the WCS info
    # FIXME: DOESN'T WORK - REJECTS ANY COMMENT OR HISTORY KEYWORD!
#    wcs = psfproduct.get_fits_wcs( 'PRIMARY' )
#    psfproduct.set_fits_wcs( wcs, 'SCI')
            
    if verb:
        print(psfproduct)
    if makeplot:
        psfproduct.plot()

    psfproduct.save( outputfile, overwrite=overwrite)
    if verb:
        print("Data saved to %s\n" % outputfile)
        
    del oldproduct, psfproduct
