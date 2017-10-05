#!/usr/bin/env python
#
# :History:
# 
# 13 Feb 2013: Created. Make a copy of the data product before modifying it.
# 15 Feb 2013: Updated the default field-def table.
# 15 Feb 2013: Added modifications to keywords in header of output file
#              to conform with Calibration Data Product delivery 1 
#              specifications. Added import of os. (Vincent Geers, DIAS) 
# 22 Aug 2013: columnnames renamed to fieldnames.
# 03 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
#              - a solution to the "Types in column 0 do not match" problem
#              suggested by Michael Droettboom at STScI.
# 04 Oct 2013: Standard bad pixel flags obtained from MIRI reserved flags
#              and the table declared with MiriBadPixelMask.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Steven Beard (UKATC)
#

"""

Script `convert_bad_pixel_mask` creates a FITS file in standard MIRI
format from a bad pixel mask in DHAS format. The script relies on the
MIRI bad pixel mask data product, which is based on the STScI data
model.

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

#import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.dqflags import reserved_flags
from miri.datamodels.miri_badpixel_model import \
    bad_pixel_flags, MiriBadPixelMaskModel

def load_dhas_mask( filename ):
    """
    
    Reads a header and mask from a DHAS format bad pixel mask FITS file.
    
    """
    # Read the FITS header and primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        fitsheader = hdulist[0].header
        maskdata = hdulist[0].data
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return fitsheader, maskdata

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a DHAS format bad pixel mask into "
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
    (header,maskdata) = load_dhas_mask( inputfile )
    
    # The DETECTOR keyword is compulsory.
    detector = header['DETECTOR']
    # Convert the DHAS mask into a MIRI CDP mask
    with MiriBadPixelMaskModel( init=inputfile, detector=detector ) as oldproduct:
        # Make a copy of the data product just created.
        maskproduct = oldproduct.copy()
        maskproduct.mask = maskdata
        maskproduct.field_def = reserved_flags + bad_pixel_flags
        
        # Apply modifications required to meet CDP-1 specifications
        if 'READOUT' in header and not maskproduct.meta.exposure.readpatt:
            maskproduct.meta.exposure.readpatt = header['READOUT']
            
        maskproduct.meta.filename_original = os.path.basename(inputfile)
        maskproduct.meta.filename = os.path.basename(outputfile)
        maskproduct.meta.origin = 'MIRI European Consortium'
        
        if verb:
            print(maskproduct)
            print(maskproduct.flags_table)
        if makeplot:
            maskproduct.plot()

        maskproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del oldproduct, maskproduct
