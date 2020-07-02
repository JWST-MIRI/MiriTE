#!/usr/bin/env python
#
# :History:
# 
# 02 Jul 2013: Created.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Steven Beard (UKATC)
#

"""

Script `convert_slope_data` converts a LVL2 file of slope data in DHAS
format into standard MIRI format. The script relies on the MIRI slope
data product, which is based on the STScI data model.

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out"
        appended.

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

from miri.datamodels.miri_measured_model import MiriSlopeModel

def load_dhas_lvl2( filename ):
    """
    
    Reads a header and data arrays from a DHAS format LVL2 slope data
    FITS file. A DHAS file contains the following HDUs:
    
    * A primary HDU containing the FITS header and a 3 x rows x columns
      array with planes containing the slope, uncertainty and data
      quality averaged over all the integrations.
    * An additional HDU for each integration containing (usually) a 8 x
      rows x columns array with planes containing the slope, uncertainty,
      data quality, zero point, number of good reads, read number of first
      saturated read, number of good segments and fit error.
      The number of planes can vary, but there are always at least 3.
      
    This function extracts the planes from the extension HDUs and combines
    them together into a collection of 3-D arrays of dimension integrations
    x rows x columns.
    
    """
    try:
        # Read the FITS header and average slope data from the FITS file.
        hdulist = pyfits.open(filename)
        header = hdulist[0].header
        avgdata = hdulist[0].data
        
        # The data arrays are expected to have the same first 2 dimensions
        # as the averaged data and 3rd dimension the same as the number
        # of additional HDUs.
        nints = len(hdulist) - 1
        if nints > 0:
            slopeshape = [nints, avgdata.shape[1], avgdata.shape[2]]
            
            signal = np.zeros(slopeshape)
            err = np.zeros(slopeshape)
            dq = np.ones(slopeshape, dtype=np.uint16)
            zeropt = np.zeros(slopeshape)
            nreads = np.zeros(slopeshape, dtype=np.uint16)
            readsat = np.zeros(slopeshape, dtype=np.uint16)
            ngoodseg = np.zeros(slopeshape, dtype=np.uint16)
            fiterr = np.zeros(slopeshape)
            
            intnum = 0
            for hdunum in range(1, len(hdulist)):
                hdu = hdulist[hdunum]
                narrays = hdu.data.shape[0]
                signal[intnum, :, :] = hdu.data[0, :, :]
                err[intnum, :, :] = hdu.data[1, :, :]
                dq[intnum, :, :] = hdu.data[2, :, :]
                if narrays > 3:
                    zeropt[intnum, :, :] = hdu.data[3, :, :]
                    if narrays > 4:
                        nreads[intnum, :, :] = hdu.data[4, :, :]
                        if narrays > 5:
                            readsat[intnum, :, :] = hdu.data[5, :, :]
                            if narrays > 6:
                                ngoodseg[intnum, :, :] = hdu.data[6, :, :]
                                if narrays > 7:
                                    fiterr[intnum, :, :] = hdu.data[7, :, :]      
                intnum += 1
        else:
            signal = avgdata[0]
            err = avgdata[1]
            dq = avgdata[2]
            zeropt = None
            nreads = None
            readsat = None
            ngoodseg = None
            fiterr = None

    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return header, signal, err, dq, zeropt, nreads, readsat, ngoodseg, fiterr

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a DHAS format LVL2 slope data into "
    usage += "standard MIRI slope data format."
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot the slope data"
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
    (header, signal, err, dq, zeropt, nreads, readsat, ngoodseg, fiterr) = \
        load_dhas_lvl2( inputfile )


    # Convert the DHAS mask into a MIRI CDP mask
    with MiriSlopeModel( init=inputfile ) as oldproduct:
        # Make a copy of the data product just created.
        slopeproduct = oldproduct.copy()
        slopeproduct.data = signal
        slopeproduct.err = err
        slopeproduct.dq = dq
        slopeproduct.zeropt = zeropt
        slopeproduct.nreads = nreads
        slopeproduct.readsat = readsat
        slopeproduct.ngoodseg = nreads
        slopeproduct.fiterr = fiterr
        
        # Apply modifications required to meet CDP-1 specifications
        if 'READOUT' in header and not slopeproduct.meta.exposure.readpatt:
            slopeproduct.meta.exposure.readpatt = header['READOUT']
            
        slopeproduct.meta.filename_original = os.path.basename(inputfile)
        slopeproduct.meta.filename = os.path.basename(outputfile)
        slopeproduct.meta.origin = 'MIRI European Consortium'
        
        if verb:
            print(slopeproduct)
        if makeplot:
            slopeproduct.plot()

        slopeproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del oldproduct, slopeproduct
