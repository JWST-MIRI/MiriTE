#!/usr/bin/env python
#
# :History:
# 
# 13 Feb 2013: Created. Make a copy of the data product before modifying it.
# 15 Feb 2013: Updated the default field-def table.
# 15 Feb 2013: Added modifications to keywords in header of output file
#              to conform with Calibration Data Product delivery 1 
#              specifications. Added import of os. (Vincent Geers, DIAS) 
# 26 Feb 2013: BIAS image removed on advice from Jane Morrison.
# 28 Feb 2013: Fixed range used for extracting the error planes from lindata
#              to be put into errdata (range was off by 1 plane).
# 01 Mar 2013: Updated to match Jane Morrison's latest DHAS flags.
# 03 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
#              - a solution to the "Types in column 0 do not match" problem
#              suggested by Michael Droettboom at STScI.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @authors: Steven Beard (UKATC), Vincent Geers (DIAS)
#

"""

Script `convert_linearity` creates a FITS file in standard MIRI
format from a linearity file in DHAS format. The script relies on the
MIRI linearity data product, which is based on the STScI data
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

import optparse
import os, sys, time

import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_linearity_model import MiriLinearityModel

def load_fm_linearity( filename, order=None, ignore_bias=True,
                          convert_quality=False ):
    """
        
    Read FM non-linearity data in the format written by the FM test software
    and return copies of the primary FITS header and the SCI, ERR, DQ,
    MIN, MAX and FITERR arrays contained.
    
    NOTE: If the input data contains a BIAS plane it is ignored.
    
    :Parameters:
    
    filename: str
        The name of a FITS file in DHAS format from which to extract the
        data.
    order: int, optional
        The order of the polynomial fit. By default, this is obtained
        from the 'ORDER' FITS header keyword.
    ignore_bias, optional
        Set True to indicate the input data contains a BIAS plane
        which must be taken into account and ignored.
    convert_quality: bool, optional
        Set to True to modify the data quality array by ignoring the
        "non-saturated data only" flag. The default is False, which
        keeps the flags intact.
        
    :Returns:
    
    datacollection: tuple of pyfits Header + 7 numpy ndarrays
        (header, scidata, errdata, dqdata, mindata, maxdata, fitdata)
    
    """
    # Read a FITS file in DHAS format. The file contains a primary
    # HDU containing linearity data in various planes.
    try:
        # Read the primary header and primary data array.
        # The data is expected to be 3-D
        hdulist = pyfits.open(filename)
        fitsheader = hdulist[0].header
        lindata = hdulist[0].data
        if lindata.ndim != 3:
            raise TypeError("3-D data expected")
        nplanes = lindata.shape[0]
  
        # Obtain the polynomial order and check the data array has
        # the expected size.
        if order is None:
            if 'ORDER' in fitsheader:
                order = fitsheader['ORDER']
            else:
                if ignore_bias:
                    order = ((nplanes - 5)/2.0) - 1
                else:
                    order = ((nplanes - 4)/2.0) - 1

        if ignore_bias:
            expected = 5 + 2 * (order+1)
        else:
            expected = 4 + 2 * (order+1)
            
        if nplanes != expected:
            strg = "LinearityCorrection - %d planes are " % expected
            strg += "expected for order %d " % order
            strg += "(%d given)." % nplanes
            raise TypeError(strg)

        # Planes 1 and 2 contains the MIN and MAX data.        
        mindata = lindata[0,:,:]
        maxdata = lindata[1,:,:]
                
        if ignore_bias:
            # Skip the BIAS data in pkane 3.
            #biasdata = lindata[2,:,:]
        
            # Plane 4 contains the quality data.
            if convert_quality:
                # The fit is valid when the quality flag is 0 (valid) or
                # 4 (valid but determined from non-saturated data only).
                # TODO: There should be another array preserving the original flags.
                qualityflags = lindata[3,:,:]
                dqdata = qualityflags.copy()
                nonsat = np.where(qualityflags == 4)
                dqdata[nonsat] = 0
            else:
                dqdata = lindata[3,:,:]
        
            # Planes 5 to 5+order+1 contain the correction terms.
            scidata = lindata[4:4+order+1,:,:]
        
            # The next order+1 planes contain the uncertainty in the
            # correction terms.
            errdata = lindata[4+order+1:4+order+order+2,:,:]
        else:
            # Without a bias plane, plane 3 contains the quality data.
            if not convert_quality:
                # The fit is valid when the quality flag is 0 (valid) or
                # 4 (valid but determined from non-saturated data only).
                # TODO: There should be another array preserving the original flags.
                qualityflags = lindata[2,:,:]
                dqdata = qualityflags.copy()
                nonsat = np.where(qualityflags == 4)
                dqdata[nonsat] = 0
            else:
                dqdata = None
        
            # Planes 3 to 4+order+1 contain the correction terms.
            scidata = lindata[3:3+order+1,:,:]
        
            # The next order+1 planes contain the uncertainty in the
            # correction terms.
            errdata = lindata[3+order+1:3+order+order+2,:,:]
        
        # The last plane contains the error on the fit
        fitdata = lindata[-1,:,:]
          
    except Exception as e:
        # If the file could not be read re-raise the exception
        # with a more meaningful error message.
        strg = \
            "%s: Could not read DHAS format non-linearity FITS file.\n   %s" % \
            (e.__class__.__name__, e)
        raise IOError(strg)
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return (fitsheader, scidata, errdata, dqdata, mindata, maxdata, fitdata)

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a DHAS format linearity file into "
    usage += "standard MIRI CDP format."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--order", dest="fitorder", type="int",
                     default=None, help="Order of the fit (if not in the file header)"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot bad pixel mask"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )
    parser.add_option("-i", "--ignore", dest="ignore", action="store_true",
                      help="Ignore warning flags in the quality array"
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

    fitorder = options.fitorder
    verb = options.verb
    makeplot = options.makeplot
    overwrite = options.overwrite
    ignore  = options.ignore
    if ignore is None:
        ignore = False

    # Read the header and data from given file.
    if verb:
        print("Reading %s" % inputfile)
    (header, scidata, errdata, dqdata, mindata, maxdata, fitdata) = \
        load_fm_linearity( inputfile, order=fitorder,
                              convert_quality=ignore )

    # TODO: Could this info be extracted from the header comments or text file?
    # But if the values are always the same it could be defined like this.
    dqdef = [(0,  'bad',          'Bad data'),
             (1,  'badpixel',     'Bad Pixel'),
             (2,  'refpixel',     'Reference Pixel'),
             (3,  'notsaturated', 'Correction from non-saturated subset only'),
             (4,  'higherror',    'High errors'),
             (5, 'negslope',     'Negative slope')
            ]

    # Convert the DHAS product into a MIRI CDP product
    with MiriLinearityModel( init=inputfile ) as oldproduct:
        # Make a copy of the data product just created.
        linproduct = oldproduct.copy()
        linproduct.data = scidata
        linproduct.err = errdata
        linproduct.dq = dqdata
        linproduct.field_def = dqdef
        linproduct.minimage = mindata
        linproduct.maximage = maxdata
        linproduct.fiterr = fitdata
        
        # Apply modifications required to meet CDP-1 specifications
        if 'READOUT' in header and not linproduct.meta.exposure.readpatt:
            linproduct.meta.exposure.readpatt = header['READOUT']
            
        linproduct.meta.filename_original = os.path.basename(inputfile)
        linproduct.meta.filename = os.path.basename(outputfile)
        linproduct.meta.origin = 'MIRI European Consortium'
        
        if verb:
            print(linproduct)
        if makeplot:
            linproduct.plot()

        linproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del oldproduct, linproduct
