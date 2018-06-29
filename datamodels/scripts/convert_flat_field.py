#!/usr/bin/env python
#
# Script 'convert_flat_field'
#
# :History:
# 
# 14 Feb 2013: Created
# 18 Feb 2013: Updated to accept inputfiles that are either single array
#              (assumed to be flat data), or a cube, where first plane is
#              flat data, second plane is error data, and optional third 
#              plane is the data quality array.
# 19 Feb 2013: Changed input parameter flattype to inputtype. The inputtype
#              determines the flattype, how to copy over flatdata (plane or
#              cube), err data and data quality (DQ) array, and how to 
#              define field_def.
#
# 11 Jul 2013: Modified by R. Azzollini (DIAS) to accept CDPs that already 
#              have ERR, DQ and FIELD_DEF planes (LRS Pixel Flat-Field and 
#              MRS Sky Flat-Fields).
#
# 22 Jul 2013: Indentation problems corrected.
# 11 Sep 2013: Modified by V. Geers to expect MRS Pixel Flats to already 
#              have SCI, ERR, DQ, and FIELD_DEF extensions.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS), R. Azzollini (DIAS)
#

"""

Script `convert_flat_field` creates a FITS file in standard MIRI format 
from an FITS file. The script relies on the MIRI flatfield model data 
product, which is based on the STScI data model.

The following command arguments are defined by position:

    inputtype[0]
        The type of flat field. Valid values are IMPIXFLAT, LRSPIXFLAT, MRSPIX,
        MRSFRINGEFLAT, MRSSKYFLAT. This parameter determines what shape is 
        assumed for input flat data, and how err, dq, and field_def are defined.
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
import astropy.io.fits as pyfits
from miri.datamodels.miri_flatfield_model import MiriFlatfieldModel

def mainx(inputfile,verb,makeplot,overwrite):

    # Read the given file.
    if verb:
        print("Reading %s" % inputfile)
    
    hdulist = pyfits.open(inputfile)
    
    header = hdulist[0].header
    
    # By default, assume no error plane, dq plane, and field_def available,
    # but use them if available
    
    if len(hdulist) == 1:
        errdata = None
        dqdata = None
        field_def = None
    
    # Based on input type, set flatfield type, field definition, data 
    # plane/cube for flat data, error, and data quality.
    if inputtype == 'IMPIXFLAT':
        # input data expected to be a cube, first plane containing
        # the SCI array, 2nd plane containing the ERR array
        flattype = 'PIXELFLAT'
        inputdata = hdulist[0].data
        flatdata = inputdata[0]
        errdata = inputdata[1]
    elif inputtype == 'LRSPIXFLAT':
        
        # input data expected to be a 5-extensions fits
        flattype = 'PIXELFLAT'
        header = hdulist[0].header
        flatdata = hdulist[1].data
        errdata = hdulist[2].data
        dqdata = hdulist[3].data
        field_def = hdulist[4].data

    elif inputtype == 'MRSPIXFLAT':
        # input data is expected to be correctly formatted, with 
        # empty primary extension, and SCI, ERR, DQ, and FIELD_DEF in
        # in subsequent extensions.
        flattype = 'PIXELFLAT'
        flatdata = hdulist[1].data
        errdata = hdulist[2].data
        dqdata = hdulist[3].data
        field_def = hdulist[4].data

    elif inputtype == 'MRSFRINGEFLAT':
        # input data expected to be a cube, first plane containing
        # the SCI array, 2nd plane containing the ERR array, 3rd plane
        # containing the DQ array. Creating a missing FIELD_DEF here.
        flattype = 'FRINGEFLAT'
        field_def = [(0,  'good',    'Good data')]
        inputdata = hdulist[0].data
        flatdata = inputdata[0]
        errdata = inputdata[1]
        dqdata = inputdata[2]

    elif inputtype == 'MRSSKYFLAT':
        # input data is expected to be correctly formatted, with 
        # empty primary extension, and SCI, ERR, DQ, and FIELD_DEF in
        # in subsequent extensions.
        flattype = 'SKYFLAT'
        flatdata = hdulist[1].data
        errdata = hdulist[2].data
        dqdata = hdulist[3].data
        field_def = hdulist[4].data        

    else:
        print('Invalid input type, cannot continue.')
        sys.exit(1)

    # Convert the flatfield FITS into a MIRI CDP flatfield
    with MiriFlatfieldModel( init=inputfile, flattype=flattype ) as oldproduct:
        # Make a copy of the data product just created.
        flatproduct = oldproduct.copy()

        # Copy over flat data
        flatproduct.data = flatdata
        
        # Copy over err and/or dq array, if available
        if not errdata is None:
            flatproduct.err = errdata
        if not dqdata is None:
            flatproduct.dq = dqdata
            if field_def is None:
                print('Warning: DQ array included, but no field definitions provided!')
        
        # Copy over field definitions for DQ array, if available
        if not field_def is None:
            flatproduct.field_def = field_def
        
        # Apply modifications required to meet CDP specifications
        if 'READOUT' in header and not flatproduct.meta.exposure.readpatt:
            flatproduct.meta.exposure.readpatt = header['READOUT']
        flatproduct.meta.filename_original = os.path.basename(inputfile)
        flatproduct.meta.filename = os.path.basename(outputfile)
        flatproduct.meta.origin = 'MIRI European Consortium'
        
        if verb:
            print(flatproduct)
        if makeplot:
            flatproduct.plot()

        flatproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del flatproduct


if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a flat field FITS file into "
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
    
    mainx(inputfile,verb,makeplot,overwrite)
    
