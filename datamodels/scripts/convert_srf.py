#!/usr/bin/env python
#
# Script 'convert_rsrf'
#
# :History:
# 
# 15 Oct 2014: Created, based off convert_rsrf.py, intended for converting
#              MIRI LRS spectral response function CDP.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS)
#

"""

Script `convert_srf` creates a FITS file in standard MIRI format 
from input ASCII table. The script relies on the MiriLrsFluxconversionModel 
data product, which is based on the STScI data model. 

Optionally, it can read in meta data from a separate ASCII table. Keywords, 
taken from a pre-determined list of required keywords, that are present in 
the metadata are inserted into the output product. 

The following command arguments are defined by position:

    inputtype[0]
        The type of distortion file. Valid values are LRSSRF.
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
    --meta or -m:
        The path+name of an ASCII file containing the meta data to 
        insert into the product.

"""

import optparse
import os, sys, time

import astropy.io.ascii as ascii
import astropy.io.fits as pyfits
import numpy as np

from miri.datamodels.miri_fluxconversion_models import \
    MiriLrsFluxconversionModel

def load_lrs_srf_file( filename, metafile ):
    """
    
    Reads SRF parameters from an ascii file.
    
    """
    # Read the SRF data from file
    srfdata = np.array(ascii.read(filename, data_start=3, delimiter='\s', \
        names=MiriLrsFluxconversionModel.fieldnames))
    
    # If available, load meta data from file
    if not metafile is None:
        metatable = ascii.read(metafile, data_start=1, delimiter='|')
        metadata = dict(list(zip(metatable['col1'].tolist(),metatable['col2'].tolist())))
    else:
        metadata = None
    
    return metadata, srfdata

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputtype inputfile outputfile\n"
    usage += "Converts an SRF file into "
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
    parser.add_option("-m", "--metafile", dest="metafile", action="store",
                      help="Filename for meta data."
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
    metafile = options.metafile

    # Read the given file.
    if verb:
        print("Reading %s" % inputfile)  
    
    if inputtype == 'LRSSRF':
        metadata, inputdata = load_lrs_srf_file( inputfile, metafile )
        
        # Convert the RSRF data into a MIRI CDP product
        dataproduct = MiriLrsFluxconversionModel( flux_table=inputdata )
        
        # Copy over missing meta data
        if not metadata is None:
            dataproduct.set_fits_keyword('MODELNAM', metadata['MODELNAM'])
            dataproduct.set_fits_keyword('DETECTOR', metadata['DETECTOR'])
            dataproduct.set_fits_keyword('DETSETNG', metadata['DETSETNG'])
            dataproduct.set_fits_keyword('READPATT', metadata['READPATT'])
            dataproduct.set_fits_keyword('SUBARRAY', metadata['SUBARRAY'])
            dataproduct.set_fits_keyword('FILTER', metadata['FILTER'])
            dataproduct.set_fits_keyword('VERSION', metadata['VERSION'])
            dataproduct.set_fits_keyword('AUTHOR', metadata['AUTHOR'])
            dataproduct.set_fits_keyword('ORIGIN', metadata['ORIGIN'])
            dataproduct.set_fits_keyword('DESCRIP', metadata['DESCRIP'])
    else:
        raise ValueError("Unrecognised input type %s" % inputtype)
         
    # default modifications for CDP-3 specifications
    dataproduct.meta.filename_original = os.path.basename(inputfile)
    dataproduct.meta.filename = os.path.basename(outputfile)

    if verb:
        print(dataproduct)
    if makeplot:
        dataproduct.plot()

    dataproduct.save( outputfile, overwrite=overwrite)
    if verb:
        print("Data saved to %s\n" % outputfile)
    
    del dataproduct
