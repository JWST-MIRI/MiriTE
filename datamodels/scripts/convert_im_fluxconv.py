#!/usr/bin/env python
#
# Script 'convert_im_fluxconv'
#
# :History:
# 
# 20 Feb 2013: Created
# 26 Feb 2013: Added SUBARRAY as a header keyword to copy over.
# 25 Jun 2013: Added astropy.io.ascii as an alternative to asciitable.
# 22 Aug 2013: columnnames renamed to fieldnames.
# 31 Oct 2013: Corrected import.
# 24 Feb 2014: Instrument name (INSTRUME) changed from meta.instrument.type to
#              meta.instrument.name.
# 03 Jun 2014: Using astropy.io.ascii instead of asciitable.
# 15 Oct 2014: Minor CDP-3 update to setting the required keywords.
#
# @author: Vincent Geers (DIAS), Steven Beard (UKATC)
#

"""

Script `convert_im_fluxconv` creates a FITS file in standard MIRI format 
from an input ASCII table. The script relies on the 
MiriImagingFluxconversionModel data product, which is based on the STScI 
data model. 

Optionally, it can read in meta data from a separate ASCII table. Keywords, 
taken from a pre-determined list of required keywords, that are present in 
the metadata are inserted into the output product. 

The following command arguments are defined by position:

    inputtype[0]
        The type of flux conversion file. Valid values are IMFLUXCAL.
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
import numpy as np

from miri.datamodels.miri_fluxconversion_models import \
    MiriImagingFluxconversionModel

def load_fluxconv_file( filename, metafile ):
    """
    
    Reads flux conversion data from an ascii file.
    
    """
    # Read the flux conversion data from file
    fluxconvdata = np.array(ascii.read(filename, data_start=0, delimiter='\s', \
        names=MiriImagingFluxconversionModel.fieldnames))
    
    # If available, load meta data from file
    if not metafile is None:
        metatable = ascii.read(metafile, data_start=1, delimiter='|')
        metadata = dict(list(zip(metatable['col1'].tolist(),metatable['col2'].tolist())))
    else:
        metadata = None
    
    return metadata, fluxconvdata

def set_meta_key( product, metadata, keyword ):
    try:
        product = metadata[keyword]
    except KeyError:
        pass
    return product

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputtype inputfile outputfile\n"
    usage += "Converts an Imager flux conversion file into "
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
    if inputtype == 'IMFLUXCAL':
        metadata, inputdata = load_fluxconv_file( inputfile, metafile )

    with MiriImagingFluxconversionModel( flux_table=inputdata ) as fluxconvproduct:
        # default modifications for CDP-1 specifications
        fluxconvproduct.meta.filename_original = os.path.basename(inputfile)
        fluxconvproduct.meta.filename = os.path.basename(outputfile)

        # Based on inputtype, set required keywords
        if inputtype == 'IMFLUXCAL' and not metadata is None:
            fluxconvproduct.meta.instrument.model = set_meta_key( fluxconvproduct.meta.instrument.model, metadata, 'MODELNAM' )
            fluxconvproduct.meta.instrument.detector = set_meta_key( fluxconvproduct.meta.instrument.detector, metadata, 'DETECTOR' )
            fluxconvproduct.meta.instrument.detector_settings = set_meta_key( fluxconvproduct.meta.instrument.detector_settings, metadata, 'DETSETNG' )
            fluxconvproduct.meta.exposure.readpatt = set_meta_key( fluxconvproduct.meta.exposure.readpatt, metadata, 'READPATT' )
            fluxconvproduct.meta.subarray.name = set_meta_key( fluxconvproduct.meta.subarray.name, metadata, 'SUBARRAY' )
            fluxconvproduct.meta.instrument.filter = set_meta_key( fluxconvproduct.meta.instrument.filter, metadata, 'FILTER' )
            fluxconvproduct.meta.version = set_meta_key( fluxconvproduct.meta.version, metadata, 'VERSION' )
            fluxconvproduct.meta.author = set_meta_key( fluxconvproduct.meta.author, metadata, 'AUTHOR' )
            fluxconvproduct.meta.origin = set_meta_key( fluxconvproduct.meta.origin, metadata, 'ORIGIN' )
            fluxconvproduct.meta.description = set_meta_key( fluxconvproduct.meta.description, metadata, 'DESCRIP' )

        if verb:
            print(fluxconvproduct)
        if makeplot:
            fluxconvproduct.plot()

        fluxconvproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
    
        del fluxconvproduct
