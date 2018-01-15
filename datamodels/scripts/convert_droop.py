#!/usr/bin/env python
#
# Script 'convert_droop'
#
# :History:
# 
# 20 Feb 2013: Created
# 26 Feb 2013: Removed "inputtype" input parameter. Added SUBARRAY as a 
#              header keyword to copy over.
# 25 Jun 2013: Added astropy.io.ascii as an alternative to asciitable.
# 22 Aug 2013: columnnames renamed to fieldnames.
# 03 Sep 2013: DETECTOR field removed from table in anticipation of CDP-2
#              delivery. Ensure there is a DETECTOR identification in the
#              header.
# 11 Sep 2013: Modified to no longer set keywords BAND and CHANNEL (not 
#              relevant for Droop); now skipping first line of inputfile, 
#              presumed to hold column names.
# 24 Feb 2014: Instrument name (INSTRUME) changed from meta.instrument.type to
#              meta.instrument.name.
# 03 Jun 2014: Using astropy.io.ascii instead of asciitable.
# 15 Oct 2014: Minor update to setting the required keywords
#
# @author: Vincent Geers (DIAS), Steven Beard (UKATC)
#

"""

Script `convert_droop` creates a FITS file in standard MIRI format 
from an input ASCII table. The script relies on the  MiriDroopModel
data product, which is based on the STScI data model. 

Optionally, it can read in meta data from a separate ASCII table. Keywords, 
taken from a pre-determined list of required keywords, that are present in 
the metadata are inserted into the output product. 

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
    --meta or -m:
        The path+name of an ASCII file containing the meta data to 
        insert into the product.

"""

from __future__ import division, print_function

import optparse
import os, sys, time

import astropy.io.ascii as ascii
import numpy as np

from miri.datamodels.miri_droop_model import MiriDroopModel

def load_droop_file( filename, metafile ):
    """
    
    Reads flux conversion data from an ascii file.
    
    """
    # Read the flux conversion data from file
    droopdata = np.array(ascii.read(filename, data_start=1, delimiter='\s',\
        names=MiriDroopModel.fieldnames))
    
    # If available, load meta data from file
    if not metafile is None:
        metatable = ascii.read(metafile, data_start=1, delimiter='|')
        metadata = dict(list(zip(metatable['col1'].tolist(),metatable['col2'].tolist())))
    else:
        metadata = None

    return metadata, droopdata

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
    usage += "Converts a droop table file into "
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
    metafile = options.metafile

    # Read the given file.
    if verb:
        print("Reading %s" % inputfile)  
    metadata, inputdata = load_droop_file( inputfile, metafile )

    # The DETECTOR keyword is compulsory.
    detector = metadata['DETECTOR']
    with MiriDroopModel( droop_table=inputdata, detector=detector ) as droopproduct:
        # default modifications for CDP specifications
        droopproduct.meta.filename_original = os.path.basename(inputfile)
        droopproduct.meta.filename = os.path.basename(outputfile)

        # Set required keywords, if metadata was provided
        if not metadata is None:
            droopproduct.meta.instrument.model = set_meta_key( droopproduct.meta.instrument.model, metadata, 'MODELNAM' )
            droopproduct.meta.instrument.detector_settings = set_meta_key( droopproduct.meta.instrument.detector_settings, metadata, 'DETSETNG' )
            droopproduct.meta.exposure.readpatt = set_meta_key( droopproduct.meta.exposure.readpatt, metadata, 'READPATT' )
            droopproduct.meta.subarray.name = set_meta_key( droopproduct.meta.subarray.name, metadata, 'SUBARRAY' )
            if detector == 'MIRIMAGE':
                droopproduct.meta.instrument.filter = set_meta_key( droopproduct.meta.instrument.filter, metadata, 'FILTER' )
            droopproduct.meta.version = set_meta_key( droopproduct.meta.version, metadata, 'VERSION' )
            droopproduct.meta.author = set_meta_key( droopproduct.meta.author, metadata, 'AUTHOR' )
            droopproduct.meta.origin = set_meta_key( droopproduct.meta.origin, metadata, 'ORIGIN' )
            droopproduct.meta.description = set_meta_key( droopproduct.meta.description, metadata, 'DESCRIP' )

        if verb:
            print(droopproduct)
        if makeplot:
            droopproduct.plot()

        droopproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
    
        del droopproduct
