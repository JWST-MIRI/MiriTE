#!/usr/bin/env python
#
# :History:
# 
# 11 Mar 2014: Created.
# 19 Mar 2014: Minor clean-up and moved from asciitable to astropy.io.ascii.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS)
#

"""

Script `convert_mrs_straylight` creates a FITS file in standard MIRI
format from an MRS Straylight pixel mask. The script relies on the
MIRI MRS Straylight data product, which is based on the STScI data
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
        Plot the straylight pixel mask.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""

import optparse
import os, sys, time

import astropy.io.ascii as ascii
import astropy.io.fits as pyfits

from miri.datamodels.dqflags import reserved_flags
from miri.datamodels.miri_straylight_model import MiriMrsStraylightModel

def load_mrs_stray_mask( filename, metafile ):
    """
    
    Reads straylight mask from a FITS file and metadata from an ascii file.
    
    """
    # Read the primary data array from the FITS file.
    try:
        hdulist = pyfits.open(filename)
        maskdata = hdulist[0].data
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    # If available, load meta data from file
    if not metafile is None:
        metatable = ascii.read(metafile, data_start=1, delimiter='|')
        metadata = dict(list(zip(metatable['col1'].tolist(),metatable['col2'].tolist())))
    else:
        metadata = None
    #stop()
    return metadata, maskdata

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts an MRS Straylight pixel mask into "
    usage += "standard MIRI CDP format."
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot straylight pixel mask"
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

    # Read the header and data from given file.
    if verb:
        print("Reading %s" % inputfile)
    metadata, maskdata = load_mrs_stray_mask( inputfile, metafile )
        
    # Convert the DHAS mask into a MIRI CDP mask
    with MiriMrsStraylightModel( init=inputfile ) as oldproduct:
        
        # Make a copy of the data product just created.
        dataproduct = oldproduct.copy()
        
        # set the DQ flags
        dataproduct.field_def = reserved_flags

        # copy over the mask data and metadata (if available)
        dataproduct.mask = maskdata
        if not metadata is None:
            dataproduct.set_housekeeping_metadata( metadata['ORIGIN'], author=metadata['AUTHOR'], version=metadata['VERSION'] )
            dataproduct.set_exposure_metadata(metadata['READPATT'], None, None)
            dataproduct.set_instrument_metadata( metadata['DETECTOR'], band=metadata['BAND'] )

        # set original and new filename
        dataproduct.meta.filename_original = os.path.basename(inputfile)
        dataproduct.meta.filename = os.path.basename(outputfile)
        
        if verb:
            print(dataproduct)
            print(dataproduct.flags_table)
        if makeplot:
            dataproduct.plot()

        dataproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del oldproduct, dataproduct
