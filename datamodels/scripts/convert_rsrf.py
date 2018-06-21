#!/usr/bin/env python
#
# Script 'convert_rsrf'
#
# :History:
# 
# 19 Feb 2013: Created
# 20 Feb 2013: Fixed errors in metadata names. Metadata keywords are now
#              set with set_meta_key function, which tries to set 
#              metadata but will pass without error if the keyword was 
#              not found in the provided metadata file.
# 21 Feb 2013: Updates to correctly handle the binary FITS table input 
#              for MRS RSRF. For MRS RSRF files, make a copy of the data
#              product before modifying it.
# 26 Feb 2013: Added SUBARRAY as a header keyword to copy over.
# 04 Mar 2013: Check for REFRSRF, REFERR and REFWAVEL keywords in the metadata.
# 25 Jun 2013: Added astropy.io.ascii as an alternative to asciitable.
# 22 Aug 2013: columnnames renamed to fieldnames.
# 31 Oct 2013: Corrected import.
# 24 Feb 2014: Instrument name (INSTRUME) changed from meta.instrument.type to
#              meta.instrument.name.
# 04 Mar 2014: Updated LRS RSRF to LRS SRF, renamed rsrfproduct to dataproduct.
#              Setting metadata with new convenience function; removed 
#              set_meta_key( product, metadata, keyword ) method.
# 03 Jun 2014: Using astropy.io.ascii instead of asciitable.
# 21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Vincent Geers (DIAS), Steven Beard (UKATC)
#

"""

Script `convert_distortion` creates a FITS file in standard MIRI format 
from input ASCII table. The script relies on the MiriDistortionModel data 
product, which is based on the STScI data model. 

Optionally, it can read in meta data from a separate ASCII table. Keywords, 
taken from a pre-determined list of required keywords, that are present in 
the metadata are inserted into the output product. 

The following command arguments are defined by position:

    inputtype[0]
        The type of distortion file. Valid values are LRSRSRF, MRSRSRF.
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
    MiriLrsFluxconversionModel, MiriMrsFluxconversionModel

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

def load_mrs_rsrf_file( filename ):
    """
    
    Reads RSRF parameters from binary FITS table.
    
    """
    # Read the RSRF data from file
    try:
        hdulist = pyfits.open(filename)
        fits_hdr_prim = hdulist[0].header
        fits_hdr_ext = hdulist[1].header
        #fitsdata = hdulist[1].data[0]
        
        fitsdata = np.core.records.fromarrays([ 
                        hdulist[1].data[0][0],
                        hdulist[1].data[0][1],
                        hdulist[1].data[0][2]], 
                        names=MiriMrsFluxconversionModel.fieldnames)
    finally:
        try:
            hdulist.close()
            del hdulist
        except Exception:
            pass

    return fits_hdr_prim, fits_hdr_ext, fitsdata

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputtype inputfile outputfile\n"
    usage += "Converts an (R)SRF file into "
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
            dataproduct.set_housekeeping_metadata(metadata['ORIGIN'], author=metadata['AUTHOR'], version=metadata['VERSION'] )
            dataproduct.set_instrument_metadata('MIRIMAGE', filt=metadata['FILTER'] )
            dataproduct.set_exposure_metadata(metadata['READPATT'], None, None)
            dataproduct.set_subarray_metadata(metadata['SUBARRAY'])

    elif inputtype == 'MRSRSRF':
        header_prim, header_ext, inputdata  = load_mrs_rsrf_file( inputfile )

        # Convert the RSRF data into a MIRI CDP product
        oldproduct = MiriMrsFluxconversionModel( init=inputfile )

        # Copy over RSRF data
        dataproduct = oldproduct.copy()
        dataproduct.flux_table = inputdata
    else:
        raise ValueError("Unrecognised input type %s" % inputtype)
         
    # default modifications for CDP-1 specifications
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
    
    if inputtype == 'MRSRSRF':
        del oldproduct
