#!/usr/bin/env python
#
# :History:
# 
# 18 Feb 2013: Created (with memory issues)
# 27 Feb 2013: Updated to handle new DHAS format, where primary 
#              HDU is empty, and SCI, ERR, and FITERR are now located
#              in extensions 1, 2, and 4; a DQ array is now present 
#              in extension 3, and a binary table FIELD_DEF in extension
#              5. (Vincent Geers, DIAS).
# 04 Mar 2013: Check for INTNUM and DARKFIT keywords in the metadata.
# 03 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
#              - a solution to the "Types in column 0 do not match" problem
#              suggested by Michael Droettboom at STScI.
# 21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
# 27 May 2015: Replaced pyfits with astropy.io.fits
#
# @author: Steven Beard (UKATC), Vincent Geers (DIAS)
#

"""

Script `convert_dark_reference` creates a FITS file in standard MIRI
format from a dark reference file in DHAS format. The script relies on the
MIRI dark reference data product, which is based on the STScI data
model.

The following command arguments are defined by position:

    inputfile[0]
        The path+name of the file to be read. Compulsory.
    outputfile[1]
        The path+name of the file to be written.
        Optional. Defaults to the same name as inputfile with "_out" appended.

The command also takes the following options:

    --maxgroups
        Limit the number of groups read from the input file.
        ONLY SPECIFY THIS LIMIT TO RECOVER FROM A LACK OF MEMORY.
    --verbose or -v:
        Generate more output.
    --plot or -p:
        Plot the dark data.
    --overwrite or -o:
        Overwrite any existing FITS file.
    --ignore or -i
        Do not generate the data quality information.
        SET TRUE ONLY TO RECOVER FROM A LACK OF MEMORY.     

"""



import optparse
import os, sys, time

import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.miri_dark_reference_model import MiriDarkReferenceModel

def load_dhas_dark( filename, maxgroups=None, ignore_quality=False ):
    """
        
    Read FM dark correction data in the format written by DHAS and return
    copies of the primary FITS header, the SCI, ERR, DQ and FIT arrays, 
    and the FIELD_DEF table contained.
                
    :Parameters:
    
    filename: str
        The name of a FITS file in DHAS format from which to extract the
        data.
    maxgroups: int, optional
        This parameter can be used to minimise the memory usage by limiting
        the data size to a maximum number of groups. By default, there is
        no limit.
    ignore_quality: bool, optional
        Set to True to skip the generation of the data quality array,
        which might be needed for very large data sets where memory is
        a premium. The default is False.
        
    :Returns:
    
    datacollection: tuple of pyfits Header + 4 numpy ndarrays
        (header, scidata, errdata, dqdata, fitdata, dqtable)
    
    """
    # Read a FITS file in DHAS format. The file contains an empty
    # primary HDU, SCI data in extension 1, the standard deviation 
    # (ERR data) in extension 2, the data quality flags (DQ) in
    # extension 3, the error on fit for a single plane in extension
    # 4, and a binary table (FIELD_DEF) explaining the DQ flag in 
    # extension 5.

    try:
        hdulist = pyfits.open(filename)
  
        # Read the primary header
        fitsheader = hdulist[0].header
        
        # Read the SCI data array
        if maxgroups is not None:
            # Only extract the first maxgroups groups from the file.
            scidata = hdulist[1].data[:maxgroups,:,:]
        else:
            # Extract all the data from the file
            scidata = hdulist[1].data

        # Read the ERR data array
        try:
            if len(hdulist) > 2:
                if maxgroups is not None:
                    errdata = hdulist[2].data[:maxgroups,:,:]
                else:
                    errdata = hdulist[2].data                
            else:
                errdata = None
        except MemoryError as e:
            # Trap a memory error and continue with the SCI array only.
            # TODO: use Python warnings module
            print("+++ERR array ignored due to MemoryError!")
            del errdata
            errdata = None

        # Read the DQ array
        if len(hdulist) > 3:
            dqdata = hdulist[3].data
        elif not ignore_quality:
            # Generate a data quality array from the information embedded
            # in the DHAS dark data. Zones where there is no dark data are
            # indicated by -1. Reference pixels are indicated by 0.
            firstplane = scidata[0,:,:]
            dqdata = np.zeros(firstplane.shape, dtype=np.uint8)
            nodark = np.where(firstplane == -1)
            dqdata[nodark] = 2
            refpix = np.where(firstplane == 0)
            dqdata[refpix] = 1
        else:
            dqdata = None
            
        # Read the FITERR array
        if len(hdulist) > 4:
            fitdata = hdulist[4].data
        else:
            fitdata = None
        
        # Read the FIELD_DEF table
        if len(hdulist) > 5:
            dqtable = hdulist[5].data
        else:
            dqtable = None
            
    except Exception as e:
        # If the file could not be read re-raise the exception
        # with a more meaningful error message.
        strg = \
            "%s: Could not read DHAS format dark FITS file.\n   %s" % \
            (e.__class__.__name__, e)
        raise IOError(strg)
        
    finally:
        try:
            hdulist.close()
            del hdulist
        except:
            pass
        
    # The main data array must be 4-D, not 3-D.
    if scidata.ndim < 4:
        newshape = list(scidata.shape)
        newshape = [1] + newshape
        scidata.shape = newshape

    return (fitsheader, scidata, errdata, dqdata, fitdata, dqtable)

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile\n"
    usage += "Converts a DHAS format dark reference file into "
    usage += "standard MIRI CDP format."
    parser = optparse.OptionParser(usage)
    parser.add_option("", "--maxgroups", dest="maxgroups", type="int",
                     default=None, help="Maximum number of groups to read"
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
                      help="If data quality array is missing, do not create one."
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

    maxgroups = options.maxgroups
    verb = options.verb
    makeplot = options.makeplot
    overwrite = options.overwrite
    ignore  = options.ignore
    if ignore is None:
        ignore = False

    # Read the header and data from given file.
    if verb:
        strg = "Reading %s" % inputfile
        if maxgroups is not None:
            strg += " (limited to the first %d groups)" % maxgroups
        if ignore:
            strg += " (ignoring the quality data)"
        print(strg)
    (header, scidata, errdata, dqdata, fitdata, dqtable) = \
        load_dhas_dark( inputfile, maxgroups=maxgroups, ignore_quality=ignore )

    # TODO: Delivered CDPs already contain a FIELD_DEF table, but with only 2 
    # columns instead of the expected 3; so redefining here:
    dqdef = [(0,  'bad',          'Bad data'),
             (1,  'highfiterr',    'High fit error'),
             (2,  'largeslope',    'Large slope'),
             (3,  'largenegslope', 'Large negative slope')
             ]

    # Convert the DHAS mask into a MIRI CDP mask
    # FIXME: Creating a MiriDarkReferenceModel with the existing dark data
    # file as a template creates a massive data structure which causes memory
    # problems.
    #with MiriDarkReferenceModel( init=inputfile ) as oldproduct:
    # Instead we create a blank product with the required size. It is
    # necessary to copy all the metadata over by hand.
    with MiriDarkReferenceModel( init=scidata.shape ) as darkproduct:
        # Make a copy of the data product just created.
        # FIXME: copy generates a memory error but not copying generates
        # a storage space error! Catch 22. Can't use input file as template.
        #darkproduct = oldproduct.copy()
        darkproduct.data = scidata
        darkproduct.err = errdata
        darkproduct.dq = dqdata
        if dqtable is None:
            darkproduct.field_def = dqdef
        else:
            darkproduct.field_def = dqtable
        darkproduct.fiterr = fitdata
        
        # Copy over the metadata
        #darkproduct.meta.date = header['DATE']
        darkproduct.meta.version = header['VERSION']
        darkproduct.meta.author = header['AUTHOR']
        darkproduct.meta.instrument.model = header['MODELNAM']
        darkproduct.meta.instrument.detector = header['DETECTOR']
        
        if 'READPATT' in header:
            darkproduct.meta.exposure.readpatt = header['READPATT']
        elif 'READOUT' in header and not darkproduct.meta.exposure.readpatt:
            darkproduct.meta.exposure.readpatt = header['READOUT']

        # Apply modifications required to meet CDP-3 specifications
        if darkproduct.meta.instrument.detector == 'IC':
            darkproduct.meta.instrument.detector = 'MIRIMAGE'
        elif darkproduct.meta.instrument.detector == 'IM':
            darkproduct.meta.instrument.detector = 'MIRIMAGE'
        elif darkproduct.meta.instrument.detector == 'SW':
            darkproduct.meta.instrument.detector = 'MIRIFUSHORT'
        elif darkproduct.meta.instrument.detector == 'LW':
            darkproduct.meta.instrument.detector = 'MIRIFUSHORT'

        if darkproduct.meta.instrument.detector == 'MIRIFUSHORT' or \
           darkproduct.meta.instrument.detector == 'MIRIFULONG':
            if 'BAND' in header:
                darkproduct.meta.instrument.band = header['BAND']
            else:
                darkproduct.meta.instrument.band = 'ANY'
            if 'CHANNEL' in header:
                darkproduct.meta.instrument.channel = header['CHANNEL']
            else:
                darkproduct.meta.instrument.channel = 'ANY'

        if 'SUBARRAY' in header:
            darkproduct.meta.subarray.name = header['SUBARRAY']
        else:
            darkproduct.meta.subarray.name = 'FULL'
 
        darkproduct.meta.filename_original = os.path.basename(inputfile)
        darkproduct.meta.filename = os.path.basename(outputfile)
        darkproduct.meta.origin = 'MIRI European Consortium'
        
        if 'INTNUM' in header:
            darkproduct.meta.integration_number = header['INTNUM']
        if 'DARKFIT' in header:
            darkproduct.meta.fitted_after_frame = header['DARKFIT']
        
        if verb:
            print(darkproduct)
        if makeplot:
            darkproduct.plot()

        darkproduct.save( outputfile, overwrite=overwrite)
        if verb:
            print("Data saved to %s\n" % outputfile)
            
        del darkproduct
