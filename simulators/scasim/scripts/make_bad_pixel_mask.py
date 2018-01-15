#!/usr/bin/env python
#
# Script `make_bad_pixel_mask` reads bad pixel mask data either from
# an ASCII file or from a standard image format file (JPEG, TIFF, GIF, PNG,
# PPM, etc...) and writes a FITS file in the format readable by the MIRI SCA
# simulator.
#
# :History:
# 
# 07 Sep 2010: Created
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: Unwanted option strings removed.
# 10 Jan 2011: Added optional scaid parameter.
# 21 Sep 2011: Corrected some typos in comments.
# 13 Nov 2012: Major restructuring of the package folder.
#              Import statements updated.
# 18 Jun 2013: Modified to match the new parameter names used
#              by data_maps.
# 27 May 2015: Replaced pyfits with astropy.io.fits
# 08 Sep 2015: Removed dependency on legacy data maps. Made compatible
#              with Python 3.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
#
#@author: Steven Beard (UKATC)

"""

The make_bad_pixel_mask script reads bad pixel data from a file in the
specified format (IMAGE or ASCII) and writes a FITS file in the format
readable by the MIRI SCA simulator. This utility is useful for creating
test data for the SCA simulator or for creating masks with specific
pixels blacked out. The compulsory command arguments are:

    inputfile
        The path+name of the file to be read.
    outputfile
        The path+name of the bad pixel FITS file to be created.

The following optional parameters may be provided by keyword:

    --filetype
        The type of input file provided.

        * If 'ASCII', the bad pixel map is read from an ASCII file.
          ASCII files describing a bad pixel map only need to define the
          non-zero pixels. The ASCII file should contain one row per
          pixel, and each row contains 3 columns:
          
          * The row number of the pixel (starting at 0).
          * The column number of the pixel (starting at 0).
          * The value to be stored in the pixel at (row, column).
          
          Pixels can be defined in any order. The size of the bad pixel
          map is inferred from the range of values contained in the
          first two columns. The bad pixel map size can be defined
          explicitly by setting the first line to "rows-1, ncolumns-1,
          0". For example, if a data set is 1024 rows and 1032 columns,
          the first line would be "1023, 1031, 0".  
        * If 'IMAGE', bad pixel data are read from a standard image
          format file (such as JPEG, GIF or PPM) whose contents are
          converted into a composite image and then thresholded into
          discrete levels to make the bad pixel mask. This is useful if
          the image contains a picture of a bad pixel mask cut from a
          graphic or created using an image editor.
        * A file type of 'FITS' is also accepted, but converting FITS
          format data into FITS format data is not useful.
          
        If a file type is not specified, the default is 'ASCII'.
    --levels
        For an 'IMAGE' file only, the number of thresholding levels
        used to create the mask. The default is 2, which will create a
        mask full of 0s and 1s..
    --scaid
        The identify of the sensor chip assembly module with which
        the bad pixel mask is associated (e.g. 'MIRIMAGE', 'MIRIFULONG'
        or 'MIRIFUSHORT'). The name will be written to the SCA_ID
        keyword in the FITS header of the output file. If not provided,
        no keyword will be written. 

The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.
    --overwrite or -o:
        Overwrite any existing bad pixel FITS file.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time
import re
import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels.cdp import MiriBadPixelMaskModel

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile"
    usage += "\n\t[--filetype] [--levels]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--filetype", dest="filetype", type="string",
                     default='ASCII', help="File type"
                     )
    parser.add_option("", "--levels", dest="levels", type="int",
                     default=2, help="Number of threshold levels"
                     )
    parser.add_option("", "--scaid", dest="scaid", type="string",
                     default='', help="Sensor Chip Assembly ID"
                     )
     
    # Boolean options (short and long option strings).
    parser.add_option("-d", "--debug", dest="debug", action="store_true",
                      help="Debugging mode"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-s", "--silent", dest="silent", action="store_true",
                      help="Silent mode"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite existing file"
                     )

    (options, args) = parser.parse_args()
    
    # Compulsory arguments
    try:
        inputfile = args[0]
        outputfile = args[1]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    filetype = options.filetype
    if filetype == 'ASCII':
        levels = int(options.levels)
        flagvals = list(range(0,levels+1))
    else:
        flagvals = None
    scaid = options.scaid

    # Boolean flags.
    overwrite = options.overwrite
    verb = options.verb
    debug = options.debug
    silent = options.silent
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb: verbose = 2
    if debug: verbose = 4
    if silent: verbose = 0

    # Note: This header info is only added when data are read from an
    # ASCII or JPEG file. If read from a FITS file, the header contained
    # in the FITS file is used.
    keyw01 = pyfits.Card('HISTORY',
                         'Written by the make_bad_pixel_mask script.',
                         '')
    str02 = "Test data %s file %s" % (filetype, inputfile)
    kws02 = "%.64s" % str02
    keyw02 = pyfits.Card('HISTORY', kws02, '')
    pheader = pyfits.Header(cards=[keyw01, keyw02])
    
    if scaid:
        # Extract only digits from the scaid to make an integer.
        digits = re.compile(r'([0-9]+)')
        scaid_int = digits.findall(scaid)
        if len(scaid_int) > 0 and scaid_int[0]:
            pheader.update("SCA_ID", int(scaid_int[0]),
                           "Sensor Chip Assembly ID")

    if verbose > 0:
        print("Converting the %s file %s into FITS file %s." % \
            (filetype, inputfile, outputfile))
        if filetype == 'ASCII':
            print("Data will be thresholded into %d levels." % levels)

    # Create a new bad pixel map object from the data
    mask_object = MiriBadPixelMaskModel(inputfile) 
    if verbose > 1:
        print(mask_object)

    # Save the bad pixel map to a FITS file.
    mask_object.save(outputfile, overwrite=overwrite)
    
    if verbose > 0:
        print("New bad pixel mask map saved.")
