#!/usr/bin/env python
#
# Script 'make_sca_calibration' generates either detector illumination
# data or calibration data using one of a set of defined data patterns
# and writes a FITS file in the format readable by the MIRI SCA
# simulator.
#
# :History:
#
# 09 Jun 2020: Created
# 
# @author: Steven Beard (UKATC)
#

"""

Script 'make_sca_subarray' generates detector illumination data
for a defined subarray. The data contains a cross shape to test
the proper alignment of the subarray onto detector coordinates.

    rows
        The number of rows in the data.
    columns
        The number of columns in the data.
    outputfile
        The path+name of the SCA format FITS file to be created.

The following optional parameters may be provided by keyword:

    --min
        The baackground level of the illumination data.
        If not given, the default is 0.1
    --max
        The forground level of the illumination data.
        If not given, the default is 1.0.
    --subarray
        If this option is specified, it overrides the given number
        of rows and columns and creates data the same size as the
        specified subarray mode ('FULL', 'MASK1550', 'MASK1140',
        'MASK1065', 'MASKLYOT', 'BRIGHTSKY',  'SUB256', 'SUB128',
        'SUB64', 'SLITNESSPRISM')
   
The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.
    --overwrite or -o:
        Overwrite any existing SCA format FITS file.

"""

import optparse
import os, sys, time
import re
import numpy as np
#import astropy.io.fits as pyfits

#from miri.datamodels import MiriMeasuredModel
#from miri.datamodels.cdp import MiriBadPixelMaskModel
from miri.datamodels.sim import MiriIlluminationModel

def create_cross_data( nrows, ncolumns, minvalue=0.1, maxvalue=1.0,
                       addfifth=False ):
   """

   Create a numpy array containing a cross shape to test
   coordinate system alignment.

   """
   # Start with an array full of minvalue (background).
   data = minvalue * np.ones( [nrows, ncolumns] )

   # Incribe the 4 edges with maxvalue (foreground)
   data[0,:] = maxvalue
   data[-1,:] = maxvalue
   for row in range(0, nrows):
      data[row,0] = maxvalue
      data[row,-1] = maxvalue

   # Inscribe a cross from corner to corner
   for row in range(0, nrows):
      col = int( row * float(ncolumns)/float(nrows) )
      data[ row,  col] = maxvalue
      data[ row, -col] = maxvalue
      data[-row,  col] = maxvalue
      data[-row, -col] = maxvalue

      # If requested, identify the 5th column from the
      # left edge with a dotted line every 8th row.
      if addfifth:
         if (ncolumns > 4) and ((row % 8) == 0):
            data[ row, 4 ] = maxvalue

   return data

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] rows columns outputfile"
    usage +="\n\t[--min] [--max] [--subarray]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--min", dest="min", type="float",
                     default=0.1, help="Minimum normal DARK multipler"
                     )
    parser.add_option("", "--max", dest="max", type="float",
                     default=1.0, help="Maximum normal DARK multipler"
                     )
    parser.add_option("", "--subarray", dest="subarray", type="string",
                     default='', help="Detector subarray mode"
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
        rows = int(args[0])
        columns = int(args[1])
        outputfile = args[2]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    # Optional arguments.
    minvalue = float(options.min)
    maxvalue = float(options.max)
    subarray_str = options.subarray
    overwrite = options.overwrite
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if options.verb: verbose = 2
    if options.debug: verbose = 4
    if options.silent: verbose = 0

    # If a subarray has been specified, override the rows and columns
    # provided on the command line.
    addfifth = False
    if subarray_str:
        if subarray_str != "BRIGHTSKY" and subarray_str != "SUB256":
            # All subarrays except BRIGHTSKY and SUB256 touch the left
            # edge. Identify the fifth column.
            addfifth = True
        from miri.simulators.scasim import detector_properties
        try:
            subarray = detector_properties.SUBARRAY[subarray_str]
            if subarray is not None:
                rows = subarray[2]
                columns = subarray[3]
                if verbose > 0:
                    print( "Defining test data for subarray %s." % subarray_str)
                    print( "Setting columns=%d, rows=%d." % (rows, columns) )
        except (KeyError, IndexError):
            strg = "Unrecognised or badly defined subarray mode: %s" % \
                subarray_str
            raise AttributeError(strg)

    # Create a new illumination data
    datavalues = create_cross_data( rows, columns, minvalue=minvalue, maxvalue=maxvalue,
                                    addfifth=addfifth )
    wdata = 5.0 * np.ones_like(datavalues)
    map_object = MiriIlluminationModel(intensity=datavalues, wavelength=wdata)
    map_object.set_instrument_metadata( "MIRIMAGE" )
    if subarray_str is not None:
        map_object.set_subarray_metadata( subarray_str )
 
    if verbose > 1:
        print( map_object )
        map_object.plot()

    # Save the illumination map to a FITS file.
    map_object.save(outputfile, overwrite=overwrite)
    
    if verbose > 0:
        print( "New illumination map saved to %s." % outputfile )
