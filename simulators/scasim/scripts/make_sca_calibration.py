#!/usr/bin/env python
#
# Script 'make_sca_calibration' generates either detector illumination
# data or calibration data using one of a set of defined data patterns
# and writes a FITS file in the format readable by the MIRI SCA
# simulator.
#
# :History:
#
# 04 Aug 2010: Created
# 06 Aug 2010: Documentation formatting problems corrected.
# 03 Sep 2010: illumination_maps renamed to data_maps
# 07 Sep 2010: BADPIXEL option added.
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: Unwanted option strings removed.
# 11 Nov 2010: DARKMAP pattern added and clone image maps created.
#              nbad, min, max and badvalue parameters added.
# 24 Nov 2010: Added the ability to create test data in a particular
#              subarray mode.
# 10 Jan 2011: Added optional scaid parameter.
# 13 Nov 2012: Major restructuring of the package folder.
#              Import statements updated.
# 18 Jun 2013: Modified to match the new parameter names used
#              by data_maps.
# 27 May 2015: Replaced pyfits with astropy.io.fits.
#              Corrected detector names.
# 08 Sep 2015: Removed dependency on legacy data maps.
#              create_test_data function moved here. Made compatible
#              with Python 3.
# 18 Sep 2015: Corrected typo in import and in some global variables.
#              Added BOX test data pattern.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 
# @author: Steven Beard (UKATC)
#

"""

Script 'make_sca_calibration' generates detector illumination data from
a defined pattern and writes a FITS file in the format readable by the
MIRI SCA simulator. This utility is useful for creating test data for
the SCA simulator - for example the CONSTANT pattern can generate a
perfect dark frame or flat-field. The compulsory command arguments are:

    rows
        The number of rows in the data.
    columns
        The number of columns in the data.
    outputfile
        The path+name of the SCA format FITS file to be created.

The following optional parameters may be provided by keyword:

    --pattern
        The type of pattern to be used

        * 'CONSTANT' generates a flat image with a constant value.
        * 'SLOPE' generates planar sloping data.
        * 'BOX' generates an image with a central bright box.
        * 'TESTIMAGE' generates a grid of test images
        * 'DARKMAP' generates a dark current multiplier map
        * 'BADPIXEL' generates a random bad pixel mask.

        If not provided the default is 'CONSTANT'.
    --constant
        Gives a constant value.
        For 'CONSTANT' the value used to fill the data array.
        For 'SLOPE' the constant to be added to all values.
        If not given, the default is 1.0.
    --rowslope
        For 'SLOPE', the data array slope along the rows.
        If not given, the default is 0.0.
    --colslope
        For 'SLOPE', the data array slope along the columns.
        If not given, the default is 0.0.
    --rowstep
        For 'TESTIMAGE', the test image row spacing.
        If not given, the default is rows/4.
    --colstep
        For 'TESTIMAGE', the test image column spacing.
        If not given, the default is columns/4.
    --nbad
        For DARKMAP, the number of random hot pixels to generate.
        For 'BADPIXEL' the number of random bad pixels to generate.
        If not given, the default is 0.
    --min
        For DARKMAP, the minimum DARK multipler.
        For TESTIMAGE or BOX, the background level.
        If not given, the default is 0.5.
    --max
        For 'DARKMAP', the maximum DARK multipler.
        For TESTIMAGE or BOX, the maximum brightness.
        If not given, the default is 2.0.
    --badvalue
        For 'DARKMAP', the dark multiplier for hot pixels.
        If not given, the default is 10000.0
    --wavmin
        Optionally, the minimum wavelength in microns associated
        with the data. A wavelength extension be created if wavmin
        and wavmax are both specified. Not used for 'DARKMAP' or
        'BADPIXEL' data.
    --wavmax
        Optionally, the maximum wavelength in microns associated
        with the data. A wavelength extension be created if wavmin
        and wavmax are both specified. Not used for 'DARKMAP' or
        'BADPIXEL' data.
    --scaid
        The identify of the sensor chip assembly module with which
        the calibration data are associated (e.g. 'MIRIMAGE', 'MIRIFULONG' or
        'MIRIFUSHORT'). Any digits provided will be written to an integer
        SCA_ID keyword in the FITS header of the output file. If
        not provided, no keyword will be written. 
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
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time
import re
import numpy as np
import astropy.io.fits as pyfits

from miri.datamodels import MiriMeasuredModel
from miri.datamodels.cdp import MiriBadPixelMaskModel
from miri.datamodels.sim import MiriIlluminationModel

def create_test_data(rows, columns, pattern, values, metadata=None,
                     cloneimage=None, add_wavelength=True,
                     wavmin=1.0, wavmax=30.0, seedvalue=None):
    """
    
    Create an data object object containing a test pattern,
    which can be useful for generating artificial dark and flat-field data,
    for example.
                 
    :Parameters:
    
    rows: int
        Number of rows
    columns: int
        Number of columns
    pattern: string
        The type of test pattern.
        CONSTANT - A constant value. useful for flat-field tests.
        SLOPE - A flat surface of constant slope. Useful for test data.
        BOX - An image with a central bright box.
        TESTIMAGE - A regular grid of bright images on a constant background.
        DARKMAP - A dark current multiplier map containing random hot pixels.
        BADPIXEL - A bad pixel map containing random dead pixels.
        Other types TBD.
    values: float or tuple of floats
        Values to be used to create the requested type of data.
        For CONSTANT:
        
            * values contains the constant value
            
        For SLOPE
        
            * values[0] contains the constant
            * values[1] contains the row slope
            * values[2] contains the column slope
            
        For BOX
        
            * values[0] contains the constant background
            * values[1] contains the peak brightness
            
        For TESTIMAGE
        
            * values[0] contains the constant background
            * values[1] contains the peak brightness
            * values[2] contains the row spacing in pixels
            * values[3] contains the column spacing in pixels
            
        For DARKMAP
        
            * values[0] contains the minimum dark multipler (normally <=1.0)
            * values[1] contains the maximum dark multipler (normally >=1.0)
            * values[2] contains the number of random hot pixels
            * values[3] contains the hot pixel multipler (normally >1000)
            
        For BADPIXEL
        
            * values contains the number of random dead pixels
            
    metadata: Metadata object, optional
        The primary metadata describing the data. If None, the metadata
        is ignored.
    cloneimage: array_like, optional
        A image to be cloned within the TESTIMAGE, DARKMAP or
        BADPIXEL patterns. If specified, some of the abnormal pixels
        will be stamped with this image. If not specified only single
        pixels are modified.
    add_wavelength: boolean, optional, default=True
        Add a test wavelength map to the data. Not valid for DARKMAP
        or BADPIXEL.
    wavmin: float, optional, default=1.0
        Minimum wavelength associated with data in microns.
    wavmax: float, optional, default=30.0
        Maximum wavelength associated with data in microns.
    seedvalue: int, optional, default=None
        The seed to be sent to the np.random number generator before
        generating the test data.
        If not specified, a value of None will be sent, which
        randomises the seed.
        
    :Raises:
        
    TypeError
        Raised if any of the parameters are of the wrong type, size
        or shape.
        
    :Returns:
    
    For the DARKMAP pattern, a new MeasuredModel object.
    For the BADPIXEL pattern, a new MiriBadPixelModel object.
    For other patterns, a new MiriIlluminationModel object.
        
    """
    # Set the seed for the np.random function.
    np.random.seed(seedvalue)

    # Test the requested data size
    try:
        rows = int(rows)
        columns = int(columns)
    except (TypeError, ValueError):
        strg = "Row and column values must be integers."
        raise TypeError(strg)

    # Define an intensity array containing the defined pattern.
    if pattern == 'CONSTANT':
        # A constant level.
        try:
            datavalues = float(values) * np.ones([rows,columns],
                                                 dtype=np.float32)
        except (TypeError, ValueError) as e:
            strg = "CONSTANT pattern needs a single floating point value"
            strg += "\n %s" % e
            raise TypeError(strg)
        
    elif pattern == 'SLOPE':
        # A flat, sloping surface.
        try:
            datavalues = \
                np.fromfunction(
                            lambda i,j: values[0] + i*values[1] + j*values[2],
                            [rows,columns])
        except (TypeError, ValueError, IndexError) as e:
            strg = "SLOPE pattern needs a tuple of 3 floats"
            strg += "\n %s" % e
            raise TypeError(strg)
        
    elif pattern == "BOX":
        # A box image.
        try:
            dmin = float(values[0])
            dmax = float(values[1])
            djump = dmax - dmin
        except (TypeError, ValueError, IndexError) as e:
            strg = "BOX pattern needs a tuple of (float, float)"
            strg += "\n %s" % e
            raise TypeError(strg)
        
        # Initialise the array full of background values and set a
        # central box to the maximum value.
        datavalues = dmin * np.ones([rows,columns], dtype=np.float32)
        
        rmin = rows/3
        rmax = rmin * 2
        cmin = columns/3
        cmax = cmin * 2
        datavalues[rmin:rmax, cmin:cmax] += djump
        
    elif pattern == "TESTIMAGE":
        # A grid of test images.
        try:
            dmin = float(values[0])
            dmax = float(values[1])
            djump = dmax - dmin
            rowstep = int(values[2])
            colstep = int(values[3])
        except (TypeError, ValueError, IndexError) as e:
            strg = "TESTIMAGE pattern needs a tuple of (float, float, int, int)"
            strg += "\n %s" % e
            raise TypeError(strg)
        
        # Initialise the array full of background values.
        datavalues = dmin * np.ones([rows,columns], dtype=np.float32)

        if cloneimage is not None:
            cloneimage = np.asarray(cloneimage)
            cloneimage *= djump
            lastrow = rows - cloneimage.shape[0] + 1
            lastcol = columns - cloneimage.shape[1] + 1
        else:
            lastrow = rows
            lastcol = columns

        for row in range(0, lastrow, rowstep):
            for col in range(0, lastcol, colstep):
                # Stamp either the clone image or a single pixel.
                if (cloneimage is not None):
                
                    for xx in range(0, cloneimage.shape[0]):
                        rr = row + xx
                        for yy in range(0, cloneimage.shape[1]):
                            cc = col + yy
                            datavalues[rr,cc] += cloneimage[xx,yy]
                else:
                    datavalues[row,col] = djump
        
    elif pattern == 'DARKMAP':
        # A dark current multipler map. There is never a wavelength array
        add_wavelength = False
        # The data contains 1.0 unless otherwise specified.
        datavalues = np.ones([rows,columns], dtype=np.float32)
        try:
            dmin = float(values[0])
            dmax = float(values[1])
            nhot = int(values[2])
            dhot = float(values[3])
        except (TypeError, ValueError, IndexError) as e:
            strg = "DARKMAP pattern needs a tuple of (float, float, int, float)"
            strg += "\n %s" % e
            raise TypeError(strg)
        
        # Rather than set every individual pixel (which would be inefficient)
        # create a map containing zones of different dark current multiplers.
        DARK_ZONES = 4
        rstep = datavalues.shape[0] // DARK_ZONES
        cstep = datavalues.shape[1] // DARK_ZONES
        
        for row in range(0, rows-rstep, rstep):
            for column in range(0, columns-rstep, cstep):
                if dmin < dmax:
                    rvalue = float(np.random.uniform(dmin, dmax))
                else:
                    # No variation
                    rvalue = dmin
                datavalues[row:row+rstep, column:column+rstep] = rvalue

        if cloneimage is not None:
            cloneimage = np.asarray(cloneimage)
            cloneimage *= dhot
            lastrow = rows - cloneimage.shape[0] + 1
            lastcol = columns - cloneimage.shape[1] + 1
        else:
            lastrow = rows
            lastcol = columns

        for zap in range(0, nhot):
            hit_row = int(np.random.uniform(0, lastrow))
            hit_column = int(np.random.uniform(0, lastcol))
            # Stamp only a small percent of the hot pixels with the clone image.
            MAX_CLONED_IMAGE_PERCENT = 5.0
            if (cloneimage is not None) and \
               (np.random.uniform(0, 100) < MAX_CLONED_IMAGE_PERCENT):
                
                for xx in range(0, cloneimage.shape[0]):
                    rr = hit_row + xx
                    for yy in range(0, cloneimage.shape[1]):
                        cc = hit_column + yy
                        datavalues[rr,cc] += cloneimage[xx,yy]
            else:
                datavalues[hit_row,hit_column] = dhot

    elif pattern == 'BADPIXEL':
        # A bad pixel map. There is never a wavelength array.
        add_wavelength = False
        datavalues = np.zeros([rows,columns], dtype=np.uint16)
        try:
            ndead = int(values)
        except (TypeError, ValueError) as e:
            strg = "BADPIXEL pattern needs an integer bad pixel count"
            strg += "\n %s" % e
            raise TypeError(strg)
        
        if cloneimage is not None:
            cloneimage = np.asarray(cloneimage)
            cloneimage *= _BAD_PIXEL
            lastrow = rows - cloneimage.shape[0] + 1
            lastcol = columns - cloneimage.shape[1] + 1
        else:
            lastrow = rows
            lastcol = columns

        for zap in range(0,ndead):
            hit_row = int(np.random.uniform(0, lastrow))
            hit_column = int(np.random.uniform(0, lastcol))
            # Stamp only a small percent of the bad pixels with the clone image.
            if (cloneimage is not None) and \
               (np.random.uniform(0, 100) < _MAX_CLONED_IMAGE_PERCENT):
                # TODO: I can't get this slicing to work.
#                uptorow = hit_row + cloneimage.shape[0] - 1
#                uptocol = hit_column + cloneimage.shape[1] - 1
#                datavalues[hit_row:uptorow,hit_column:uptocol] = cloneimage
                for xx in range(0, cloneimage.shape[0]):
                    rr = hit_row + xx
                    for yy in range(0, cloneimage.shape[1]):
                        cc = hit_column + yy
                        datavalues[rr,cc] = cloneimage[xx,yy]
            else:
                datavalues[hit_row,hit_column] = _BAD_PIXEL
                
#         flagvalues=(_GOOD_PIXEL,_BAD_PIXEL)
#         flagnames=('GOOD','BAD')
        
    # TODO: Add FRINGE pattern?
    else:
        # Other data types are not supported
        strg = "Sorry, test data pattern %s is not supported." % pattern
        raise ValueError(strg)
    
    # Add some test wavelength data if required. The wavelength will
    # increase linearly from bottom to top over the range specified.
    if add_wavelength:
        wavelength = np.empty_like(datavalues)
        wav = np.linspace(wavmin, wavmax, wavelength.shape[1])
        wavelength[:,:] = wav
#        wavelength = np.transpose(wavelength)
    else:
        wavelength = None
    
    data_metadata = None
    wavelength_metadata = None

    # Finally, create and return an appropriate object
    if pattern == 'BADPIXEL':
        data_object = MiriBadPixelMaskModel(dq=datavalues)
#         data_object = BadPixelMap(datavalues, flagvalues=flagvalues,
#                                   flagnames=flagnames, metadata=metadata,
#                                   dq_metadata=data_metadata)
    elif pattern == 'DARKMAP':
        data_object = MiriMeasuredModel(data=datavalues)          
#     elif pattern == 'FRINGEMAP':
#         data_object = MiriMeasuredModel(data=datavalues)          
    else:
        data_object = MiriIlluminationModel(intensity=datavalues,
                                            wavelength=wavelength)  
    return data_object

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] rows columns outputfile"
    usage +="\n\t[--pattern] [--constant] [--rowslope] [--colslope]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--pattern", dest="pattern", type="string",
                     default='CONSTANT', help="Data pattern"
                     )
    parser.add_option("", "--constant", dest="constant", type="float",
                     default=1.0, help="Constant data value"
                     )
    parser.add_option("", "--rowslope", dest="rowslope", type="float",
                     default=0.0, help="Slope in row direction"
                     )
    parser.add_option("", "--colslope", dest="colslope", type="float",
                     default=0.0, help="Slope in column direction"
                     )
    parser.add_option("", "--rowstep", dest="rowstep", type="int",
                     default=None, help="Number of row steps"
                     )
    parser.add_option("", "--colstep", dest="colstep", type="int",
                     default=None, help="Number of column steps"
                     )
    parser.add_option("", "--nbad", dest="nbad", type="int",
                     default=0, help="Number of bad or hot pixels"
                     )
    parser.add_option("", "--min", dest="min", type="float",
                     default=0.5, help="Minimum normal DARK multipler"
                     )
    parser.add_option("", "--max", dest="max", type="float",
                     default=2.0, help="Maximum normal DARK multipler"
                     )
    parser.add_option("", "--badvalue", dest="badvalue", type="float",
                     default=10000.0, help="DARK multipler for hot pixels"
                     )
    parser.add_option("", "--wavmin", dest="wavmin", type="float",
                     default=None, help="Minimum wavelength (microns)"
                     )
    parser.add_option("", "--wavmax", dest="wavmax", type="float",
                     default=None, help="Maximum wavelength (microns)"
                     )
    parser.add_option("", "--scaid", dest="scaid", type="string",
                     default='', help="Sensor Chip Assembly ID"
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
    pattern = options.pattern
    constant = float(options.constant)
    rowslope = float(options.rowslope)
    colslope = float(options.colslope)
    if options.rowstep is not None:
        rowstep = int(options.rowstep)
    else:
        rowstep = rows/4
    if options.colstep is not None:
        colstep = float(options.colstep)
    else:
        colstep = columns/4
    nbad = int(options.nbad)
    dmin = float(options.min)
    dmax = float(options.max)
    badvalue = float(options.badvalue)
    wavmin = options.wavmin
    wavmax = options.wavmax
    scaid = options.scaid
    subarray_str = options.subarray
        
    if (pattern != 'BADPIXEL') and \
       (wavmin is not None) and (wavmax is not None):
        add_wavelength = True
        wavmin = float(wavmin)
        wavmax = float(wavmax)
    else:
        add_wavelength = False

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

    # If a subarray has been specified, override the rows and columns
    # provided on the command line.
    if subarray_str:
        from miri.simulators.scasim import detector_properties
        try:
            subarray = detector_properties.SUBARRAY[subarray_str]
            if subarray is not None:
                rows = subarray[2]
                columns = subarray[3]
                if verbose > 0:
                    print( "Setting columns=%d, rows=%d." % (rows, columns) )
        except (KeyError, IndexError):
            strg = "Unrecognised or badly defined subarray mode: %s" % \
                subarray_str
            raise AttributeError(strg)

    keyw01 = pyfits.Card('HISTORY',
                         'Written by the make_sca_calibration script.',
                         '')
    str02 = "Test data pattern %s" % pattern
    kws02 = "%.64s" % str02
    keyw02 = pyfits.Card('HISTORY', kws02, '')
    if pattern == 'CONSTANT':
        str03 = "Constant value applied %.2f" % constant
        kws03 = "%.64s" % str03
        values = constant
    elif pattern == 'SLOPE':
        str03 = "Const=%.2f Rowslope=%.3g Colslope=%.3g" % \
            (constant, rowslope, colslope)
        kws03 = "%.64s" % str03
        values = (constant, rowslope, colslope)
    elif pattern == 'BADPIXEL':
        str03 = "Number of random bad pixels: %d" % nbad
        kws03 = "%.64s" % str03
        values = nbad
    elif pattern == 'BOX':
        str03 = "Min=%.2f Max=%.2f" % \
            (dmin, dmax)
        kws03 = "%.64s" % str03
        values = (dmin, dmax)
    elif pattern == 'TESTIMAGE':
        str03 = "Min=%.2f Max=%.2f test images spaced by %d,%d" % \
            (dmin, dmax, rowstep, colstep)
        kws03 = "%.64s" % str03
        values = (dmin, dmax, rowstep, colstep)
    elif pattern == 'DARKMAP':
        str03 = "Min=%.2f Max=%.2f %d hot pixels with %.1f" % \
            (dmin, dmax, nbad, badvalue)
        kws03 = "%.64s" % str03
        values = (dmin, dmax, nbad, badvalue)
    else:
        strg = "Unknown pattern name: %s" % pattern
        raise AttributeError(strg)
        
    keyw03 = pyfits.Card('HISTORY', kws03, '')
    pheader = pyfits.Header(cards=[keyw01, keyw02, keyw03])
    
    if scaid:
        # Extract only digits from the scaid to make an integer.
        digits = re.compile(r'([0-9]+)')
        scaid_int = digits.findall(scaid)
        if len(scaid_int) > 0 and scaid_int[0]:
            pheader.update("SCA_ID", int(scaid_int[0]),
                           "Sensor Chip Assembly ID")
    if subarray_str:
        pheader.update("SUBMODE", subarray_str, "Detector subarray mode")

    if verbose > 0:
        print( "Generating a %s illumination pattern into FITS file %s." % \
            (pattern, outputfile) )
        print( str03 )
        if add_wavelength:
            print( "A wavelength extension will be added covering the range " \
                "%.2f to %.2f microns." % (wavmin,wavmax) )

    # These particular image patterns look like the words
    # "bad" and "hot" when stamped on the data.
    if pattern == 'BADPIXEL':
        cloneimage = [[1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],
                      [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                      [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1],
                      [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
                      [1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0]]
    elif pattern == 'DARKMAP':
        cloneimage = [[1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0],
                      [1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0],
                      [1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0],
                      [1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1]]
    elif pattern == 'TESTIMAGE':
        # This pattern makes a "*" shape.
        cloneimage = [[0,   0,   0,   1,   0,   0,   0],
                      [0,   0.1, 0,   1,   0,   0.1, 0],
                      [0,   0,   0.2, 1,   0.2, 0,   0],
                      [1,   1,   1,   1,   1,   1,   1],
                      [0,   0,   0.2, 1,   0.2, 0,   0],
                      [0,   0.1, 0,   1,   0,   0.1, 0],
                      [0,   0,   0,   1,   0,   0,   0]]
    else:
        cloneimage = None
    
    # Create a new data map object from the data,
    map_object = create_test_data(rows, columns, pattern, values,
                                  metadata=pheader, cloneimage=cloneimage,
                                  add_wavelength=add_wavelength,
                                  wavmin=wavmin, wavmax=wavmax)
    if verbose > 1:
        print( map_object )

    # Save the illumination map to a FITS file.
    map_object.save(outputfile, overwrite=overwrite)
    
    if verbose > 0:
        print( "New illumination map saved." )
