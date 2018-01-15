#!/usr/bin/env python
#
# Script 'make_fringe_map' generates either detector illumination
# data or calibration data using one of a set of defined data patterns
# and writes a FITS file in the format readable by the MIRI SCA
# simulator.
#
# The working version of this script will use information on the refractive
# index and thickness of the detector, together with information on the
# intensity, wavelength and direction data from an illumination map to
# generate a detector fringe map.
#
# :History:
#
# 16 Oct 2014: Original dummy version.
# 27 May 2015: Replaced pyfits with astropy.io.fits
# 14 Aug 2015: Added path_difference functions.
# 08 Sep 2015: Made compatible with Python 3.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 17 Mar 2017: Replaced obsolete FringeMap data model with MiriFlatfieldModel.
# 
# @author: Steven Beard (UKATC)
#

"""

Script 'make_fringe_map' generates a MIRI fringe map and writes a FITS
file in the format readable by the MIRI SCA simulator.
The compulsory command arguments are:

    inputfile
        The path+name of the input file containing detector illumination
        information.
    outputfile
        The path+name of the fringe map FITS file to be created.

The following optional parameters may be provided by keyword:

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

import math

pi2 = math.pi * 2.0

from miri.datamodels.cdp import MiriFlatfieldModel

def refraction_angle( theta1, index2, index1=1.0):
    """
    
    Determine the exit angle (theta2) of a beam arriving at angle
    theta1 at a boundary between two media of refractive indices
    index1 and index2.
    
    The angle of refraction is determined by Snell's law,
    
        sin theta1 / sin theta2 = n2 / n1
        
    where n is the refractive index (named "index" in this function).
    
    :Parameters:
    
    theta1: float
        Arrival angle of incidence (radians)
    index2: float
        Refractive index of second medium.
    index1: float (optional)
        Refractive index of first medium. Defaults to a vacuum (1.0)
        
    :Returns:
    
    theta2: float
        Exit angle of incidence (radians)
    
    """
    assert index2 > 0.0
    sin_theta2 = math.sin(theta1) * index1 / index2
    theta2 = math.asin(sin_theta2)
    return theta2
    
def path_difference( theta1, index2, thickness, index1=1.0 ):
    """
    
    Determine the optical path difference between a ray reflected at
    the front detector surface and a ray which passes through and is
    reflected by the back surface.
    
    The optical path difference is given by the Fabry Perot formula,
    
        p = 2 n d cos theta
        
    where n is the refractive index (named "index" in this function)
    and d is the thickness of the medium (named "thickness" in this
    function).
    
    :Parameters:
    
    theta1: float
        Arrival angle of incidence (radians)
    index2: float
        Refractive index of detector material (second medium).
    thickness: float
        Thickness of the detector (microns)
    index1: float (optional)
        Refractive index of first medium. Defaults to a vacuum (1.0)
        
    :Returns:
    
    path_diff: float
        The path difference (in wavelength units). Fringes are maximum
        when this is a whole number of wavelengths.
        
    """
    assert index1 > 0.0
    # First calculate the angle of refraction, theta2, inside the medium.
    theta2 = refraction_angle(theta1, index2, index1=index1)
    
    # Then calculate the Fabry Perot path difference
    path_diff = 2.0 * index2 * thickness * math.cos(theta2)
    return path_diff

def phase_difference( wavelength, path_diff):
    """
    
    Calculate the phase difference from the wavelength and the optical
    path difference, given by the formula,
    
        phase = (2 pi / wavelength) * path_difference
        
    There is constructive interference when the phase difference is positive
    and destructive interference when the phase difference is negative.

    :Parameters:
    
    wavelength: float
        Wavelength (microns)
    path_diff: float
        The path difference (microns).
        
    :Returns:
    
    phase_diff: float
        The phase difference.
    
    """
    assert wavelength > 0.0
    phase_diff = pi2 * path_diff / wavelength
    return phase_diff

def generate_fringe_map( rows, columns ):
    """
    
    Generate a test fringe map.
    
    """
    
    wavelength = 10.0 # microns
    thickness = 25.0 # microns
    index2 = 3.6
    
    # The intensity incident angle, and location of each source
    source_data = [[0.1,85.0,500,500]]
#                    [0.1,85.0,550,550]]
#     source_angle = 0.0
#     source_flux = 1.0
    
    fringe_data = np.ones([rows,columns])
    
    for row in range(0,rows):
        for column in range (0,columns):
            for source in source_data:
                source_amplitude = source[0]
                source_angle = math.radians(source[1])
                source_row = source[2]
                source_column = source[3]
#                 print( "Adding source", source )
                sourcediff = math.sqrt((row-source_row)**2 + (column-source_column)**2)
                path_diff = sourcediff - path_difference(source_angle, index2, thickness)
                phase_diff = phase_difference(wavelength, path_diff)
                amplitude = source_amplitude * phase_diff
            
                fringe_data[row][column] += amplitude
            
    return fringe_data

def generate_test_map( rows, columns ):
    """
    
    Generate a test map containing a made up set of fringes.
    The test map is a beat pattern between a set of sources.
    
    """
    
#     wavelength = 10.0 # microns
#     thickness = 25.0 # microns
    
    # The intensity and location of each source
    source_data = [[0.02,400,400],
                   [0.05,450,450],
                   [0.1,500,500],
                   [0.1,550,550],
                   [0.05,600,600],
                   [0.02,650,650]]
#     source_angle = 0.0
#     source_flux = 1.0
    
    fringe_data = np.ones([rows,columns])
    
    for row in range(0,rows):
        for column in range (0,columns):
            for source in source_data:
                source_amplitude = source[0]
                source_row = source[1]
                source_column = source[2]
#                 print( "Adding source", source )
                sourcediff = math.sqrt((row-source_row)**2 + (column-source_column)**2)
                amplitude = source_amplitude * math.sin(sourcediff / 10.0)
            
                fringe_data[row][column] += amplitude
            
    return fringe_data

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile"
    parser = optparse.OptionParser(usage)
    
#     # Optional arguments (long option strings only).
#     parser.add_option("", "--subarray", dest="subarray", type="string",
#                      default='', help="Detector subarray mode"
#                      )

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
        
#     # Optional arguments.
#     subarray_str = options.subarray
        
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

#     # If a subarray has been specified, override the rows and columns
#     # provided in the input file.
#     if subarray_str:
#         from miri.scasim import detector_properties
#         try:
#             subarray = detector_properties.SUBARRAY[subarray_str]
#             if subarray is not None:
#                 rows = subarray[2]
#                 columns = subarray[3]
#                 if verbose > 0:
#                     print( "Setting columns=%d, rows=%d." % (rows, columns) )
#         except KeyError or IndexError:
#             strg = "Unrecognised or badly defined subarray mode: %s" % \
#                 subarray_str
#             raise AttributeError(strg)

    keyw01 = pyfits.Card('HISTORY',
                         'Written by the make_fringe_map script.',
                         '')
    pheader = pyfits.Header(cards=[keyw01])
    
#     if subarray_str:
#         pheader.update("SUBMODE", subarray_str, "Detector subarray mode")

    if verbose > 0:
        print( "Using %s to generate a fringe pattern into FITS file %s." % \
            (inputfile, outputfile) )

    # Create a new data map object from the data,
    fringe_data = generate_fringe_map( 1024, 1024 )
    fringe_map = MiriFlatfieldModel( fringe_data )
    #fringe_map = MiriFlatfieldModel( fringe_data, metadata=pheader )
    if verbose > 1:
        print( fringe_map )

    # Save the fringe map to a FITS file.
    fringe_map.save(outputfile, overwrite=overwrite)
    
    print( "New fringe map saved to %s." % outputfile )
