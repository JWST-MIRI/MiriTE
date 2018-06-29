#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module amplifier_properties - Defines the properties of the MIRI SCA
detector readout amplifiers.

The properties are stored in a list of Python dictionaries.

NOTE: This legacy module is not currently used by SCASim. It may be used
in the future to simulate a detector where the differences between the
zones read out by different amplifiers is significant. MIRI CDPs don't
distinguish between amplifiers and are applied to the whole detector.

:History:

15 Jul 2010: Created
20 Jul 2010: Bias and reference levels changed to integer.
             Gain polynomial changed to more sensible coefficients.
             Maximum DN value added to simulate amplifier saturation.
27 Jul 2010: Amplifier area added for cosmic ray target calculation.
16 Aug 2010: There can now be several different sets of amplifiers,
             each associated with a different focal plane module.
             The amplifier list associated with each FPM can now be
             looked up in a dictionary.
24 Aug 2010: Read noise measurements now read from FITS files.
30 Aug 2010: Locate the configuration files using an absolute file path
             so that scasim may be run from any directory.
             Normal and reference amplifiers are now counted per
             detector.
02 Sep 2010: Cosmic ray parameters updated.
27 Sep 2010: Python environment for windows verified.
30 Sep 2010: Added GAINTYPE and reduced gain by a factor of 2.7 (following
             comments by Scott Friedman). What are the real gain values?
07 Jan 2011: Amplifier properties classified and looked up by SCA_ID
             rather than by FPM_ID.
04 Feb 2011: Polynomial gain coefficients reversed to make them compatible
             with the numpy.poly1d function.
03 Mar 2011: Linear gain examples added but commented out.
24 Oct 2011: Modified to use new "filesearching" module, which should
             make it easier for users to substitute their own configuration
             files.
25 Oct 2013: Amplifier properties classified and looked up by DETECTOR
             rather than by SCA_ID.
24 Feb 2014: Use find_simulator_file function to find auxiliary data files
             rather than searching PYTHONPATH (which takes a long time when
             the software is installed using MIRICLE/ureka).
01 Jul 2014: Corrected the amplifier gain and made it much more linear.
09 Mar 2015: Note that this file simulates a variation in amplifier gain,
             but the CDP-3 GAIN reference file at the moment contains the
             same value (5.5 = 1/1.1818...)
11 Mar 2015: Local detector naming convention updated from 493,494,495 to IM,LW,SW.
08 Sep 2015: Made compatible with Python 3

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


from miri.simulators.find_simulator_file import find_simulator_file

# 
# MIRI contains three focal plane modules, each of which contains a
# sensor chip array. Each MIRI SCA uses 5 amplifiers. 4 of the amplifiers
# are assigned to columns on the detector (including the reference
# columns at the left and right edges) while a 5th amplifier is assigned
# to a bank of reference pixels.
#
# The following kinds of gain function (GAINTYPE) are supported. The
# meaning of the GAIN coefficients varies with GAINTYPE:
#
# LINEAR - A single gain coefficient containing the linear gain.
#
# POLYNOMIAL - A set of polynomial coefficients in the order
# accepted by the numpy.poly1d function, e.g. (2nd order, 1st order,
# zeroth order).
#
# POLYCUBE - A set of polynomial coefficients which vary from pixel
# to pixel. NOT FULLY IMPLEMENTED YET.
#
# The read noise measurement for each amplifier is obtained from a file.
# At the moment all the amplifiers belonging to the same FPM share the
# same file and have the same read noise, but it is possible to specify
# a different read noise for every amplifier (if needed).
#
#
# The find_simulator_file function searches for a named file within a
# search path of simulator data files(starting with the current working
# directory) and returns the absolute path.
#
# TODO: Replace with the measured FM/JPL gain and non-linearity settings.
_amp493_1 = {}
_amp493_1['ID'] = 'SCA493_1'
_amp493_1['NAME'] = "Amplifier 1"
_amp493_1['COMMENTS'] = "Reads every 4th detector column starting at column 0"
_amp493_1['TYPE'] = "illuminated"
_amp493_1['REGION'] = "0"   # Starting column
_amp493_1['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp493_1['BIAS'] = 2       # Bias in electrons
_amp493_1['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp493_1['GAIN'] = (-1.6e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp493_1['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp493_1['GAIN'] = 0.1818               # Coefficients for gain function
_amp493_1['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp493_1['READ_NOISE_FILE'] = find_simulator_file('read_noiseIM.fits')
_amp493_1['READ_NOISE_COL'] = 1  # Column number where to find read noise

_amp493_2 = {}
_amp493_2['ID'] = 'SCA493_2'
_amp493_2['NAME'] = "Amplifier 2"
_amp493_2['COMMENTS'] = "Reads every 4th detector column starting at column 1"
_amp493_2['TYPE'] = "illuminated"
_amp493_2['REGION'] = "1"   # Starting column
_amp493_2['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp493_2['BIAS'] = 2       # Bias in electrons
_amp493_2['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp493_2['GAIN'] = (-1.6e-8, 0.1820, 0.0)   # Coefficients for gain function
#_amp493_2['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp493_2['GAIN'] = 0.1818               # Coefficients for gain function
_amp493_2['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp493_2['READ_NOISE_FILE'] = find_simulator_file('read_noiseIM.fits')
_amp493_2['READ_NOISE_COL'] = 2  # Column number where to find read noise

_amp493_3 = {}
_amp493_3['ID'] = 'SCA493_3'
_amp493_3['NAME'] = "Amplifier 3"
_amp493_3['COMMENTS'] = "Reads every 4th detector column starting at column 2"
_amp493_3['TYPE'] = "illuminated"
_amp493_3['REGION'] = "2"   # Starting column
_amp493_3['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp493_3['BIAS'] = 2       # Bias in electrons
_amp493_3['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp493_3['GAIN'] = (-1.0e-8, 0.1746, 0.0)   # Coefficients for gain function
#_amp493_3['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp493_3['GAIN'] = 0.1818               # Coefficients for gain function
_amp493_3['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp493_3['READ_NOISE_FILE'] = find_simulator_file('read_noiseIM.fits')
_amp493_3['READ_NOISE_COL'] = 3  # Column number where to find read noise

_amp493_4 = {}
_amp493_4['ID'] = 'SCA493_4'
_amp493_4['NAME'] = "Amplifier 4"
_amp493_4['COMMENTS'] = "Reads every 4th detector column starting at column 3"
_amp493_4['TYPE'] = "illuminated"
_amp493_4['REGION'] = "3"   # Starting column
_amp493_4['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp493_4['BIAS'] = 2       # Bias in electrons
_amp493_4['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp493_4['GAIN'] = (-1.6e-8, 0.1809, 0.0)   # Coefficients for gain function
#_amp493_4['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp493_4['GAIN'] = 0.1818               # Coefficients for gain function
_amp493_4['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp493_4['READ_NOISE_FILE'] = find_simulator_file('read_noiseIM.fits')
_amp493_4['READ_NOISE_COL'] = 4  # Column number where to find read noise

_amp493_5 = {}
_amp493_5['ID'] = 'SCA493_5'
_amp493_5['NAME'] = "Amplifier 5"
_amp493_5['COMMENTS'] = "Reads from a separate bank of reference pixels"
_amp493_5['TYPE'] = "reference"
_amp493_5['REGION'] = "top" # Arranged at the top of the data (level 1 FITS)
_amp493_5['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp493_5['BIAS'] = 2       # Bias in electrons
_amp493_5['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp493_5['GAIN'] = (-1.0e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp493_5['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp493_5['GAIN'] = 0.1818               # Coefficients for gain function
_amp493_5['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp493_5['READ_NOISE_FILE'] = find_simulator_file('read_noiseIM.fits')
_amp493_5['READ_NOISE_COL'] = 5  # Column number where to find read noise

#----------------------------------------------------------------------------
_amp494_1 = {}
_amp494_1['ID'] = 'SCA494_1'
_amp494_1['NAME'] = "Amplifier 1"
_amp494_1['COMMENTS'] = "Reads every 4th detector column starting at column 0"
_amp494_1['TYPE'] = "illuminated"
_amp494_1['REGION'] = "0"   # Starting column
_amp494_1['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp494_1['BIAS'] = 2       # Bias in electrons
_amp494_1['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp494_1['GAIN'] = (-1.6e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp494_1['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp494_1['GAIN'] = 0.1818               # Coefficients for gain function
_amp494_1['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp494_1['READ_NOISE_FILE'] = find_simulator_file('read_noiseLW.fits')
_amp494_1['READ_NOISE_COL'] = 1  # Column number where to find read noise

_amp494_2 = {}
_amp494_2['ID'] = 'SCA494_2'
_amp494_2['NAME'] = "Amplifier 2"
_amp494_2['COMMENTS'] = "Reads every 4th detector column starting at column 1"
_amp494_2['TYPE'] = "illuminated"
_amp494_2['REGION'] = "1"   # Starting column
_amp494_2['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp494_2['BIAS'] = 2       # Bias in electrons
_amp494_2['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp494_2['GAIN'] = (-1.6e-8, 0.1820, 0.0)   # Coefficients for gain function
#_amp494_2['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp494_2['GAIN'] = 0.1818               # Coefficients for gain function
_amp494_2['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp494_2['READ_NOISE_FILE'] = find_simulator_file('read_noiseLW.fits')
_amp494_2['READ_NOISE_COL'] = 2  # Column number where to find read noise

_amp494_3 = {}
_amp494_3['ID'] = 'SCA494_3'
_amp494_3['NAME'] = "Amplifier 3"
_amp494_3['COMMENTS'] = "Reads every 4th detector column starting at column 2"
_amp494_3['TYPE'] = "illuminated"
_amp494_3['REGION'] = "2"   # Starting column
_amp494_3['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp494_3['BIAS'] = 2       # Bias in electrons
_amp494_3['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp494_3['GAIN'] = (-1.0e-8, 0.1746, 0.0)   # Coefficients for gain function
#_amp494_3['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp494_3['GAIN'] = 0.1818               # Coefficients for gain function
_amp494_3['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp494_3['READ_NOISE_FILE'] = find_simulator_file('read_noiseLW.fits')
_amp494_3['READ_NOISE_COL'] = 3  # Column number where to find read noise

_amp494_4 = {}
_amp494_4['ID'] = 'SCA494_4'
_amp494_4['NAME'] = "Amplifier 4"
_amp494_4['COMMENTS'] = "Reads every 4th detector column starting at column 3"
_amp494_4['TYPE'] = "illuminated"
_amp494_4['REGION'] = "3"   # Starting column
_amp494_4['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp494_4['BIAS'] = 2       # Bias in electrons
_amp494_4['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp494_4['GAIN'] = (-1.6e-8, 0.1809, 0.0)   # Coefficients for gain function
#_amp494_4['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp494_4['GAIN'] = 0.1818               # Coefficients for gain function
_amp494_4['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp494_4['READ_NOISE_FILE'] = find_simulator_file('read_noiseLW.fits')
_amp494_4['READ_NOISE_COL'] = 4  # Column number where to find read noise

_amp494_5 = {}
_amp494_5['ID'] = 'SCA494_5'
_amp494_5['NAME'] = "Amplifier 5"
_amp494_5['COMMENTS'] = "Reads from a separate bank of reference pixels"
_amp494_5['TYPE'] = "reference"
_amp494_5['REGION'] = "top" # Arranged at the top of the data (level 1 FITS)
_amp494_5['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp494_5['BIAS'] = 2       # Bias in electrons
_amp494_5['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp494_5['GAIN'] = (-1.0e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp494_5['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp494_5['GAIN'] = 0.1818               # Coefficients for gain function
_amp494_5['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp494_5['READ_NOISE_FILE'] = find_simulator_file('read_noiseLW.fits')
_amp494_5['READ_NOISE_COL'] = 5  # Column number where to find read noise

#----------------------------------------------------------------------------
_amp495_1 = {}
_amp495_1['ID'] = 'SCA495_1'
_amp495_1['NAME'] = "Amplifier 1"
_amp495_1['COMMENTS'] = "Reads every 4th detector column starting at column 0"
_amp495_1['TYPE'] = "illuminated"
_amp495_1['REGION'] = "0"   # Starting column
_amp495_1['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp495_1['BIAS'] = 2       # Bias in electrons
_amp495_1['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp495_1['GAIN'] = (-1.6e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp495_1['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp495_1['GAIN'] = 0.1818               # Coefficients for gain function
_amp495_1['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp495_1['READ_NOISE_FILE'] = find_simulator_file('read_noiseSW.fits')
_amp495_1['READ_NOISE_COL'] = 1  # Column number where to find read noise

_amp495_2 = {}
_amp495_2['ID'] = 'SCA495_2'
_amp495_2['NAME'] = "Amplifier 2"
_amp495_2['COMMENTS'] = "Reads every 4th detector column starting at column 1"
_amp495_2['TYPE'] = "illuminated"
_amp495_2['REGION'] = "1"   # Starting column
_amp495_2['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp495_2['BIAS'] = 2       # Bias in electrons
_amp495_2['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp495_2['GAIN'] = (-1.6e-8, 0.1720, 0.0)   # Coefficients for gain function
#_amp495_2['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp495_2['GAIN'] = 0.1818               # Coefficients for gain function
_amp495_2['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp495_2['READ_NOISE_FILE'] = find_simulator_file('read_noiseSW.fits')
_amp495_2['READ_NOISE_COL'] = 2  # Column number where to find read noise

_amp495_3 = {}
_amp495_3['ID'] = 'SCA495_3'
_amp495_3['NAME'] = "Amplifier 3"
_amp495_3['COMMENTS'] = "Reads every 4th detector column starting at column 2"
_amp495_3['TYPE'] = "illuminated"
_amp495_3['REGION'] = "2"   # Starting column
_amp495_3['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp495_3['BIAS'] = 2       # Bias in electrons
_amp495_3['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp495_3['GAIN'] = (-1.0e-8, 0.1646, 0.0)   # Coefficients for gain function
#_amp495_3['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp495_3['GAIN'] = 0.1818               # Coefficients for gain function
_amp495_3['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp495_3['READ_NOISE_FILE'] = find_simulator_file('read_noiseSW.fits')
_amp495_3['READ_NOISE_COL'] = 3  # Column number where to find read noise

_amp495_4 = {}
_amp495_4['ID'] = 'SCA495_4'
_amp495_4['NAME'] = "Amplifier 4"
_amp495_4['COMMENTS'] = "Reads every 4th detector column starting at column 3"
_amp495_4['TYPE'] = "illuminated"
_amp495_4['REGION'] = "3"   # Starting column
_amp495_4['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp495_4['BIAS'] = 2       # Bias in electrons
_amp495_4['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp495_4['GAIN'] = (-1.6e-8, 0.1609, 0.0)   # Coefficients for gain function
#_amp495_4['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp495_4['GAIN'] = 0.1818               # Coefficients for gain function
_amp495_4['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp495_4['READ_NOISE_FILE'] = find_simulator_file('read_noiseSW.fits')
_amp495_4['READ_NOISE_COL'] = 4  # Column number where to find read noise

_amp495_5 = {}
_amp495_5['ID'] = 'SCA495_5'
_amp495_5['NAME'] = "Amplifier 5"
_amp495_5['COMMENTS'] = "Reads from a separate bank of reference pixels"
_amp495_5['TYPE'] = "reference"
_amp495_5['REGION'] = "top" # Arranged at the top of the data (level 1 FITS)
_amp495_5['AREA'] = 5.0e5   # Area of amplifier surface in square microns
_amp495_5['BIAS'] = 2       # Bias in electrons
_amp495_5['GAINTYPE'] = 'POLYNOMIAL'         # Type of gain function
_amp495_5['GAIN'] = (-1.0e-8, 0.1818, 0.0)   # Coefficients for gain function
#_amp495_5['GAINTYPE'] = 'LINEAR'         # Type of gain function
#_amp495_5['GAIN'] = 0.1818               # Coefficients for gain function
_amp495_5['MAX_DN'] = 65535.0 # Maximum DN output by this amplifier
_amp495_5['READ_NOISE_FILE'] = find_simulator_file('read_noiseSW.fits')
_amp495_5['READ_NOISE_COL'] = 5  # Column number where to find read noise


# Lists of defined amplifiers for each FPM.
_amplist494 = (_amp494_1, _amp494_2, _amp494_3, _amp494_4, _amp494_5)
_amplist495 = (_amp495_1, _amp495_2, _amp495_3, _amp495_4, _amp495_5)
_amplist493 = (_amp493_1, _amp493_2, _amp493_3, _amp493_4, _amp493_5)


# The following constants may be imported by SCA simulator software modules.
#
# Dictionary of known focal plane modules. The list of amplifiers for
# each detector can be obtained by looking up its unique detector name
# in this dictionary (i.e. the overall data structure is a dictionary
# of lists of dictionaries).
AMPLIFIERS_DICT = {'MIRIMAGE'    : _amplist493,
                   'MIRIFULONG'  : _amplist494,
                   'MIRIFUSHORT' : _amplist495}

# Record the number of amplifiers reading out the illuminated zone
# of the detector and the number of amplifiers reading out the
# reference pixels (as defined by the amplifier type), for each detector.
N_ILLUM_AMPLIFIERS = {}
N_REF_AMPLIFIERS = {}
for detid in AMPLIFIERS_DICT:
    detector = AMPLIFIERS_DICT[detid]
    _icount = 0
    _rcount = 0
    for amp in detector:
        if amp['TYPE'] == 'illuminated':
            _icount += 1
        else:
            _rcount += 1
    N_ILLUM_AMPLIFIERS[detid] = _icount
    N_REF_AMPLIFIERS[detid] = _rcount

# The initial reference level in electrons and the maximum random
# drift in the reference level expected after each readout.
# (These parameters control how short term drifts in the amplifier
# output are simulated.)
INITIAL_REF_LEVEL = 5
MAX_REF_DRIFT = 3

# The following parameters describe how the amplifiers respond when
# hit by a cosmic ray which would have released NE electrons in the
# detector. The following possible effects are considered:
#
# (1) If the cosmic ray strikes an amplifier during an A/D conversion
# it may cause a glitch in the conversion. A high energy cosmic ray
# might affect the next few conversions made by that A/D and generate
# a streak of bad readings. These A/D glitches affect a single reading
# only and the amplifier returns to normal on the next reading.
#
# The COSMIC_RAY_GLITCH_SENSITIVITY parameter defines the sensitivity
# of the A/D converters to a cosmic ray strike. A strike with an energy
# equivalent to NE electrons will cause a streak of glitches
# COSMIC_RAY_GLITCH_SENSITIVITY x NE pixels long. Set this parameter to
# 0.0 to turn off glitches altogether.
COSMIC_RAY_GLITCH_SENSITIVITY = 1.0e-4 # A glitch for every 10000e energy.
#
# The COSMIC_RAY_GLITCH_IPC parameter defines how much cosmic ray energy
# communicated from one A/D conversion to the next. Set this parameter to
# 0.0 to generate single pixel glitches only or to 1.0 to make streaks
# of glitches with no decay in their effect.
COSMIC_RAY_GLITCH_IPC = 0.1 # Each subsequent glitch has 10% of the effect.
#
# (2) A cosmic ray may heat the amplifier slightly or induce unwanted
# currents that cause a temporary increase in read noise. This increased
# noise might linger for the next few readings but will go away the next
# time the detector is reset.
#
# The COSMIC_RAY_NOISE_SENSITIVITY parameter defines how much the read
# noise is affected. A strike with an energy equivalent to NE electrons
# will increase  the read noise by COSMIC_RAY_NOISE_SENSITIVITY x NE
# electrons. Set this parameter to 0.0 to turn off the read noise effect.
COSMIC_RAY_NOISE_SENSITIVITY = 1.0e-5  # 1e of extra noise per 100000e energy.
#
# The COSMIC_RAY_NOISE_DECAY parameter defines by what factor the excess
# read noise decays with each new reading (assuming an exponential decay).
# Set this parameter to 0.0 to make a cosmic ray affect a single reading
# only, or set it to zero to make the excess noise stay the same until
# the next detector reset.
COSMIC_RAY_NOISE_DECAY = 0.5 # Excess noise reduced by x0.5 with every readout.

if __name__ == '__main__':
    print( "NOTE: The AmplifierProperties module is supposed to be " \
        "imported by another module, not run as a main program." )
    print( "Initial reference level = %d and max drift = %d" % \
        (INITIAL_REF_LEVEL, MAX_REF_DRIFT) )
    print( "The following amplifiers are defined:" )
    for detid in AMPLIFIERS_DICT:
        print( "          DETECTOR %s\n-----------------------------" % detid )
        detector = AMPLIFIERS_DICT[detid]
        print( "Detector %s has %d normal and %d reference amplifiers." % \
            (detid, N_ILLUM_AMPLIFIERS[detid],
             N_REF_AMPLIFIERS[detid]) )
        for amp in detector:
            for key in amp:
                print( "%24s = %s" % (key, amp[key]) )
            print( "- - - - - - - - - - - - - - - - - - - -" )
    print( "Finished." )
