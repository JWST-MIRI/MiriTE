#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module detector_properties - Defines the properties of the MIRI SCA
detector.

NOTE: Other JWST detectors can be simulated as long as their
properties can be described by the parameters contained here.
The main differences are the material, pixel size and detector
thickness.

NOTE: These properties have been defined to match the data obtained
during FM testing of the detectors. For the best simulation, the
parameters should be kept up to date as detector knowledge improves.

Each set of properties is stored in a Python dictionary and the
properties belonging to a particular focal plane module can be
looked up using a dictionary of dictionaries.

Sources:

(1) JPL D-25632, MIRI Operational Concept Document, 4 March 2010.
(2) JPL D-46944, MIRI Flight Focal Plane Module End Item Data Package
    (FPM EIDP), 10 May 2010.
(3) UA_MIRI_006, MIRI Test Report, Description of Bad Pixel Mask,
    Version 4, 28 August 2011.
(4) The Mid-Infrared Instrument for the James Webb Space Telescope,
    VIII: The MIRI Focal Plane System; M. E. Ressler et al.,
    Publications of the Astronomical Society of Pacific, Volume 127,
    Issue 953, pp. 675 (2015)
(5) JPL MIRI DFM 478 04.02, MIRI FPS Exposure Time Calculations (SCE FPGA2),
    M. E. Ressler, October 2014.

:History:

15 Jul 2010: Created
22 Jul 2010: Pixel size and thickness added.
25 Jul 2010: Well depth increased to 100000 (ref 2).
27 Jul 2010: Subarray options added.
01 Aug 2010: Added CHIP_ID and DEFAULTs.
09 Aug 2010: Renamed. QE data file renamed to qe_measurement.txt
16 Aug 2010: Detector properties compiled into a list, since
             MIRI contains three focal plane modules which can behave
             slightly differently. There are now separate dark current
             and QE measurements for each FPM.
24 Aug 2010: QE and dark current measurements now read from FITS files.
             Detector thickness increased to 470 microns.
30 Aug 2010: Locate the configuration files using an absolute file path
             so that scasim may be run from any directory.
03 Sep 2010: Added COSMIC_RAY_LEAKAGE_FRACTION. According to Mike Ressler,
             a small fraction of cosmic ray events can cause charge
             leakage rather than a charge jump.
             Bad pixel maps added.
27 Sep 2010: Python environment for windows verified.
01 Oct 2010: Dropped groups added to readout modes.
12 Oct 2010: Number of frames per group added to readout modes (always 1
             for MIRI but will allow SCAsim to be extended for other JWST
             detectors).
18 Oct 2010: Group gap for GAP modes changed from 4 to 8 to increase the
             exposure time that can be fitted into a reasonable sized
             output file.
11 Nov 2010: DARK_MAP parameters added.
15 Nov 2010: Mistake in CHIP_ID corrected.
22 Nov 2010: Corrections to FPM and SCA IDs in FITS header.
24 Nov 2010: Added test Subarray modes.
15 Dec 2010: SCA_ID values changed from 104,105,106 to 493,494,495
             (as expected by miri_cube)
07 Jan 2011: ID, SCA_ID and DETECTOR values updated to reflect the
             values reported by Tim Grundy (RAL) on 20 Dec 2010.
             Detector properties classified and looked up by SCA_ID
             rather than by FPM_ID.
10 Mar 2011: Subarray modes MASK1065, MASK1550 and SLITLESSPRISM
             corrected to match definitions in OCD Revision C
             (4 March 2010). Subarray modes LRS3, AXIS64, AXIS128 and
             AXIS256 added for test purposes. Subarray modes TEST32
             and TEST64 removed.
25 Mar 2011: PERSISTENCE parameter added, which can be a linear factor
             or a set of polynomial coefficients.
05 Apr 2011: Pixel response function added (not used yet).
21 Sep 2011: The bad pixel masks and dark maps derived from FM tests (3)
             are now used by default.
05 Oct 2011: PERSISTENCE altered for more realistic simulation of
             persistence effects.
             NOTE: Should FRAME_RESETS be increased to simulate the
             post-FM testing adjustment of detector controller parameters?
             NOTE: Which is the correct detector naming convention:
             Tim Grundy's or  Jane Morrison's?
24 Oct 2011: Modified to use new "filesearching" module, which should
             make it easier for users to substitute their own configuration
             USE_FM_MEASUREMENTS flag removed (too complicated).
             files. Pixel response function removed (never used).
13 Jun 2013: Changed detector and subarray names to match new data model.
25 Oct 2013: Make the primary keyword for determining the detector ID
             'DETECTOR' rather than 'SCA_ID'. Detector properties are
             looked up by detector name rather than sca_id.
24 Feb 2014: Use find_simulator_file function to find auxiliary data files
             rather than searching PYTHONPATH (which takes a long time when
             the software is installed using MIRICLE/ureka).
07 May 2014: Added zeropoint drift, linearity and latency factors.
             Added charge trapping and slope drift factors.
08 May 2014: Charge trapping parameter added.
05 Jun 2014: Removed charge trapping parameters and added slow and fast
             latency parameters. Decay parameters now given as a timescale.
19 Jun 2014: Slow and fast zeropoint drift parameters included.
27 Feb 2015: Subarray parameters brought up to date.
09 Mar 2015: Switch over to CDP-3 bad pixel masks. Well depth adjusted to
             match CDP-3 pixel saturation results (approximately).
11 Mar 2015: Switch over to CDP-3 dark masks. Local detector naming
             convention updated from 493,494,495 to IM,LW,SW.
21 May 2015: Amplifier level removed and gain and read noise defined by
             calibration data products.
06 Aug 2015: Noise calibration factor added.
08 Sep 2015: Made compatible with Python 3.
02 Oct 2015: Removed readout modes which are no longer available.
08 Dec 2015: Correction to subarray parameters to match Mike Ressler's PASP
             paper.
20 Mar 2016: Pedestal values added, to simulate electronics bias and keep
             the final simulated DN values away from zero. Bad pixels can
             be given a different pedestal value.
23 Mar 2016: Obtain MASK, PIXELFLAT and GAIN calibration from CDP files
             found by searching the CDP repository rather than by
             looking up a reference defined here. (The DARK map needs to
             remain here because it doesn't directly match any CDP file.)
05 May 2016: Detector sensitivity coefficients updated.
14 Jul 2016: Pedestal values adjusted.
02 Nov 2016: Legacy PERSISTENCE coefficients zeroed. Legacy DARK and bad
             pixel masks also removed.
14 Feb 2017: Zeropoint, pedestal level and noise factor adjusted.
20 Jun 2017: Mean gain is (e/DN), not (DN/e).
26 Jun 2017: DARK_CURRENT parameter added, to define the expected level of
             dark current in e/s.
19 Jul 2017: Well depth adjusted to match expected saturation levels.
             Nominal dark current changed from 0.12 to 0.21 e/s.
13 Oct 2017: New frame time calculation from Mike Ressler.
             SLOW mode now uses 8 out of 9 samples. READOUT_MODE now defines
             samplesum, sampleskip and refpixsampleskip parameters separately.
01 Nov 2017: Note that the detector drift and latent coefficients are only
             valid for FAST mode.
13 Dec 2017: DARK_MAP property removed.
29 Jun 2018: Global parameters moved to miri.parameters.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.

# Import global properties
from miri.parameters import READOUT_MODE
from miri.parameters import SUBARRAY as MIRI_SUBARRAY

from miri.simulators.find_simulator_file import find_simulator_file

#
# MIRI contains three focal plane modules, each of which contains a
# detector with a 1024 x 1024 zone illuminated by the instrument. Each
# focal plane module is identified by a unique Sensor Chip Assembly ID
# and a unique Focal Plane Module ID defined as follows. The SCA ID is
# used in MIRI documentation and the FPM ID is used in the FPM EIDP
# (see email from Tim Grundy, RAL):
#    * For the MIRI imager:
#         SCA 493 containing FPM S/N 106.
#    * For the long wavelength arm of the MIRI MRS:
#         SCA 494 containing FPM S/N 104.
#    * For the short wavelength arm of the MIRI MRS:
#         SCA 495 containing FPM S/N 105.
# Each 1024x1024 pixel detector has 4 extra non-illuminated reference
# columns just off the left and right edges of the illuminated zone.
# In addition, there is a separate bank of non-illuminated reference
# pixels ganged together known as reference outputs. These reference
# outputs are not contiguous with the illuminated zone, but the data
# is rearranged in level 1 FITS files so these reference outputs appear
# as extra rows on top of the normal detector image.
#
# Note that DARK_CURRENT_FILE describes how the dark current varies
# with detector temperature. A 2-D DARK_MAP describes how the dark current
# varies over the detector surface (including hot pixels which have
# excessive dark current). A 3 or 4-D DARK_MAP also describes how the
# dark current changes with group and integration.
#
# Note that detector drifts and latency have only been modelled in FAST
# mode. The parameters defined here are not valid in SLOW mode.
#
# The find_simulator_file function searches for a named file within a
# search path of simulator data files(starting with the current working
# directory) and returns the absolute path.
#
# TODO: Is there a better model to describe detector drifts?
_sca493 = {}
_sca493['SCA_ID'] = 493             # Numerical SCA ID
_sca493['FPM_ID'] = "FPMSN106"      # Unique FPM ID
_sca493['NAME'] = "Sensor Chip Assembly 493 with Focal Plane Module 106"
_sca493['DETECTOR'] = "MIRIMAGE"    # ASCII detector ID (previously "IM")
_sca493['CHIP'] = 'SiAs'            # Type of detector chip
_sca493['COMMENTS'] = "Describes MIRI FPM S/N 106 detector data with ref pixels"
_sca493['ILLUMINATED_ROWS'] = 1024
_sca493['ILLUMINATED_COLUMNS'] = 1024
_sca493['LEFT_COLUMNS'] = 4    # Reference columns on detector
_sca493['RIGHT_COLUMNS'] = 4   # Reference columns on detector
_sca493['BOTTOM_ROWS'] = 0     # There are no extra rows at the bottom
_sca493['TOP_ROWS'] = 256      # Reference rows in level 1 FITS image
_sca493['PIXEL_SIZE'] = 25.0   # Pixel size in microns
_sca493['THICKNESS'] = 470.0   # Detector thickness in microns
_sca493['WELL_DEPTH'] = 354720 # Well depth in electrons
_sca493['PEDESTAL'] = 10000    # Pedestal value in electrons
_sca493['BAD_PEDESTAL'] = 1000 # Pedestal value for bad pixels in electrons
_sca493['MEAN_GAIN'] = 5.5     # Mean gain (e/DN)  (Superceded by GAIN CDP)
_sca493['PERSISTENCE'] = 0.0   # Linear persistence factor (0.0 to 1.0)
# _sca493['PERSISTENCE'] = [1.0e-8, 0.03, 0.0]  # Persistence coefficients [2nd,1st,0th]
# NOTE: The following 4 parameters are valid for FAST mode only.
_sca493['LATENCY_SLOW'] = [1.67e-9, 136000.0] # Slow latency parameters [gain(1/e),decay]
_sca493['LATENCY_FAST'] = [0.002, 300.0]      # Fast latency parameters [gain,decay]
_sca493['ZP_SLOW'] = [45000.0, 0.0084]        # Slow zeropoint drift [const(e),scale(e/s)]
# Fast zeropoint jumps as a function of integration number and flux [[const(e),scale(s)]
_sca493['ZP_FAST'] = [[0.0, -2.917], [0.0, -2.292], [0.0, -2.396], [0.0, -2.408]]
# [1.0, 0.0] means the integrator is perfectly linear
# [1.0, -1.0] means the counter starts at 100% sensitivity and reduces to 0% at full well.
# _sca493['SENSITIVITY'] = [1.0, 0.0]  # Linearity sensitivity coeffs [const,slope]
_sca493['SENSITIVITY'] = [1.1, -0.4]  # Linearity sensitivity coeffs [const,slope]
_sca493['CLOCK_TIME'] = 1.0e-5 # Detector clock time in seconds
_sca493['RESET_WIDTH'] = 4 # The width of the reset pulse in clock cycles
_sca493['RESET_OVERHEAD'] = 3 # Number of clock cycles per reset
_sca493['FRAME_RESETS'] = 0    # Extra resets between integrations
_sca493['TARGET_TEMPERATURE'] = 6.7   # Target temperature in K
_sca493['DARK_CURRENT_FILE'] = find_simulator_file('dark_currentIM.fits')
_sca493['DARK_CURRENT'] = 0.21 # Nominal dark current level (electrons/s)
_sca493['QE_FILE'] = find_simulator_file('qe_measurementIM.fits')
_sca493['NOISEFACTOR'] = 1.0   # Noise adjustment factor

_sca494 = {}
_sca494['SCA_ID'] = 494             # Numerical SCA ID
_sca494['FPM_ID'] = "FPMSN104"      # Unique FPM ID
_sca494['NAME'] = "Sensor Chip Assembly 494 with Focal Plane Module 104"
_sca494['DETECTOR'] = "MIRIFULONG"  # ASCII detector ID (previously "LW")
_sca494['CHIP'] = 'SiAs'            # Type of detector chip
_sca494['COMMENTS'] = "Describes MIRI FPM S/N 104 detector data with ref pixels"
_sca494['ILLUMINATED_ROWS'] = 1024
_sca494['ILLUMINATED_COLUMNS'] = 1024
_sca494['LEFT_COLUMNS'] = 4    # Reference columns on detector
_sca494['RIGHT_COLUMNS'] = 4   # Reference columns on detector
_sca494['BOTTOM_ROWS'] = 0     # There are no extra rows at the bottom
_sca494['TOP_ROWS'] = 256      # Reference rows in level 1 FITS image
_sca494['PIXEL_SIZE'] = 25.0   # Pixel size in microns
_sca494['THICKNESS'] = 470.0   # Detector thickness in microns
_sca494['WELL_DEPTH'] = 359950 # Well depth in electrons
_sca494['PEDESTAL'] = 10000    # Pedestal value in electrons
_sca494['BAD_PEDESTAL'] = 1000 # Pedestal value for bad pixels in electrons
_sca494['MEAN_GAIN'] = 5.5     # Mean gain (e/DN) (Superceded by GAIN CDP)
_sca494['PERSISTENCE'] = 0.0   # Linear persistence factor (0.0 to 1.0)
# _sca494['PERSISTENCE'] = [1.0e-8, 0.03, 0.0]  # Persistence coefficients [2nd,1st,0th]
# NOTE: The following 4 parameters are valid for FAST mode only.
_sca494['LATENCY_SLOW'] = [1.67e-9, 136000.0] # Slow latency parameters [gain(1/e),decay(s)]
_sca494['LATENCY_FAST'] = [0.002, 300.0]      # Fast latency parameters [gain,decay(s)]
_sca494['ZP_SLOW'] = [45000.0, 0.0084]        # Slow zeropoint drift [const(e),scale(e/s)]
# Fast zeropoint jumps as a function of integration number and flux [[const(e),scale(s)]
_sca494['ZP_FAST'] = [[0.0, -2.917], [0.0, -2.292], [0.0, -2.396], [0.0, -2.408]]
# [1.0, 0.0] means the integrator is perfectly linear
# [1.0, -1.0] means the counter starts at 100% sensitivity and reduces to 0% at full well.
# _sca494['SENSITIVITY'] = [1.0, 0.0]  # Linearity sensitivity coeffs [const,slope]
_sca494['SENSITIVITY'] = [1.1, -0.4]  # Linearity sensitivity coeffs [const,slope]
_sca494['CLOCK_TIME'] = 1.0e-5 # Detector clock time in seconds
_sca494['RESET_WIDTH'] = 4 # The width of the reset pulse in clock cycles
_sca494['RESET_OVERHEAD'] = 3 # Number of clock cycles per reset
_sca494['FRAME_RESETS'] = 0    # Extra resets between integrations
_sca494['TARGET_TEMPERATURE'] = 6.7   # Target temperature in K
_sca494['DARK_CURRENT_FILE'] = find_simulator_file('dark_currentLW.fits')
_sca494['DARK_CURRENT'] = 0.21 # Nominal dark current level (electrons/s)
_sca494['QE_FILE'] = find_simulator_file('qe_measurementLW.fits')
_sca494['NOISEFACTOR'] = 1.0   # Noise adjustment factor

_sca495 = {}
_sca495['SCA_ID'] = 495             # Numerical SCA ID
_sca495['FPM_ID'] = "FPMSN105"      # Unique FPM ID
_sca495['NAME'] = "Sensor Chip Assembly 495 with Focal Plane Module 105"
_sca495['DETECTOR'] = "MIRIFUSHORT" # ASCII detector ID (previously "SW")
_sca495['CHIP'] = 'SiAs'            # Type of detector chip
_sca495['COMMENTS'] = "Describes MIRI FPM S/N 105 detector data with ref pixels"
_sca495['ILLUMINATED_ROWS'] = 1024
_sca495['ILLUMINATED_COLUMNS'] = 1024
_sca495['LEFT_COLUMNS'] = 4    # Reference columns on detector
_sca495['RIGHT_COLUMNS'] = 4   # Reference columns on detector
_sca495['BOTTOM_ROWS'] = 0     # There are no extra rows at the bottom
_sca495['TOP_ROWS'] = 256      # Reference rows in level 1 FITS image
_sca495['PIXEL_SIZE'] = 25.0   # Pixel size in microns
_sca495['THICKNESS'] = 470.0   # Detector thickness in microns
_sca495['WELL_DEPTH'] = 358190 # Well depth in electrons
_sca495['PEDESTAL'] = 10000    # Pedestal value in electrons
_sca495['BAD_PEDESTAL'] = 1000 # Pedestal value for bad pixels in electrons
_sca495['MEAN_GAIN'] = 5.5     # Mean gain (e/DN) (Superceded by GAIN CDP)
_sca495['PERSISTENCE'] = 0.0   # Linear persistence factor (0.0 to 1.0)
# _sca495['PERSISTENCE'] = [1.0e-8, 0.03, 0.0]  # Persistence coefficients [2nd,1st,0th]
# NOTE: The following 4 parameters are valid for FAST mode only.
_sca495['LATENCY_SLOW'] = [1.67e-9, 136000.0] # Slow latency parameters [gain(1/e),decay]
_sca495['LATENCY_FAST'] = [0.002, 300.0]      # Fast latency parameters [gain,decay]
_sca495['ZP_SLOW'] = [45000.0, 0.0084]        # Slow zeropoint drift [const(e),scale(e/s)]
# Fast zeropoint jumps as a function of integration number and flux [[const(e),scale(s)]
_sca495['ZP_FAST'] = [[0.0, -2.917], [0.0, -2.292], [0.0, -2.396], [0.0, -2.408]]
# [1.0, 0.0] means the integrator is perfectly linear
# [1.0, -1.0] means the counter starts at 100% sensitivity and reduces to 0% at full well.
# _sca495['SENSITIVITY'] = [1.0, 0.0]  # Linearity sensitivity coeffs [const,slope]
_sca495['SENSITIVITY'] = [1.1, -0.4]  # Linearity sensitivity coeffs [const,slope]
_sca495['CLOCK_TIME'] = 1.0e-5 # Detector clock time in seconds
_sca495['RESET_WIDTH'] = 4 # The width of the reset pulse in clock cycles
_sca495['RESET_OVERHEAD'] = 3 # Number of clock cycles per reset
_sca495['FRAME_RESETS'] = 0    # Extra resets between integrations
_sca495['TARGET_TEMPERATURE'] = 6.7   # Target temperature in K
# The pixel response function is the sensitivity variation across the
# surface of each individual pixel.
_sca495['PIXEL_RESPONSE'] = None      # No pixel response function
_sca495['DARK_CURRENT_FILE'] = find_simulator_file('dark_currentSW.fits')
_sca495['DARK_CURRENT'] = 0.21 # Nominal dark current level (electrons/s)
_sca495['QE_FILE'] = find_simulator_file('qe_measurementSW.fits')
_sca495['NOISEFACTOR'] = 1.0   # Noise adjustment factor

# Other detector descriptions (e.g. for other JWST instruments) could be
# added here.

# The following constants may be imported by SCA simulator software modules.
#
# Dictionary of known focal plane modules. The detector properties for
# each FPM can be obtained by looking up its unique detector name in this
# dictionary and using the result as another dictionary (i.e. the overall
# data structure is a dictionary of dictionaries).
#
DETECTORS_DICT = {'MIRIMAGE'    : _sca493,
                  'MIRIFULONG'  : _sca494,
                  'MIRIFUSHORT' : _sca495}

# The readout modes have been obtained from miri.parameters.
# This is the default mode.
DEFAULT_READOUT_MODE = 'FAST'


def flip_subarray_params(input_params):
    """
    
    Helper function to switch the row and column entries in a SUBARRAY tuple
    
    """
    assert isinstance(input_params, (tuple,list))
    assert len(input_params) == 4
    output_params = [input_params[1],
                     input_params[0],
                     input_params[3],
                     input_params[2]]
    return output_params

#
# Convert subarray parameters to the ordering needed by SCASim.
# TODO: Change SCASim to use the same ordering as the MIRI CDP software?
# Row and column numbers start at 1.
# The tuple contains (firstrow, firstcol, subrows, subcolumns)
#
SUBARRAY = {}
SUBARRAY['FULL'] = None
SUBARRAY['MASK1140'] =      flip_subarray_params( MIRI_SUBARRAY['MASK1140'] )
SUBARRAY['MASK1550'] =      flip_subarray_params( MIRI_SUBARRAY['MASK1550'] )
SUBARRAY['MASK1065'] =      flip_subarray_params( MIRI_SUBARRAY['MASK1065'] )
SUBARRAY['MASKLYOT'] =      flip_subarray_params( MIRI_SUBARRAY['MASKLYOT'] )
SUBARRAY['BRIGHTSKY'] =     flip_subarray_params( MIRI_SUBARRAY['BRIGHTSKY'] )
SUBARRAY['SUB256'] =        flip_subarray_params( MIRI_SUBARRAY['SUB256'] )
SUBARRAY['SUB128'] =        flip_subarray_params( MIRI_SUBARRAY['SUB128'] )
SUBARRAY['SUB64'] =         flip_subarray_params( MIRI_SUBARRAY['SUB64'] )
SUBARRAY['SLITLESSPRISM'] = flip_subarray_params( MIRI_SUBARRAY['SLITLESSPRISM'] )
# The following additional subarray options can be uncommented for testing
# special cases. The detectors are never actually read out using these modes.
# SUBARRAY['LRS1'] =          (   1,    1, 1024,  420)
# SUBARRAY['LRS2'] =          (   1,  292, 1024,  128)
# SUBARRAY['LRS3'] =          (   1,  292,  512,   64)
# SUBARRAY['AXIS256'] =       ( 384,  388,  256,  256)
# SUBARRAY['AXIS128'] =       ( 448,  452,  128,  128)
# SUBARRAY['AXIS64'] =        ( 480,  484,   64,   64)
# SUBARRAY['TEST64'] =        ( 128,  128,   64,   80)
# SUBARRAY['TEST32'] =        (   8,    8,   32,   32)
# SUBARRAY['RHS256'] =        (   1,  776,  256,  256)
STANDARD_SUBARRAYS = ('FULL', 'MASK1065', 'MASK1140', 'MASK1550', 'MASKLYOT',
                      'BRIGHTSKY', 'SUB256', 'SUB128', 'SUB64', 'SLITLESSPRISM')
DEFAULT_SUBARRAY = 'FULL'

#
# This fraction of cosmic ray events will cause charge leakage (negative jump)
# rather than a charge increase.
# NOTE: The negative jumps are more likely to be caused by a cosmic ray strike
# on the readout electronics rather than a true charge leakage.
#
COSMIC_RAY_LEAKAGE_FRACTION = 0.002

if __name__ == '__main__':
    print( "NOTE: The DetectorProperties module is supposed to be " \
        "imported by another module, not run as a main program." )
    print( "The following detector properties are defined:" )
    for detid in DETECTORS_DICT:
        print( "DETECTOR %s\n--------------------" % detid )
        detector = DETECTORS_DICT[detid]
        for key in detector:
            print( "%24s = %s" % (key, detector[key]) )
    print( "READOUT_MODE\n------------" )
    for key in READOUT_MODE:
        print( "%24s = %s" % (key, READOUT_MODE[key]) )
    print( "SUBARRAY\n--------" )
    for key in SUBARRAY:
        print( "%24s = %s" % (key, SUBARRAY[key]) )
    print( "Finished." )
