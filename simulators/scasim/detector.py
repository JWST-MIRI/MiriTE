#!/usr/bin/env python

"""

Module detector - Contains the DetectorArray class.

This module simulates a detector which accumulates charge and may be
read out non-destructively. The simulation is based on around an
integrator object (ImperfectIntegrator), which simulates an integrator
with Poisson noise and latency effects. The DetectorArray class adds a
bad pixel map and simulates the gain, flat-field and dark current.

The simulation is based on the parameters contained in the module
detector_properties.py and the calibration data contained in the MIRI
bad pixel mask, gain, nonlinearity, dark current and pixel flat-field
Calibration Data Products (CDPs).

:History:

18 Jun 2010: Created.
28 Jun 2010: Reworked to recommended Python syntax standards.
05 Jul 2010: Added verbose flag.
06 Jul 2010: Added the Amplifier class to simulate gain and read noise.
14 Jul 2010: Dark current can vary with temperature.
15 Jul 2010: Modified for use with the amplifier_properties module.
20 Jul 2010: set_readout_mode method added to forward nsample to
             Amplifier objects. Amplifiers created with a maxdn argument.
             Plot the dark current and read noise when in verbose mode.
27 Jul 2010: More accurate integration time calculation added.
02 Aug 2010: Added subarray mode.
06 Aug 2010: Documentation formatting problems corrected.
11 Aug 2010: Added FITS header generation.
             Corrected mistake in load_mv arguments.
16 Aug 2010: Modified to work with different focal plane modules.
24 Aug 2010: Read detector configuration data from FITS files.
27 Aug 2010: Added more checks. Corrected bug in exposure time calculation.
             Problem with read noise calculation solved by ensuring the
             readout array is floating point.
             Corrected bug in subarray positioning.
01 Sep 2010: Improved attribute checking for public functions.
             Reset the amplifiers during a reset().
03 Sep 2010: Added the possibility of a negative cosmic ray event.
             Added bad pixel map.
07 Sep 2010: Corrected bug where bad pixel mask not correctly sliced
             when there are no reference rows.
14 Sep 2010: Added flags to turn Poisson noise and read noise on and off.
07 Oct 2010: FITS header keywords updated after reading "Definition
             of the MIRI FM ILT Keywords" by T.W.Grundy.
11 Oct 2010: Added flags to turn dark current or amplifier effects on and off.
             Added makeplot flag.
12 Oct 2010: Number of frames and dropped frames included in exposure time
             calculation. Corrected the format of some doc strings.
19 Oct 2010: Ensure flux is stored in a floating point array.
             Take the effect of nframes on the Poisson noise and read noise
             into account.
05 Nov 2010: Number of frames per group taken into account when determining
             the effect of a cosmic ray hit.
12 Nov 2010: Dark current can now vary on a pixel to pixel basis, based on
             values contained in a dark current multipler map defined for
             each detector.
15 Nov 2010: MeasuredVariable class moved to miri.miritools.
16 Nov 2010: Some FITS header content rounded to a more sensible number
             of decimal places. Added function to set the random number
             generator seed, for testing.
18 Nov 2010: Added simulate_bad_pixels flag.
22 Nov 2010: Corrections to FPM and SCA IDs in FITS header.
             Pass the good pixel value to PoissonIntegrator rather than
             hard-coding it.
24 Nov 2010: Input of subarray data now supported.
26 Nov 2010: Subarray logic in readout corrected. Minor typo corrected.
26 Jan 2010: get_header function converted to set_metadata.
03 Mar 2011: Load dark map into DarkMap object rather than DataMap.
             Also check that the SCA_ID associated with a bad pixel map
             or dark map is the same as the detector.
04 Mar 2011: Documentation tweaks to resolve bad formatting.
22 Mar 2011: Plotting functions modified to use matplotlib figures
             and axes more flexibly (unfortunately the colorbar
             function doesn't work properly). MeasuredVariable class
             plotting functions changed.
25 Mar 2011: Persistence factor defined and passed to PoissonIntegrator.
21 Apr 2011: clock_time function added. Corrected TREFROW and the REFPIXEL
             and REFIMG groups of FITS keyword. Calculation of the number
             of reference output rows in subarray data corrected.
03 May 2011: Report SLITLESSPRISM frame time.
14 Jul 2011: Use of exceptions made more consistent and documented.
14 Sep 2011: Exposure data written in uint16 format rather than float32.
14 Sep 2011: Modified to keep up with the separation of
             miri.miritools.miridata into miri.miritools.metadata and
             miri.miritools.miricombdata
21 Sep 2011: Bad pixel masks stored as uint16 rather than uint8, to be
             compatible with DHAS flag values.
23 Sep 2011: Comment about persistence effects included in development
             and testing section.
29 Sep 2011: Added the ability to preflash the detector when created,
             to allow persistence effects to be checked.
05 Oct 2011: Persistence coefficients written to FITS header.
20 Oct 2011: Destructor added to remove large objects.
25 Oct 2011: detector_properties and amplifier_properties imported using
             ParameterFileManager. This allows new parameter files to be
             substituted by putting new versions into the current working
             directory, without the need to rebuild the source code.
29 Mar 2012: Removed unwanted imports.
02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
14 May 2012: MeasuredVariable class implemented as a JwstDataProduct.
21 May 2012: Plotting modified to use the miriplot utilities.
13 Nov 2012: Major restructuring of the package folder. Import statements
             updated.
15 Nov 2012: PoissonIntegrator class moved to simulators package.
21 May 2013: Use new MIRI data model instead of DataMap.
21 May 2013: Removed references to the old MiriMetadata object.
05 Jun 2013: Tidied up constructor and moved bad pixel mask and dark map
             access to separate functions. Rationalised attribute names.
13 Jun 2013: Changed metadata to match the new data model.
25 Oct 2013: Make the primary keyword for determining the detector ID
             'DETECTOR' rather than 'SCA_ID'. Pass detectorid as a
             parameter instead of scaid.
23 Apr 2014: Changed to use the ImperfectIntegrator class. Removed
             redundant methods.
07 May 2014: Extract the zeropoint drift, linearity and latency
             parameters from the detector properties.
08 May 2014: Charge trapping parameters extracted from detector
             properties.
05 Jun 2014: Removed charge trapping parameters and added slow and fast
             latency parameters. STILL UNDER TEST.
10 Jun 2014: Corrected some problems defining detector flux and dark current
             when there are no reference pixels.
17 Jun 2014: Improved memory management. Changed verbosity definitions
             and some log messages. Added ability to define the detector
             frame time explicitly when calculating the exposure time.
             Removed the old "preflash" latency implementation.
19 Jun 2014: Slow and fast zeropoint drift parameters obtained.
             Detector data converted to uint32 format rather than uint16
             format, to prevent data being truncated when
             simulate_gain=False.
01 Jul 2014: Corrected default detector temperature from 6.5K to 6.7K.
             This makes a significant difference to the read noise.
11 Jul 2014: Modified to use the new MiriMeasurement class instead of
             the old MeasuredVariable class to obtain the dark current
             and read noise.
18 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
02 Mar 2015: Include burst mode in subarray timing calculation.
21 May 2015: Added wait() method.
21 May 2015: Amplifier level removed and gain and read noise defined by
             calibration data products.
06 Aug 2015: Bug corrected in read noise calculation (the sampling factor
             was not taken into account). Noise calibration factor added.
             Gain and read noise metadata added to the output data.
14 Aug 2015: Inconsistency in name of bad pixel flag corrected.
03 Sep 2015: Legacy bad pixel mask format removed.
08 Sep 2015: Legacy dark map format removed. Made compatible with Python 3.
02 Oct 2015: Attribute names modified to be compatible with documentation.
07 Oct 2015: Made exception catching Python 3 compatible.
10 Nov 2015: Tightened up on some of the data type conversion.
04 Dec 2015: Added a logger.
10 Mar 2016: Bad pixel simulation moved from integrator to here. Bad pixels
             are no longer zeroed but respond to the dark current. Pedestal
             values are used to keep DN readings away from 0 and (optionally)
             to distinguish bad pixels.
23 Mar 2016: Obtain MASK, PIXELFLAT and GAIN calibrations from CDP files.
             NOTE: This generates more error messages, because CDP files
             are not available for all combinations of readout modes.
30 Mar 2015: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
26 Apr 2016: Brought back the SCA_ID metadata keyword.
29 Apr 2016: Added path and credentials for read-only access to simulator-
             specific CDP files.
03 May 2016: Added reference output array size and shape calcuations.
05 May 2016: Look for alternatives if a flat-field matching the exact criteria
             cannot be found. Linearity renamed to sensitivity.
             simulate_amp_effects split into simulate_gain and
             simulate_nonlinearity. Shortened metadata comments.
06 May 2016: Use the pixel flat-field appropriate for each subarray mode
             (and filter or band, if appropriate). Optionally, plot the
             CDPs used.
06 Jun 2016: Added flags to allow the version numbers of imported CDPs to be
             specified. Protect against a null readnoise map or gain map.
20 Jun 2016: Missing flat-field logged as a warning instead of an error.
             It was messing up the unit test output.
05 Jul 2016: CDP-6 dark maps added. Integration number added to dark map
             selection. Dark maps are now used unnormalised, and the dark
             current vs temperature measurement is instead normalised at 6.7K.
14 Jul 2016: Added simulate_drifts and simulate_latency flags. Distinguish
             between dead and noisy bad pixels and flag reference pixels
             as "non-science".
05 Aug 2016: Reset the reference output array size when the subarray mode is
             FULL.
05 Sep 2016: Separated add_calibration_data from the constructor.
28 Sep 2016: miri.miritools.dataproduct renamed miri.datamodels and
             miri.miritools renamed miri.tools.
10 Jan 2017: Report the filter and band when creating a new detector object.
11 Jan 2017: When searching for a flat-field with no particular filter,
             first try to find a CDP which works for all filters before
             accepting one designed for a specific filter. If the 'P750L'
             filter is required, only accept a CDP for that particular filter.
14 Feb 2017: Log the name of the DARK calibration file used.
13 Jun 2017: Added cdp_ftp_path parameter.
16 Jun 2017: Added functions to obtain an averaged DARK CDP file, or the
             full DARK CDP file, from the CDP repository.
             Corrected mistake in the setting of ftp_path.
20 Jun 2017: Mean gain is (e/DN), not (DN/e). Multiply the dark current
             and the read noise by the gain, rather than dividing.
             Scale the DARK to the expected dark current level.
             Scale the read noise to compensate for noise already added
             by the DARK.
07 Jul 2017: Added alternative add_dark function. Added INSTALLED_DARK and
             AVERAGED_DARK flags to control which dark is used.
18 Jul 2017: Corrected a problem handling subarray flat-fields.
21 Jul 2017: Do not issue log messages when scaling something by 1.0.
26 Jul 2017: Corrected a problem in add_dark_map_cdp caused by an incorrect
             MiriCDPInterface API.
27 Jul 2017: Count the number of cosmic ray events and the number of pixels
             affected. Added an option to the wait() method to simulate a
             non-zero background flux while the detector is idling.
             Make some metadata optional.
04 Sep 2017: Explicitly delete old calibration data objects to save memory.
             Added _clear_calibration_data() method.
13 Oct 2017: New frame time calculation from Mike Ressler.
             SLOW mode now uses 8 out of 9 samples. READOUT_MODE now defines
             samplesum, sampleskip and refpixsampleskip parameters separately.
             Added frame_rti function and associated tests.
13 Dec 2017: INSTALLED_DARK option removed.
14 Feb 2018: Added table-based non-linearity correction, which can be selected
             with the NONLINEARITY_BY_TABLE flag. Added function to read the
             non-linearity CDP.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import os
import math
import numpy as np

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# MIRI data models
from miri.datamodels.miri_measurement import MiriMeasurement
from miri.datamodels.cdp import MiriDarkReferenceModel
#from miri.datamodels import MiriMeasuredModel
from miri.datamodels.cdp import MiriGainModel, MiriReadnoiseModel, \
    MiriBadPixelMaskModel
from miri.datamodels.cdplib import get_cdp, cdp_version_decode, \
    MiriCDPInterface, MIRI_SUBARRAYS

from miri.simulators import ImperfectIntegrator

# Search for the detector parameters file and parse it into a
# properties dictionary. The file is searched for in 3 places:
# (a) The current directory.
# (b) The directory where this Python file is being executed from.
# (c) The miri.simulators.scasim installation directory.
from miri.tools.filesearching import ParameterFileManager, make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
detector_properties = ParameterFileManager(
                            "detector_properties.py",
                            search_path=search_path,
                            description="detector properties",
                            logger=LOGGER)

# Initialise some CDP FTP defaults
SIM_CDP_FTP_PATH   = 'CDPSIM'
SIM_CDP_FTP_USER   = 'cdpuser'
SIM_CDP_FTP_PASSWD = 'R7ZWEXEEsAH7'

# Bad pixel mask flags
# TODO: Get from MiriBadPixelMask data model
MASK_DO_NOT_USE = 1
MASK_DEAD = 2
MASK_HOT = 4
MASK_UNRELIABLE_SLOPE = 8
MASK_RC_PIXEL = 16
MASK_NON_SCIENCE = 512

# Configure the simulation of the DARK.
# Set AVERAGED_DARK True to simulate an average 2-D DARK scaled to the expected
# dark current during each integration.
# Set AVERAGED_DARK False to add the full 4-D DARK at the end of each exposure.
AVERAGED_DARK = False

# Configure the nonlinearity simulation.
# Setting this to True simulates nonlinearity using a lookup table
# obtained from the nonlinearity CDP. The wavelength dependence of
# the effect will be included but the saturation level will be more
# uncertain.
# Setting this to False simulates nonlinearity by simulating the
# change in detector sensitivity with charge accumulated. The wavelength
# dependence will not be included, but the saturation level will be
# more accurate.
NONLINEARITY_BY_TABLE = True

# Set to True to include debugging information in the metadata
EXTRA_METADATA = False

#
# Global helper functions
#
def frame_rti(detectorrows, ampcolumns, colstart, colstop, rowstart, rowstop,
              sampleskip, refpixsampleskip, samplesum, resetwidth,
              resetoverhead, burst_mode):
    """
    
    Helper function which calculates the number of clock cycles, or RTI,
    needed to read out the detector with a given a list of FPGA parameters.
    This function can be used to test the frame time calculation over a
    wide range of scenarios.
    
    :Reference:
    
    JPL MIRI DFM 478 04.02, MIRI FPS Exposure Time Calculations,
    M. E. Ressler, October 2014

    :Parameters:
    
    detectorrows:
        Total number of detector rows.
    ampcolumns
        Total number of amplifier columns (including reference columns).
        = Total number of detector columns / Number of amplifiers.
    colstart: int
        Amplifier column at which readout starts.
        Derived from the subarray mode.
    colstop: int
        Amplifier column at which readout stops.
        Derived from the subarray mode.
    rowstart: int
        Detector row at which readout starts.
        Derived from the subarray mode.
    rowstop: int
        Detector row at which readout stops.
        Derived from the subarray mode.
    sampleskip: int
        Number of clock cycles to dwell on a pixel before reading it.
        Derived from the readout mode.
    refpixsampleskip: int
        Number of clock cycles to dwell on a reference pixel before reading it.
        Derived from the readout mode.
    samplesum: int
        Number of clock cycles to sample a pixel during readout.
    resetwidth: int
        Width of the reset pulse in clock cycles.
        Derived from the readout mode.
    resetoverhead: int
        The overhead, in clock cycles, associated with setting shift
        registers back to zero.
    burst_mode: boolean, optional, default is False
        True if a subarray is being read out in burst mode,
        which skips quickly over unwanted columns.
        
    :Returns:
    
    frame_rti: float
        Number of clock cycles (RTI) for readout.
    
    """            
    # The time taken to skip a row not contained in a subarray is the
    # row reset overhead, which is the reset pulse width plus the overhead
    # to reset the shift registers.
    row_reset = resetwidth + resetoverhead
    rti_row_non_roi = row_reset

    # The number of clock cycles needed to read a pixel is the         
    pix_clocks = sampleskip + samplesum
    refpix_clocks = refpixsampleskip + samplesum
    
    if (colstart == 1) or (not burst_mode):
        # For a row within a subarray (if burst_mode is False or the subarray
        # begins at column 1), we start with the row reset overheads, then add
        # the reference pixel time, refpix_clocks, then add the science pixel
        # time, n * pix_clocks.
        rti_row_roi = row_reset
#         strg = "rti_row_roi=%d " % row_reset
        if colstart == 1:
            rti_row_roi += refpix_clocks
#             strg += " + %d " % refpix_clocks
        else:
            rti_row_roi += resetwidth
        if colstop == ampcolumns:
            rti_row_roi += (colstop-2) * pix_clocks
#             strg += " + (%d * %d) " % ((colstop-2), pix_clocks)
            rti_row_roi += refpix_clocks
#             strg += " + %d " % refpix_clocks
        else:
            rti_row_roi += (colstop-1) * pix_clocks
#             strg += " + (%d * %d) " % ((colstop-1), pix_clocks)
#             print(strg)
    else:
        # In burst_mode, when a subarray does not touch the left hand edge,
        # each row begins with the same reset overhead, then the left hand
        # pixels are clocked at 5 times the normal rate until the subarray
        # is reached. If the number of burst columns is not an even multiple
        # of 5, the readout is paused until a normal clock boundary.
        rti_row_roi = row_reset
#         strg = "rti_row_roi=%d " % row_reset
        if colstart == 1:
            rti_row_roi += refpix_clocks
#             strg += " + %d " % refpix_clocks
        else:
            rti_row_roi += refpix_clocks + ((colstart-2)*pix_clocks+4)//5
#             strg += " + %d + %d " % (refpix_clocks,
#                                      ((colstart-2)*pix_clocks+4)//5)
            
        if colstop == ampcolumns:
            rti_row_roi += ((colstop-2)-colstart+1) * pix_clocks
#             strg += " + (%d * %d) " % (((colstop-2)-colstart+1), pix_clocks)
            rti_row_roi += refpix_clocks
            strg += " + %d " % refpix_clocks
        else:
            rti_row_roi += (colstop-colstart+1) * pix_clocks
#             strg += " + (%d * %d) " % ((colstop-colstart+1), pix_clocks)
#             print(strg)

    # The total number of clock cycles to read the frame is the sum
    # of the portion before the ROI, after the ROI and during the ROI.
#         print("rti_row_non_roi=", rti_row_non_roi, "rti_row_roi=", rti_row_roi)
    frame_rti = (rowstart-1) * rti_row_non_roi + \
        (rowstop-rowstart+1) * rti_row_roi + \
        (detectorrows-rowstop) * rti_row_non_roi
    return frame_rti


class DetectorArray(object):
    """
    
    Class DetectorArray - Simulates the behaviour of a MIRI detector which
    consists of an illuminated central zone with extra dark columns at the
    left and right edges of the detector. Note that the total size of the
    detector itself in pixels is:
    
        (left_columns + columns + right_columns columns) x (rows)
    
    In addition, this class simulates the rearrangement of the data into
    the image format expected within a level 1 FITS file, so the class
    also supports bottom_rows and top_rows additional rows of reference
    pixel data at the bottom and/or top of the detector.
    (For MIRI the reference rows are always at the top, but bottom_rows
    allows some flexibility for reuse with future detectors.)
                 
    :Parameters:
    
    detectorid: string
        The detector ID, identifying a particular detector.
        The MIRI instrument has three detectors: 'MIRIMAGE', 'MIRIFULONG'
        and 'MIRIFUSHORT'.
        NOTE: The 'IC' and 'JPL' detectors are not currently recognised.
    rows: int
        The number of illuminated detector rows.
    columns: int
        The number of illuminated detector columns.
    temperature: float
        The detector temperature in K. This is used to look up the
        dark current and readout noise multipliers.
    left_columns: int, optional, default=4
        The number of extra dark columns at the left edge of the detector.
    right_columns: int, optional, default=4
        The number of extra dark columns at the right edge of the detector.
    bottom_rows: int, optional, default=0
        The number of extra reference rows at the bottom edge of the
        detector.
    top_rows: int, optional, default=256
        The number of extra reference rows at the top edge of the detector.
    well_depth: int, optional, default=250000
        The maximum number of particles that can be contained in each
        detector pixel.
    particle: string, optional, default="electron"
        The name of the particles being counted by the detector pixels.
    time_unit: string, optional, default="seconds"
        The unit of time measurement.
    readpatt: string, optional, default='FAST'
        A default readout mode with which to search for CDP files.
    subarray: string, optional, default='FULL'
        Detector subarray mode for output. This can be one 'FULL',
        'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
        'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
        full-frame and the other modes read out portions of the detector
        as described in the MIRI Operational Concept Document.
        The simulator can also accept the other subarray modes defined
        in detector_properties.py for test purposes.
    mirifilter: string, optional, default=None
        The name of the MIRI imager filter associated with the data (if any).
        This information is only used to look up an appropriate flat-field.
    miriband: string, optional, default=None
        The name of the MRS band associated with the data (if any).
        This information is only used to look up an appropriate flat-field.
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for example
        to observe what effects in a simulation are caused by Poisson
        noise).
    simulate_read_noise: boolean, optional, default=True
        A flag that may be used to switch off read noise (for example
        to observe what effects in a simulation are caused by read
        noise).
    simulate_bad_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of bad
        pixels in the data, even if a bad pixel map containing bad pixels
        is specified in the detector properties.
    simulate_dark_current: boolean, optional, default=True
        A flag that may be used to switch off the addition of dark current
        (for example to observe what effects in a simulation are caused by
        dark current).
    simulate_flat_field: boolean, optional, default=True
        A flag that may be used to switch off the simulation of the pixel
        flat-field (for example to observe the effect of quantum efficiency
        in isolation).
    simulate_gain: boolean, optional, default=True
        A flag that may be used to switch off the bias and gain effects.
        *Note that when this flag is False the ratio of DNs to electrons
        is exactly 1.0*
    simulate_nonlinearity: boolean, optional, default=True
        A flag that may be used to switch off non-linearity effects
        (for example to experiment with jump detection on perfectly
        linear data).
    simulate_drifts: boolean, optional, default=True
        A flag that may be used to switch off the simulation of
        detector drifts, such as the zeropoint drift.
    simulate_latency: boolean, optional, default=True
        A flag that may be used to switch off the simulation of
        detector latency and persistence effects.
    cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
        If specified, a list of folders (or folders) on the SFTP host
        where the MIRI CDPs are held to be searched, consisting of a
        list of folder names separated by a ":" delimiter.
        Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
        If not specified, the default CDP repository at Leuven is used.
    readnoise_version: string, optional, default=''
        A specific readnoise CDP version number of the form 'x.y.z'.
    bad_pixels_version: string, optional, default=''
        A specific bad pixel mask CDP version number of the form 'x.y.z'.
    flat_field_version: string, optional, default=''
        A specific pixel flat-field CDP version number of the form 'x.y.z'.
    linearity_version: string, optional, default=''
        A specific nonlinearity CDP version number of the form 'x.y.z'.
    gain_version: string, optional, default=''
        A specific gain CDP version number of the form 'x.y.z'.
    makeplot: boolean, optional, default=False
        Plotting flag. Activates plotting of data when True.
    verbose: int, optional, default=2
        Verbosity level. Activates logging statements when non-zero.
        
        * 0 - no output at all except for error messages
        * 1 - warnings and errors only
        * 2 - normal output
        * 3, 4, 5 - additional commentary
        * 6, 7, 8 - extra debugging information

    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
        
    :Raises:
    
    ValueError
        Raised if any of the initialisation parameters are out of range.
    TypeError
        Raised if any of the initialisation parameters are of the
        wrong type, size or shape.
    KeyError
        Raised if any of the initialisation parameters do not contain
        a recognised keyword.
        
    """

    def __init__(self, detectorid, rows, columns, temperature,
                 left_columns=4, right_columns=4,
                 bottom_rows=0, top_rows=256, well_depth=250000,
                 particle="electron", time_unit="seconds", readpatt='FAST',
                 subarray='FULL', mirifilter=None, miriband=None,
                 simulate_poisson_noise=True, simulate_read_noise=True,
                 simulate_bad_pixels=True,
                 simulate_dark_current=True, simulate_flat_field=True,
                 simulate_gain=True, simulate_nonlinearity=True,
                 simulate_drifts=True, simulate_latency=True,
                 cdp_ftp_path=SIM_CDP_FTP_PATH,
                 readnoise_version='', bad_pixels_version='',
                 dark_map_version='', flat_field_version='',
                 linearity_version='', gain_version='',
                 makeplot=False, verbose=2, logger=LOGGER):
        """
        
        Constructor for class DetectorArray.
        
        Parameters: See class doc string.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("scasim.detector")
        # Note that all parameters are explicitly converted into the
        # expected data type, since some of the values extracted from
        # the properties dictionary can sometimes be converted by Python
        # into strings, which then upsets formatted I/O statements.
        self._makeplot = bool(makeplot)
        self._verbose = int(verbose)
        if verbose > 3:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug( "+++DetectorArray object created with" + \
                " detectorid=" + str(detectorid) + \
                " rows=" + str(rows) + \
                " columns=" + str(columns) + \
                " left_columns=" + str(left_columns) + \
                " right_columns=" + str(right_columns) + \
                " bottom_rows=" + str(bottom_rows) + \
                " top_rows=" + str(top_rows) + \
                " well_depth=" + str(well_depth) )
         
        # Extract the properties of the particular sensor chip assembly
        # being simulated. (Convert to a string in case the caller
        # supplied an integer. Only a string will work as a dictionary
        # keyword.)
        self.detectorid = detectorid
        
        try:
            self._sca = detector_properties.get('DETECTORS_DICT',
                                                self.detectorid)
        except KeyError:
            strg = "%s is not a known detector ID." % self.detectorid
            raise KeyError(strg)

        self.fpm_id   = self._sca['FPM_ID']
        self.fpm_name = self._sca['NAME']
        self.sca_id   = int(self._sca['SCA_ID'])
        self.detector_name = self._sca['DETECTOR']
        self.chip     = self._sca['CHIP']
        
        # The detector must have a sensible size and reference pixel layout.
        if int(rows) <= 0 or int(columns) <= 0:
            strg = "The detector must have a non-zero size."
            raise ValueError(strg)
        
        # The detector itself has extra unilluminated columns at the left and
        # right sides. A ImperfectIntegrator object is created to include both
        # kinds of pixels.
        if int(left_columns) < 0 or int(right_columns) < 0 or \
           int(bottom_rows) < 0 or int(top_rows) < 0:
            strg = "The number of reference pixels cannot be negative."
            raise ValueError(strg)
        
        self.illuminated_shape = (int(rows), int(columns))
        self.left_columns = int(left_columns)
        self.right_columns = int(right_columns)
        self.bottom_rows = int(bottom_rows)
        self.top_rows = int(top_rows)
        self.subrows_bottom = self.bottom_rows
        self.subrows_top = self.top_rows
        
        detector_columns = left_columns + columns + right_columns
        detector_rows = bottom_rows + rows + top_rows
        self.detector_shape = (detector_rows, detector_columns)
        min_illuminated_row = bottom_rows
        max_illuminated_row = bottom_rows + rows - 1
        
        # Determine the shape of the isolated reference output data
        self._normal_amps = 4
        refout_rows = rows
        refout_columns = top_rows * detector_columns / refout_rows
        self.refout_shape = (refout_rows, refout_columns)
        
        # Default readout mode
        self.samplesum = 1
        self.sampleskip = 0
        self.refpixsampleskip = 3
        self.nframes = 1


# NOTE: AMPLIFIER CLASS IS STILL AVAILABLE BUT NO LONGER USED FOR MIRI.
#         # Add the readout amplifiers associated with this detector
#         self._add_amplifiers(temperature, min_illuminated_row,
#                              max_illuminated_row,
#                              simulate_read_noise=simulate_read_noise,
#                              simulate_gain=simulate_gain)
     
        # Create a particle counter with the specified properties.
        self.particle = particle
        self.time_unit = time_unit
        self.well_depth = int(well_depth)
        self.pixels = ImperfectIntegrator(
                                detector_rows, detector_columns,
                                particle=particle,
                                time_unit=time_unit,
                                bucket_size=well_depth,
                                simulate_poisson_noise=simulate_poisson_noise,
                                verbose=verbose)
        self.temperature = temperature
        # Initialise the cosmic ray counters.
        self.cosmic_ray_count = 0
        self.cosmic_ray_pixel_count = 0
        
        self.simulate_drifts = simulate_drifts
        if self.simulate_drifts:
            self.pixels.set_zeropoint(self._sca['ZP_SLOW'], self._sca['ZP_FAST'])
        else:
            self.pixels.set_zeropoint(None, None)
        
        self.simulate_nonlinearity = simulate_nonlinearity
        if self.simulate_nonlinearity and not NONLINEARITY_BY_TABLE:
            self.pixels.set_sensitivity( self._sca['SENSITIVITY'] )
        else:
            self.pixels.set_sensitivity( [1.0, 0.0] )
        
        self.simulate_latency = simulate_latency
        if self.simulate_latency:
            self.pixels.set_persistence(self._sca['PERSISTENCE'])
            self.pixels.set_latency(self._sca['LATENCY_SLOW'],
                                    self._sca['LATENCY_FAST'])
        else:
            self.pixels.set_persistence(None)
            self.pixels.set_latency(None, None)

        # Add the calibration data associated with this detector
        self.simulate_bad_pixels = simulate_bad_pixels
        self.simulate_dark_current = simulate_dark_current
        self.simulate_gain = simulate_gain
        self.mean_gain = 1.0
        self.mean_dark = 0.0
        self.simulate_read_noise = simulate_read_noise
        self.simulate_flat_field = simulate_flat_field
        self.add_calibration_data(self._sca['DETECTOR'],
                            readpatt=readpatt, subarray=subarray,
                            mirifilter=mirifilter, miriband=miriband,
                            cdp_ftp_path=cdp_ftp_path,
                            bad_pixels_version=bad_pixels_version,
                            dark_map_version=dark_map_version,
                            flat_field_version=flat_field_version,
                            linearity_version=linearity_version,
                            readnoise_version=readnoise_version,
                            gain_version=gain_version)
            
    def __del__(self):
        """
        
        Destructor for class DetectorArray.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        try:
            # Objects created in the constructor
            del self.pixels
            self._clear_calibration_data()
        except Exception:
            pass

# THIS CODE ALLOWS AMPLIFIERS TO BE TREATED SEPARATELY
#     def _add_amplifiers(self, temperature,
#                         min_illuminated_row, max_illuminated_row,
#                         simulate_read_noise=True, simulate_gain=True,
#                         simulate_nonlinearity=True):
#         """
#         
#         Add the list of amplifiers associated with the detector.
#         Helper function used during __init__().
#         
#         """
#         # Create a list of amplifiers responsible for reading out
#         # particular subsets of the rows and columns of the detector.
#         # Each amplifier has its own gain and read noise.
#         #
#         # Amplifier 1 - Every 4th column starting at column 0
#         #               for the illuminated rows only.
#         # Amplifier 2 - Every 4th column starting at column 1
#         #               for the illuminated rows only.
#         # Amplifier 3 - Every 4th column starting at column 2
#         #               for the illuminated rows only.
#         # Amplifier 4 - Every 4th column starting at column 3
#         #               for the illuminated rows only.
#         self.amplifier_list = []
#         self._normal_amps = 0
#         self._ref_amps = 0
#         namp = 0
#         
#         try:
#             amplist = amplifier_properties.get('AMPLIFIERS_DICT',
#                                                self.detectorid)
#         except KeyError:
#             strg = "Detector %s is " % self.detectorid
#             strg += "not included in amplifier properties."
#             raise KeyError(strg)
# 
#         total_amps = len(amplist)
#         
#         # If required, create a new plotting window to contain the read
#         # noise plot.
#         if self._makeplot:
#             fig = mplt.new_figure(1, figsize=(10,8),
#                                   stitle="Read noise measurements")
# 
#         for amp in amplist:
#             namp += 1
#             noisemv = MiriMeasurement(amp['READ_NOISE_FILE'], name="Read noise")
#             column = amp['READ_NOISE_COL'] - 1
#             read_noise = noisemv.lookup(temperature, column=column)
# 
#             # Plot the read noise if required, using a different
#             # sub-plot for each amplifier.
#             if self._makeplot:
#                 pltstr = "for %s" % amp['NAME']
#                 ax = mplt.add_subplot(plotfig=fig, subrows=2, subcols=3,
#                                       subno=namp)
#                 ax = noisemv.plot(plotfig=fig, plotaxis=ax, columns=column,
#                                   title=pltstr)
#                 if namp == total_amps:
#                     strg = "Plotting read noise measurements. " \
#                         "Close plot window to continue."
#                     mplt.show_plot(prompt=strg)
# 
#             if amp['TYPE'] == "illuminated":
#                 self._normal_amps += 1
#                 colstart = int(amp['REGION'])
#                 colstep = \
#                     int(amplifier_properties.get('N_ILLUM_AMPLIFIERS', 
#                                                  self.detectorid))
#                 ampobject = Amplifier(
#                                 namp,
#                                 min_illuminated_row, max_illuminated_row,
#                                 None, colstart, None, colstep,
#                                 amp['BIAS'], amp['GAIN'],
#                                 read_noise, amp['MAX_DN'],
#                                 gaintype=amp['GAINTYPE'],
#                                 simulate_read_noise=simulate_read_noise,
#                                 simulate_gain=simulate_gain,
#                                 simulate_nonlinearity=simulate_nonlinearity,
#                                 verbose=self._verbose)
# 
#                 self.amplifier_list.append(ampobject)
#             else:
#                 self._ref_amps += 1
#                 if amp['REGION'] == "bottom":
#                     max_reference_row = min_illuminated_row - 1
#                     ampobject = Amplifier(
#                                     namp, None, max_reference_row, None,
#                                     None, None, None,
#                                     amp['BIAS'], amp['GAIN'],
#                                     read_noise, amp['MAX_DN'],
#                                     gaintype=amp['GAINTYPE'],
#                                     simulate_read_noise=simulate_read_noise,
#                                     simulate_gain=simulate_gain,
#                                     simulate_nonlinearity=simulate_nonlinearity,
#                                     verbose=self._verbose)
#                     self.amplifier_list.append(ampobject)
#                 else:
#                     # Assume 'top'
#                     min_reference_row = max_illuminated_row
#                     ampobject = Amplifier(
#                                     namp, min_reference_row, None, None,
#                                     None, None, None,
#                                     amp['BIAS'], amp['GAIN'],
#                                     read_noise, amp['MAX_DN'],
#                                     gaintype=amp['GAINTYPE'],
#                                     simulate_read_noise=simulate_read_noise,
#                                     simulate_gain=simulate_gain,
#                                     simulate_nonlinearity=simulate_nonlinearity,
#                                     verbose=self._verbose)
#                     self.amplifier_list.append(ampobject)
#                     
#         # Set the initial electronic bias.
#         set_reference_level(float(amplifier_properties['INITIAL_REF_LEVEL']))

    def add_calibration_data(self, detector, readpatt=None, subarray=None,
                    mirifilter=None, miriband=None,
                    cdp_ftp_path=SIM_CDP_FTP_PATH,
                    bad_pixels_version='', dark_map_version='',
                    flat_field_version='', linearity_version='',
                    readnoise_version='', gain_version=''):
        """
        
        Add calibration data which defines the behaviour of the detector.
        The following calibration data are recognised:
        
        * Bad pixel mask
        * Dark map
        * Flat-field map
        * Linearity correction
        * Readnoise map
        * Gain map
        
        :Parameters:
        
        detector: string
            The detector name with which to look up the pixel flat-field.
        readpatt: string, optional
            The readout mode with which to look up the flat-field.
        subarray: string, optional
            The subarray with which to look up the flat-field
            or the dark calibration.
        mirifilter: string, optional
            The filter name with which to look up a flat-field
        miriband: string, optional
            The band name with which to look up a flat-field
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        bad_pixels_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        dark_map_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        flat_field_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        linearity_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        readnoise_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        gain_version: string, optional, default=''
            A specific version number of the form 'x.y.z'.
        
        """
        # Get the bad pixel mask associated with this detector.
        self.bad_pixels = None
        if self.simulate_bad_pixels:
            self.add_bad_pixel_mask(self._sca['DETECTOR'],
                                    cdp_ftp_path=cdp_ftp_path,
                                    cdp_version=bad_pixels_version)
            # Distinguish the good and bad pixels by adding a different
            # fixed pedestal value.
            pedestal = np.ones_like(self.pixels.expected_count) * \
                        self._sca['PEDESTAL']
            if self.bad_pixels is not None:
                bad_zones = np.where((self.bad_pixels & MASK_DO_NOT_USE) > 0)
                pedestal[bad_zones] = self._sca['BAD_PEDESTAL']
            self.pixels.set_pedestal( pedestal )
        else:
            self.add_bad_pixel_mask(None)
            # Add a fixed good pixel pedestal to all pixels.
            pedestal = np.ones_like(self.pixels.expected_count) * \
                        self._sca['PEDESTAL']
            self.pixels.set_pedestal( pedestal )
        
        # Obtain the amplifier gain for this detector.
        self.gain_map = None
        if self.simulate_gain:
            self.add_gain_map(self._sca['DETECTOR'], cdp_ftp_path=cdp_ftp_path,
                              cdp_version=gain_version)
        else:
            self.add_gain_map(None)
            
        # Extract the bias and dark current from the DARK CDP associated
        # with this detector.
        self.dark_map = None
        if self.simulate_dark_current:
            # Read the DARK from the MIRI CDP repository.
            # The AVERAGED_DARK flag determines whether an averaged dark or
            # the full 4-D dark is obtained.
            self.add_dark_map_cdp(self._sca['DETECTOR'], readpatt=readpatt,
                              subarray=subarray, averaged=AVERAGED_DARK,
                              cdp_ftp_path=cdp_ftp_path,
                              cdp_version=dark_map_version)
        else:
            # Do not simulate the dark current.
            self.add_dark_map_cdp(None)
        
        # Look up the dark current multiplier for the given detector temperature.
        darkmv = MiriMeasurement(init=self._sca['DARK_CURRENT_FILE'],
                                 name="Dark current multiplier", mcolumns=1,
                                 interptype='LOGLIN')
        # Plot the dark current if requested.
        if self._makeplot:
            pltstr = "Dark current multipler for %s" % self._sca['NAME']
            darkmv.plot(title=pltstr)

        self.dark_mv = darkmv
        if self.simulate_dark_current:
            self.dark_current = darkmv.lookup(self.temperature)
        else:
            self.dark_current = 0.0
            
        # Get the pixel flat-field associated with this detector.
        self.flat_map = None
        if self.simulate_flat_field:
            self.add_flat_map(self._sca['DETECTOR'], readpatt=readpatt,
                              subarray=subarray, mirifilter=mirifilter,
                              miriband=miriband, cdp_ftp_path=cdp_ftp_path,
                              cdp_version=flat_field_version)
        else:
            self.add_flat_map(None)

        # Get the pixel flat-field associated with this detector.
        self.linearity_table_left = None
        self.linearity_table_right = None
        if self.simulate_nonlinearity and NONLINEARITY_BY_TABLE:
            self.add_linearity_table(self._sca['DETECTOR'], mirifilter=mirifilter,
                              miriband=miriband, cdp_ftp_path=cdp_ftp_path,
                              cdp_version=linearity_version)
        else:
            self.add_linearity_table(None)
 
        # Get the read noise associated with this detector.
        self.readnoise_map = None
        if self.simulate_read_noise:
            self.add_readnoise_map(self._sca['DETECTOR'], readpatt=readpatt,
                                   cdp_ftp_path=cdp_ftp_path,
                                   cdp_version=readnoise_version)
            # Start with the noise calibration factor contained in the
            # detector properties.
            self.noise_factor = self._sca['NOISEFACTOR']
            # Scale the read noise to avoid double-counting the noise when
            # adding a DARK map or multiplying by a FLAT which already
            # contains some read noise.
            if self.readnoise_map is not None:
                if self.simulate_dark_current and \
                   self.dark_averaged and (self.dark_map is not None):
                    # Scale the read noise to remove the contribution
                    # already contained in the DARK.
                    posdark = np.where(self.dark_map > 0.0)
                    dknoise = np.std(self.dark_map[posdark])
                    rdnoise = np.mean(self.readnoise_map)
                    nsquared = dknoise / rdnoise
                    if nsquared >= 1.0:
                        self.noise_factor = 0.0
                    else:
                        self.noise_factor = self.noise_factor * \
                            math.sqrt(1.0 - nsquared)
        else:
            self.read_noise = 0.0
            self.noise_factor = 1.0

    def _clear_calibration_data(self):
        """
        
        Clear the calibration data created by add_calibration_data.
                
        :Parameters:
        
        None
        
        """
        if self.bad_pixels is not None:
            del self.bad_pixels
        self.bad_pixels = None
        
        if self.gain_map is not None:
            del self.gain_map
        self.gain_map = None
        
        if self.dark_mv is not None:
            del self.dark_mv
        self.dark_mv = None

        if self.dark_map is not None:
            del self.dark_map
        self.dark_map = None

        if self.flat_map is not None:
            del self.flat_map
        self.flat_map = None

        if self.linearity_table_left is not None:
            del self.linearity_table_left
        self.linearity_table_left = None
        if self.linearity_table_right is not None:
            del self.linearity_table_right
        self.linearity_table_right = None

        if self.readnoise_map is not None:
            del self.readnoise_map
        self.readnoise_map = None

    def add_bad_pixel_mask(self, detector, cdp_ftp_path=SIM_CDP_FTP_PATH,
                           cdp_version=''):
        """
        
        Add a bad pixel mask associated with the detector.
        
        :Parameters:
        
        detector: string
            The detector name with which to look up the bad pixel mask.
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific bad pixel mask CDP version number of the form 'x.y.z'.
        
        """        
        # If a null detector label is given, clear the bad pixel map.
        if detector is None:
            if self.bad_pixels is not None:
                del self.bad_pixels
            self.bad_pixels = None
            self.pixels.set_pedestal( None )
            self.bad_pixel_filename = ''
            return            

        # Get the bad pixel mask, which must be the same size as the
        # illuminated part of the detector. Since there is one bad pixel
        # mask per FPM but test data can come in various sizes extract the
        # relevant window from the mask. The bad pixel mask includes the
        # normal detector pixels (including the reference columns) but does
        # not include the reference rows added to the level 1 FITS data.
        window = (1, 1, self.illuminated_shape[0], self.detector_shape[1])

        # Get the required CDP version numbers
        (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )

        # Read a MIRI bad pixel mask in standard STScI CDP format
        # and ensure the mask refers to the correct detector.
        bad_pixel_mask = get_cdp('MASK', detector=detector,
                                 ftp_path=cdp_ftp_path,
                                 ftp_user=SIM_CDP_FTP_USER, 
                                 ftp_passwd=SIM_CDP_FTP_PASSWD,
                                 cdprelease=cdprelease,
                                 cdpversion=cdpversion,
                                 cdpsubversion=cdpsubversion,
                                 logger=self.toplogger)
        if bad_pixel_mask is None:
            strg = "Could not find a bad pixel mask CDP for detector %s." % detector
            strg += " No bad pixels will be simulated."
            self.logger.error(strg)
            self.bad_pixels = None
            self.pixels.set_pedestal( None )
            self.bad_pixel_filename = ''
            return

        self.bad_pixel_filename = bad_pixel_mask.meta.filename
        bp_detectorid = bad_pixel_mask.meta.instrument.detector
        if bp_detectorid is not None:
            if bp_detectorid != self.detectorid:
                strg = "Bad pixel mask is for the wrong detector "
                strg += "(%s rather than %s)." % (bp_detectorid,
                                                  self.detectorid)
                self.logger.warn( "***%s" % strg )
                #raise ValueError(strg)

        # Define reference columns (if any) as NON-SCIENCE.
        if self.left_columns > 0:
            bad_pixel_mask.dq[:,:self.left_columns] = MASK_NON_SCIENCE
        if self.right_columns > 0:
            bad_pixel_mask.dq[:,-self.right_columns:] = MASK_NON_SCIENCE

        # Define the bad pixel array.
        if self.bad_pixels is not None:
            del self.bad_pixels
        subarray_mask = bad_pixel_mask.get_subarray(window)
        
        # If there are any reference rows, extend the bad pixel array to
        # include the additional pixels in the reference rows, otherwise
        # use the bad pixel map as is.
        if self.top_rows > 0:
            self.bad_pixels = MASK_NON_SCIENCE * np.ones([self.detector_shape[0],
                                self.detector_shape[1]],
                                dtype=np.uint32)
            self.bad_pixels[self.bottom_rows:-self.top_rows,:] = subarray_mask
        elif self.bottom_rows > 0:
            self.bad_pixels = MASK_NON_SCIENCE * np.ones([self.detector_shape[0],
                                self.detector_shape[1]],
                                dtype=np.uint32)
            self.bad_pixels[self.bottom_rows:,:] = subarray_mask
        else:
            self.bad_pixels = subarray_mask
        
        # Plot the bad pixel map if requested.
        if self._verbose > 1 and self._makeplot:
            tstrg = "Bad pixel map obtained from " + \
                os.path.basename(self.bad_pixel_filename)
            mplt.plot_image(self.bad_pixels, title=tstrg)

    def add_gain_map(self, detector, cdp_ftp_path=SIM_CDP_FTP_PATH,
                     cdp_version=''):
        """
        
        Add a gain map associated with the detector.

        :Parameters:
        
        detector: string
            The detector name with which to look up the gain map.
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific gain CDP version number of the form 'x.y.z'.
        
        """
        # If a null detector label is given, clear the gain map.
        if detector is None:
            self.gain_map = None
            self.gain_map_filename = ''
            self.mean_gain = 1.0
            return

        # Get the required CDP version numbers
        (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )
         
        # Get the gain map, which must be the same size as the
        # illuminated part of the detector. Since there is one gain map
        # mask per FPM but test data can come in various sizes extract the
        # relevant window from the map. The gain map includes the
        # normal detector pixels (including the reference columns) but does
        # not include the reference rows added to the level 1 FITS data.        
        gain_model = get_cdp('GAIN', detector=detector,
                             ftp_path=cdp_ftp_path,
                             ftp_user=SIM_CDP_FTP_USER, 
                             ftp_passwd=SIM_CDP_FTP_PASSWD,
                             cdprelease=cdprelease,
                             cdpversion=cdpversion,
                             cdpsubversion=cdpsubversion,
                             logger=self.toplogger)
        if gain_model is None:
            strg = "Could not find gain CDP for detector %s." % detector
            strg += " A gain of 1.0 (e/DN) will be assumed."
            self.logger.error(strg)
            self.gain_map = None
            self.gain_map_filename = ''
            self.mean_gain = 1.0
            return

        self.gain_map_filename = gain_model.meta.filename
        # The CDP data includes the reference columns but not the
        # reference rows.
        gain_map = gain_model.data[0:self.illuminated_shape[0],
                                   0:self.detector_shape[1]]
        mean_gain = np.mean(gain_map)

        # If there are any reference rows, extend the gain map to
        # include the additional pixels in the reference rows, otherwise
        # use the gain map as is.
        if self.top_rows > 0:
            new_map = mean_gain * np.ones([self.detector_shape[0],
                                           self.detector_shape[1]])
            new_map[self.bottom_rows:-self.top_rows,:] = gain_map
            gain_map = new_map
        elif self.bottom_rows > 0:
            new_map = mean_gain * np.ones([self.detector_shape[0],
                                           self.detector_shape[1]])
            new_map[self.bottom_rows:,:] = gain_map
            gain_map = new_map

        # Plot the gain map if requested.
        if self._verbose > 1 and self._makeplot:
            tstrg = "Gain map obtained from " + \
                os.path.basename(self.gain_map_filename)
            mplt.plot_image(gain_map, title=tstrg)

        if self.gain_map is not None:
            del self.gain_map
        self.gain_map = gain_map
        self.mean_gain = mean_gain
        del gain_model

# # Original version which used a dark supplied with the SCASim release.
#     def add_dark_map_fixed(self, filename, readpatt=None):
#         """
#         
#         Add a fixed dark map associated with the detector.
#         The dark current map is scaled to the expected dark current
#         in electrons/s.
#        
#         This version loads a fixed averaged dark which is supplied
#         with the SCASim release.
# 
#         :Parameters:
#         
#         filename: string
#             The name of the file from which to load the dark map.
#         readpatt: string
#             Detector readout mode to which the dark map applies.
#        
#         """
#         # If a null filename is given, clear the dark map.
#         if filename is None or not filename:
#             self.dark_map = None
#             self.dark_map_filename = ''
#             self.dark_averaged = True
#             self.mean_dark = 0.0
#             return
#         
#         # The DARK calibration used here is generated from the CDP DARK
#         # and supplied with SCASim.
#         if self._verbose > 1:
#             strg = "Reading \'DARK\' model from \'%s\'" % filename
#             self.logger.info( strg )
# 
#         # Get the dark current map, which must be the same size as the
#         # illuminated part of the detector. Since there is one dark map
#         # mask per FPM but test data can come in various sizes extract the
#         # relevant window from the mask. The dark map includes the
#         # normal detector pixels (including the reference columns) but does
#         # not include the reference rows added to the level 1 FITS data.
#         #dark_model = MiriMeasuredModel( filename )
#         dark_model = MiriDarkReferenceModel( filename )
#         nints = dark_model.data.shape[0]
#         # Multiply by the gain to convert the read noise from DN into electrons.
#         # The CDP data includes the reference columns but not the
#         # reference rows.
#         dark_map = dark_model.data[:, 0:self.illuminated_shape[0],
#                                    0:self.detector_shape[1]] * self.mean_gain
#         self.dark_map_filename = filename
# 
#         # Scale the averaged dark map to the expected level in electrons/s.
#         expected_level = self._sca['DARK_CURRENT']
#         pos_values = np.where(dark_map > 0.0)
#         actual_level = np.median(dark_map[pos_values])
# # CANNOT TEST THE FOLLOWING - MemoryError
# #         actual_std = np.std(dark_map[pos_values])
# #         good_values = np.where(pos_values < (actual_mean + 3.0*actual_std))
# #         good_mean = np.mean(dark_map[good_values])
#         if self._verbose > 1:
#             self.logger.debug("Scaling the DARK by %.4g" % \
#                              (expected_level / actual_level))
#         dark_map = dark_map * expected_level / actual_level
# 
#         # Remove negative values from the dark map.
#         neg_values = np.where(dark_map < 0.0)
#         dark_map[neg_values] = 0.0
#         self.mean_dark = np.mean(dark_map)
# 
#         # If there are any reference rows, extend the dark map to
#         # include the additional pixels in the reference rows, otherwise
#         # use the dark map as is.
#         if self.top_rows > 0:
#             new_map = self.mean_dark * np.ones([nints, self.detector_shape[0],
#                                                 self.detector_shape[1]])
#             new_map[:,self.bottom_rows:-self.top_rows,:] = dark_map
#             dark_map = new_map
#         elif self.bottom_rows > 0:
#             new_map = self.mean_dark * np.ones([nints, self.detector_shape[0],
#                                                 self.detector_shape[1]])
#             new_map[:,self.bottom_rows:,:] = dark_map
#             dark_map = new_map
# 
#         # Plot the dark map if requested.
#         if self._verbose > 1 and self._makeplot:
#             tstrg = "Dark map obtained from " + \
#                 os.path.basename(filename)
#             mplt.plot_image(dark_map, title=tstrg)
#             #mplt.plot_hist(dark_map.ravel(), title=tstrg)
# 
#         if self.dark_map is not None:
#             del self.dark_map
#         self.dark_map = dark_map
#         self.dark_averaged = True
#         del dark_model

# This version fetches a DARK CDP from the repository.
    def add_dark_map_cdp(self, detector, readpatt=None, subarray=None,
                         averaged=False, cdp_ftp_path=SIM_CDP_FTP_PATH,
                         cdp_version=''):
        """
         
        Add a dark current map associated with the detector.
         
        This version loads a suitable dark from the MIRI CDP repository.
         
        NOTE: add_gain_map should be called first to obtain the mean gain,
        otherwise a value of 1.0 will be assumed.
 
        :Parameters:
         
        detector: string
            The detector name with which to look up the dark CDP.
        readpatt: string
            Detector readout mode to which the dark CDP applies.
        subarray: string, optional
            The subarray with which to look up the dark CDP
        averaged: bool, optional
            If True, an averaged 3-D version of the DARK model is obtained,
            otherwise the full 4-D model is obtained. 
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific dark map CDP version number of the form 'x.y.z'.
        
        """
        # If a null detector label is given, clear the dark map.
        if detector is None:
            self.dark_map = None
            self.dark_map_filename = ''
            self.dark_averaged = averaged
            self.mean_dark = 0.0
            return
 
#         # Get the required CDP version numbers
#         (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )
         
        # Some of the DARK files do not have standard CDP file names, so we
        # need to access them directly using the CDP interface class
        CDPInterface = MiriCDPInterface(ftp_path=cdp_ftp_path,
                                        ftp_user=SIM_CDP_FTP_USER, 
                                        ftp_passwd=SIM_CDP_FTP_PASSWD,
                                        logger=self.toplogger)
         
        # Refresh the interface if any parameters are different from when the
        # class was first created (necessary when the class is a singleton).
        CDPInterface.refresh(ftp_path=cdp_ftp_path,
                             ftp_user=SIM_CDP_FTP_USER, 
                             ftp_passwd=SIM_CDP_FTP_PASSWD)
         
        # The DARK maps must match the detector and readout pattern
        # and be the averaged versions of the CDP originals.
        must_contain = []
        must_not_contain = []
        if detector:
            must_contain.append(detector)
        if readpatt:
            must_contain.append(readpatt)
        if subarray and subarray != 'FULL':
            must_contain.append(subarray)
        else:
            for suba in MIRI_SUBARRAYS:
                must_not_contain.append(suba)                      
        must_contain.append('DARK')
        if averaged:
            must_contain.append('averaged')
        else:
            must_not_contain.append('averaged')
        if cdp_version:
            must_contain.append('cdp_version')
        (matched_files, matched_folders) = \
            CDPInterface.match_cdp_substrings(mustcontain=must_contain,
                                              mustnotcontain=must_not_contain)
         
        if len(matched_files) > 1:
            # More than one match. Take the last file.
            filename = matched_files[-1]
            ftp_path = matched_folders[-1]
        elif len(matched_files) > 0:
            # Exactly one match. Take the first (and only) file.
            filename = matched_files[0]
            ftp_path = matched_folders[0]
        else:
            # No CDPs were matched
            filename = None
             
        if filename:
            # Update the local CDP cache to make sure it contains the specified file,
            # and obtain the local file path and name.
            local_filename = CDPInterface.update_cache(filename, ftp_path)
            if averaged:
                strg = "Reading averaged DARK model from \'%s\'" % local_filename
            else:
                strg = "Reading DARK model from \'%s\'" % local_filename
            self.logger.info(strg)
            dark_model = MiriDarkReferenceModel( init=local_filename )
        else:
            dark_model = None
 
        del CDPInterface
        if dark_model is None:
            if averaged:
                strg = "Could not find averaged DARK CDP for detector %s" % detector
            else:
                strg = "Could not find DARK CDP for detector %s" % detector
            if readpatt:
                strg += " with %s readout mode" % readpatt
            if subarray:
                strg += " for %s subarray" % subarray
            strg += ". No dark current will be simulated."
            self.logger.error(strg)
            self.dark_map = None
            self.dark_map_filename = ''
            self.mean_dark = 0.0
            self.dark_averaged = averaged
            return
        
        self.dark_map_filename = local_filename
        nints = dark_model.data.shape[0]
        # The dark map is managed differently when it is averaged.
        if averaged:
            #ngroups = dark_model.data.shape[0]
            # Multiply by the gain to convert the read noise from DN into electrons.
            # The CDP data includes the reference columns but not the
            # reference rows.
            dark_map = dark_model.data[:, 0:self.illuminated_shape[0],
                                       0:self.detector_shape[1]] * self.mean_gain
 
            # Scale an averaged dark map to the expected level in electrons/s.
            expected_level = self._sca['DARK_CURRENT']
            pos_values = np.where(dark_map > 0.0)
            actual_level = np.median(dark_map[pos_values])
# CANNOT TEST THE FOLLOWING - MemoryError
#             actual_std = np.std(dark_map[pos_values])
#             good_values = np.where(pos_values < (actual_mean + 3.0*actual_std))
#             good_mean = np.mean(dark_map[good_values])
            if self._verbose > 1:
                self.logger.debug("Scaling the DARK by %.4g" % \
                                 (expected_level / actual_level))
            dark_map = dark_map * expected_level / actual_level
 
            # Remove negative values from the dark map.
            neg_values = np.where(dark_map < 0.0)
            dark_map[neg_values] = 0.0
            self.mean_dark = np.mean(dark_map)
     
            # If there are any reference rows, extend the dark map to
            # include the additional pixels in the reference rows, otherwise
            # use the dark map as is.
            if self.top_rows > 0:
                new_map = self.mean_dark * np.ones([nints, self.detector_shape[0],
                                               self.detector_shape[1]])
                new_map[:,self.bottom_rows:-self.top_rows,:] = dark_map
                dark_map = new_map
            elif self.bottom_rows > 0:
                new_map = self.mean_dark * np.ones([nints, self.detector_shape[0],
                                                self.detector_shape[1]])
                new_map[:,self.bottom_rows:,:] = dark_map
                dark_map = new_map
            self.dark_averaged = True
 
            # Plot the dark map if requested.
            if self._verbose > 1 and self._makeplot:
                tstrg = "Dark map obtained from " + \
                    os.path.basename(self.dark_map_filename)
                mplt.plot_image(dark_map, title=tstrg)
                #mplt.plot_hist(dark_map.ravel(), title=tstrg)
        else:
            # The full-sized DARK model is added to the final ramp data in DN,
            # so it does not need to be scaled by the gain.
            # Since the full-sized DARK is already very large, it is not
            # extended to include reference rows.
            self.dark_averaged = False
            dark_map = dark_model.data
            self.mean_dark = np.mean(dark_map[0,0,:,:])

            # Plot the dark map if requested.
            # Only a subset of the full data can be plotted.
            if self._verbose > 1 and self._makeplot:
                tstrg = "First group of dark map obtained from " + \
                    os.path.basename(self.dark_map_filename)
                mplt.plot_image(dark_map[0,0,:,:], title=tstrg)
 
        if self.dark_map is not None:
            del self.dark_map
        self.dark_map = dark_map
        del dark_model

    def add_flat_map(self, detector, readpatt=None, subarray=None,
                     mirifilter=None, miriband=None,
                     cdp_ftp_path=SIM_CDP_FTP_PATH, cdp_version=''):
        """
        
        Add a pixel flat-field associated with the detector.

        :Parameters:
        
        detector: string
            The detector name with which to look up the pixel flat-field CDP.
        readpatt: string, optional
            The readout mode with which to look up the pixel flat-field CDP.
        subarray: string, optional
            The subarray with which to look up the pixel flat-field CDP.
        mirifilter: string, optional
            The filter name with which to look up a pixel flat-field CDP.
        miriband: string, optional
            The band name with which to look up a flat-field
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific pixel flat-field CDP version number of the form 'x.y.z'.
       
        """
        # If a null detector label is given, clear the flat field.
        if detector is None:
            self.flat_map = None
            self.flat_map_filename = ''
            return

        strg = "Find a PIXELFLAT for detector=%s, readpatt=%s, subarray=%s" % \
            (detector, str(readpatt), str(subarray))
        strg += " filter=%s, band=%s" %  (str(mirifilter), str(miriband))
        self.logger.debug(strg)

        # Get the pixel flat-field, which must either (in FULL-frame mode)
        # be the same size as the illuminated part of the detector or
        # (in subarray mode) be the same size as the subarray.
        # Since test data can come in various sizes the flat-field is
        # enlarged to full size and the relevant window extracted from
        # the data. The flat-field includes the normal detector pixels
        # (including the reference columns) but does not include the
        # reference output rows.
        if subarray is not None and subarray != 'FULL':
            (firstrow, firstcol, expected_rows, expected_cols) = \
                detector_properties.get('SUBARRAY', subarray)
        else:
            (expected_rows, expected_cols) = \
                (self.illuminated_shape[0], self.detector_shape[1])

        # Get the required CDP version numbers
        (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )

        # Flat-fields are (for historical reasons) available for only a subset
        # of detector readout modes, and have a variety of different file
        # names. The following search paths are executed to find a suitable CDP.
        type_search = ['PIXELFLAT', 'FLAT']
        readpatt_search = [readpatt, 'ANY']
        if mirifilter:
            if mirifilter == 'P750L':
                # If the 'P750L' filter has been explicitly specified, only
                # a CDP made specifically for that filter will be valid.
                filter_search = [mirifilter]
                #readpatt_search = ['ANY']
            else:
                filter_search = [mirifilter, 'GENERIC', 'ANY']
        else:
            filter_search = ['GENERIC', 'ANY']
        if miriband:
            band_search = [miriband, 'ANY']
        else:
            band_search = ['ANY']
        nattempts = 1
        logstrg = ""
        reported = False
        for try_type in type_search:
            for try_readpatt in readpatt_search:
                for try_filter in filter_search:
                    for try_band in band_search:
                        flat_model = get_cdp(try_type, detector=detector,
                                        readpatt=try_readpatt,
                                        subarray=subarray,
                                        mirifilter=try_filter,
                                        band=try_band,
                                        ftp_path=cdp_ftp_path,
                                        ftp_user=SIM_CDP_FTP_USER, 
                                        ftp_passwd=SIM_CDP_FTP_PASSWD,
                                        cdprelease=cdprelease,
                                        cdpversion=cdpversion,
                                        cdpsubversion=cdpsubversion,
                                        logger=self.toplogger,
                                        fail_message=False)
                        if flat_model is not None:
                            # Check the shape
                            if flat_model.data.shape == \
                               (expected_rows, expected_cols):
                                # Success - break out of the loop
                                break
                            else:
                                # Wrong size. Delete this model and keep looking.
                                if not reported:
                                    # Report a size problem once.
                                    logstrg = "Pixel flat has the wrong size: "
                                    logstrg += "%s instead of %s.\n" % \
                                        (str( flat_model.data.shape ),
                                         str( (expected_rows, expected_cols) ))
                                    reported = True
                                del flat_model
                                flat_model = None
                        nattempts += 1
                    if flat_model is not None:
                        break
                if flat_model is not None:
                    break
            if flat_model is not None:
                break

        if flat_model is None:
            logstrg += "Could not find suitable pixel flat-field"
            logstrg += " for detector %s" % detector
            if subarray:
                logstrg += " and subarray %s" % subarray
            logstrg += ". No flat-field will be applied."
            self.logger.warn(logstrg)
            self.flat_map = None
            self.flat_map_filename = ''
            return
        elif nattempts > 1:
            logstrg += "Could not find exact match for pixel flat-field"
            logstrg += " for detector %s" % detector
            alternative = False
            if readpatt:
                logstrg += " with %s mode" % readpatt
                alternative = True
            if mirifilter:
                logstrg += " with filter=\'%s\'" % mirifilter
                alternative = True
            if miriband:
                logstrg += " with band=\'%s\'" % miriband
                alternative = True
            if alternative:
                logstrg += ". An alternative is being used."
            self.logger.warn(logstrg)
                
        self.flat_map_filename = flat_model.meta.filename
        # Fill the masked parts of the flat-field. Do not reshape the
        # array yet, because it might be subarray-sized.
        flat_map = flat_model.data_filled
        
        # Remove any remaining NaNs which not been removed by filling
        # the masked parts of the flat-field.
        wherenan = np.where( np.isnan(flat_map) )
        flat_map[wherenan] = 1.0
     
        # Set the non-science parts of the flat-field to zero.
        # This will add a mask showing the non-illuminated areas.
        NON_SCIENCE = 2
        wherenonsci = np.where( (flat_model.dq & NON_SCIENCE) > 0 )
        flat_map[wherenonsci] = 0.0

        # If the flat-field is for a subarray mode, extend it to cover
        # the full detector.
        if subarray is not None and subarray != 'FULL':
            (firstrow, firstcol, subrows, subcols) = \
                detector_properties.get('SUBARRAY', subarray)
            new_map = np.ones([self.illuminated_shape[0],
                                self.detector_shape[1]])
            # Subarrays locations start from 1 instead of 0
            rowstart = int(firstrow) - 1
            colstart = int(firstcol) - 1
            rowend = rowstart + subrows
            colend = colstart + subcols
            new_map[rowstart:rowend, colstart:colend] = flat_map
            flat_map = new_map
        else:
            # A full-frame CDP data includes the reference columns but
            # not the reference rows.
            new_map = flat_map[0:self.illuminated_shape[0],
                               0:self.detector_shape[1]]
            flat_map = new_map
            

        # If there are any reference rows, extend the flat-field to
        # include the additional pixels in the reference rows, otherwise
        # use the flat-fild as is.
        if self.top_rows > 0:
            new_map = np.ones([self.detector_shape[0],
                               self.detector_shape[1]])
            new_map[self.bottom_rows:-self.top_rows,:] = flat_map
            flat_map = new_map
        elif self.bottom_rows > 0:
            new_map = np.ones([self.detector_shape[0],
                               self.detector_shape[1]])
            new_map[self.bottom_rows:,:] = flat_map
            flat_map = new_map

        # Plot the flat-field if requested.
        if self._verbose > 1 and self._makeplot:
            tstrg = "Flat-field obtained from " + \
                os.path.basename(self.flat_map_filename)
            mplt.plot_image(flat_map, title=tstrg)

        if self.flat_map is not None:
            del self.flat_map
        self.flat_map = flat_map
        del flat_model

    def add_linearity_table(self, detector, mirifilter=None, miriband=None,
                     cdp_ftp_path=SIM_CDP_FTP_PATH, cdp_version=''):
        """
        
        Add a linearity correction associated with the detector.

        :Parameters:
        
        detector: string
            The detector name with which to look up the linearity
            correction CDP.
        mirifilter: string, optional
            The filter name with which to look up the linearity
            correction CDP
        miriband: string, optional
            The band name with which to look up the linearity
            correction CDP.
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific pixel flat-field CDP version number of the form 'x.y.z'.
       
        """
        # If a null detector label is given, clear the linearity correction.
        if detector is None:
            self.linearity_table_left = None
            self.linearity_table_right = None
            self.linearity_filename = ''
            return

        strg = "Find a LINEARITY CDP for detector=%s" % detector
        strg += " filter=%s, band=%s" %  (str(mirifilter), str(miriband))
        self.logger.debug(strg)
        
        # Get the required CDP version numbers
        (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )

        # Get the linearity correction.        
        linearity_model = get_cdp('LINEARITY', detector=detector,
                             mirifilter=mirifilter, band=miriband,
                             ftp_path=cdp_ftp_path,
                             ftp_user=SIM_CDP_FTP_USER, 
                             ftp_passwd=SIM_CDP_FTP_PASSWD,
                             cdprelease=cdprelease,
                             cdpversion=cdpversion,
                             cdpsubversion=cdpsubversion,
                             logger=self.toplogger)
        if linearity_model is None and (mirifilter or miriband):
            # If a particular filter or band could not be matched,
            # try and find a CDP matching any filter or band.
            linearity_model = get_cdp('LINEARITY', detector=detector,
                                 mirifilter=None, band=None,
                                 ftp_path=cdp_ftp_path,
                                 ftp_user=SIM_CDP_FTP_USER, 
                                 ftp_passwd=SIM_CDP_FTP_PASSWD,
                                 cdprelease=cdprelease,
                                 cdpversion=cdpversion,
                                 cdpsubversion=cdpsubversion,
                                 logger=self.toplogger)
        if linearity_model is None:
            strg = "Could not find linearity CDP for detector %s" % detector
            if mirifilter:
                strg += " with filter=\'%s\'" % mirifilter
                alternative = True
            if miriband:
                strg += " with band=\'%s\'" % miriband
                alternative = True
            strg += ". There will be no linearity correction."
            self.logger.error(strg)
            self.linearity_table_left = None
            self.linearity_table_right = None
            self.linearity_filename = ''
            return

        self.linearity_filename = linearity_model.meta.filename
        # Extract the linearity tables from the left and right halves of 
        # the linearity CDP.
        ncolumns = linearity_model.data.shape[-1]
        nrows = linearity_model.data.shape[-2]
        leftcol = ncolumns//4
        rightcol = leftcol + ncolumns//2
        row = nrows//2
        self.linearity_table_left = \
            linearity_model.get_reverse_table(row, leftcol)
        self.linearity_table_right = \
            linearity_model.get_reverse_table(row, rightcol)

        # Plot the linearity tables if requested.
        if self._verbose > 1 and self._makeplot:
            tstrg = "linearity table obtained from " + \
                os.path.basename(self.linearity_filename)
            mplt.plot_xy( None, self.linearity_table_left,
                          xlabel='Linear DN', ylabel='Nonlinear DN',
                          title="LEFT " + tstrg)
            mplt.plot_xy( None, self.linearity_table_right,
                          xlabel='Linear DN', ylabel='Nonlinear DN',
                          title="RIGHT " + tstrg)
# Does not work because left and right tables have difference lengths
#             mplt.plot_xycolumn( None, [self.linearity_table_left,
#                                 self.linearity_table_right], title=tstrg)

        del linearity_model

    def add_readnoise_map(self, detector, readpatt=None,
                          cdp_ftp_path=SIM_CDP_FTP_PATH, cdp_version=''):
        """
        
        Add a read noise map associated with the detector.
        
        NOTE: add_gain_map should be called first to obtain the mean gain,
        otherwise a value of 1.0 will be assumed.

        :Parameters:
        
        detector: string
            The detector name with which to look up the readnoise map.
        readpatt: string, optional
            The readout mode with which to look up the readnoise map.
        cdp_ftp_path: str, optional, default=SIM_CDP_FTP_PATH
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.        
        cdp_version: string, optional, default=''
            A specific readnoise CDP version number of the form 'x.y.z'.
        
        """
        # If a null detector label is given, clear the read noise map.
        if detector is None:
            self.readnoise_map = None
            self.readnoise_map_filename = ''
            return

        # Get the required CDP version numbers
        (cdprelease, cdpversion, cdpsubversion) = cdp_version_decode( cdp_version )
        
        # Get the read noise map, which must be the same size as the
        # illuminated part of the detector. Since there is one read noise map
        # mask per FPM but test data can come in various sizes extract the
        # relevant window from the map. The read noise map includes the
        # normal detector pixels (including the reference columns) but does
        # not include the reference rows added to the level 1 FITS data.
        readnoise_model = get_cdp('READNOISE', detector=detector,
                                  readpatt=readpatt,
                                  ftp_path=cdp_ftp_path,
                                  ftp_user=SIM_CDP_FTP_USER, 
                                  ftp_passwd=SIM_CDP_FTP_PASSWD,
                                  cdprelease=cdprelease,
                                  cdpversion=cdpversion,
                                  cdpsubversion=cdpsubversion,
                                  logger=self.toplogger)
        if readnoise_model is None:
            strg = "Could not find read noise CDP for detector %s" % detector
            if readpatt:
                strg += " with %s readout mode" % readpatt
            strg += ". No read noise will be simulated."
            self.logger.error(strg)
            self.readnoise_map = None
            self.readnoise_map_filename = ''
            return

        self.readnoise_map_filename = readnoise_model.meta.filename
        # Multiply by the gain to convert the read noise from DN into electrons.
        readnoise_map = readnoise_model.data[0:self.illuminated_shape[0],
                                             0:self.detector_shape[1]] * \
                                            self.mean_gain

        rnshape = readnoise_map.shape
        if rnshape[0] == self.illuminated_shape[0] and \
           rnshape[1] == self.detector_shape[1]:
            # If there are any reference rows, extend the readnoise map to
            # include the additional pixels in the reference rows, otherwise
            # use the dark map as is.
            if self.top_rows > 0:
                new_map = np.zeros([self.detector_shape[0],
                                    self.detector_shape[1]])
                new_map[self.bottom_rows:-self.top_rows,:] = readnoise_map
                readnoise_map = new_map
            elif self.bottom_rows > 0:
                new_map = np.zeros([self.detector_shape[0],
                                         self.detector_shape[1]])
                new_map[self.bottom_rows:,:] = readnoise_map
                readnoise_map = new_map
        else:
            strg = "Read noise map has the wrong size (%d x %d). " % rnshape
            strg += "It should be %d x %d." % (self.illuminated_shape[0],
                                               self.detector_shape[1])
            raise TypeError(strg)

        # Plot the readnoise map if requested.
        if self._verbose > 1 and self._makeplot:
            tstrg = "Read noise map obtained from " + \
                os.path.basename(self.readnoise_map_filename)
            mplt.plot_image(readnoise_map, title=tstrg)

        if self.readnoise_map is not None:
            del self.readnoise_map
        self.readnoise_map = readnoise_map
        
    def get_subarray_shape(self, subarray=None):
        """
        
        Get the shape of the data generated by a subarray mode.
        
        The function starts with the basic subarray size and adds
        the necessary number of reference rows.
        
        :Parameter:
        
        subarray: tuple of 4 ints, optional, default is None
            If None a full frame readout is assumed. Otherwise this
            parameter should be set to subarray parameters
            (firstrow, firstcol, subrows, subcolumns).
            *NOTE: Rows and columns are numbered from 1 when describing
            a subarray.*
        
        :Returns:
        
        subarray_shape: tuple of 2 ints
            The shape of the subarray data in pixels (rows, columns)
            
        """
        if subarray is None:
            detector_columns = self.left_columns + self.illuminated_shape[1] + \
                               self.right_columns
            refout_rows = self.illuminated_shape[0]
            refout_columns = self.top_rows * detector_columns / refout_rows
            self.refout_shape = (refout_rows, refout_columns)
            return self.detector_shape
        else:
            subrows_illuminated = subarray[2]
            subrows_bottom = int(subrows_illuminated * \
                                 self.bottom_rows / self.illuminated_shape[0])
            subrows_top = int(subrows_illuminated * \
                              self.top_rows / self.illuminated_shape[0])
            subrows = subrows_bottom + subrows_illuminated + subrows_top
            subcols = subarray[3]
            
            refout_rows = subrows_illuminated
            refout_columns = int(subrows_top * subcols / refout_rows)
            self.refout_shape = (refout_rows, refout_columns)
           
            return (subrows, subcols)

    def get_counts(self):
        """
        
        Return the current electron count for the detector pixels.
        
        """
        return self.pixels.get_counts()

    def set_metadata(self, metadata):
        """
        
        Add detector keywords to the given MIRI metadata.
        
        :Parameters:
        
        metadata: dictionary-like object
            A keyword-addressable metadata object to which detector
            keywords should be added.
            
        :Returns:
        
        metadata: dictionary-like object
            An updated metadata object.
            
        """        
        # Find out if the metadata object is a plain dictionary or
        # an object which supports comments and history records.
        comments_possible = hasattr(metadata, 'add_comment')
#         history_possible = hasattr(metadata, 'add_history')
        
        #metadata["FPM_ID"] = self.fpm_id
        #metadata["FPM"] = self.fpm_name
        metadata["SCA_ID"] = str(self.sca_id)
        metadata["DETECTOR"] = str(self.detector_name)
        #metadata["CHIP"] = self.chip
        
        metadata["WELL"] = self.well_depth
        if self.simulate_dark_current:
            if self.dark_averaged:
                # Round the value to 4 decimal places
                combined_dark_current = self.dark_current * self.mean_dark
                metadata["DARKCURR"] = float("%.4f" % combined_dark_current)
        elif comments_possible:
            metadata.add_comment("Simulated dark current disabled.")

        if self.simulate_gain:
            # Report the average gain and read noise, rounded to 4 decimal places
            # metadata['GAINFN'] = 'LINEAR'
            if self.gain_map is None:
                metadata["GAINCF"] = 1.0
                metadata['MAXDN'] = int(min(self.well_depth, 65535))
            else:
                metadata["GAINCF"] = float("%.4f" % self.gain_map.mean())
                mingain = float(max(1.0, self.gain_map.min()))
                metadata['MAXDN'] = int(min(self.well_depth/mingain, 65535))
        elif comments_possible:
            metadata.add_comment("Simulated gain disabled.")

        if self.simulate_read_noise:
            if self.readnoise_map is None:
                metadata["RDNOISE"] = 0.0
            else:
                # Only determine the mean read noise in non-zero areas
                # (i.e. skipping non-illuminated  pixels)
                valid_rn = np.where(self.readnoise_map > 0.0)
                mean_rn = self.readnoise_map[valid_rn].mean()
                metadata["RDNOISE"] = float("%.4f" % mean_rn)
                if comments_possible and (self.noise_factor != 1.0):
                    metadata.add_comment( \
                        "Read noise scaled by %.4f to compensate for DARK noise." % \
                        self.noise_factor)       
        elif comments_possible:
            metadata.add_comment("Simulated read noise disabled.")

        if self.simulate_nonlinearity:
            if NONLINEARITY_BY_TABLE:
                metadata.add_comment("Non-linearity simulated by translation table.")
                comment = "Nonlinearity: " + \
                    os.path.basename(self.linearity_filename)
                metadata.add_comment(comment)
            else:
                metadata.add_comment("Non-linearity simulated by detector sensitivity.")          
            if EXTRA_METADATA:
                metadata["SENCOFF"] = str(self.pixels.sensitivity)
        elif comments_possible:
            metadata.add_comment("Simulated non-linearity effects disabled.")

        if self.simulate_drifts:
            if EXTRA_METADATA:
                metadata["ZPCOFFS"] = str(self._sca['ZP_SLOW'])
                metadata["ZPCOFFF"] = str(self._sca['ZP_FAST'])
        elif comments_possible:
            metadata.add_comment("Simulated detector drift effects disabled.")

        if self.simulate_latency:
            if EXTRA_METADATA:            
                metadata["PSCOFF"]  = str(self._sca['PERSISTENCE'])
                metadata["LTCOFFS"] = str(self._sca['LATENCY_SLOW'])
                metadata["LTCOFFF"] = str(self._sca['LATENCY_FAST'])
        elif comments_possible:
            metadata.add_comment("Simulated detector latency effects disabled.")
        
        if comments_possible:
            simulate_poisson_noise = self.pixels.simulate_poisson_noise
            if not simulate_poisson_noise:
                metadata.add_comment("Simulated Poisson noise is disabled.")

            if self.bad_pixels is None:
                metadata.add_comment("Bad pixel mask: None applied.")
            else:       
                comment = "Bad pixel mask: " + \
                    os.path.basename(self.bad_pixel_filename)
                metadata.add_comment(comment)
                all_locations = np.where(self.bad_pixels > 0)
                all_count = len(all_locations[0])
                bad_locations = np.where((self.bad_pixels & MASK_DO_NOT_USE) > 0)
                bad_count = len(bad_locations[0])
                strg = "Bad pixel mask defines %d pixels, of which %d are unusable." % \
                    (all_count, bad_count)
                metadata.add_comment(strg)

            if self.dark_map is None:
                strg = "Dark map: None. All pixels have equal dark current."
                metadata.add_comment(strg)
            else:
                comment = "Dark map: " + \
                    os.path.basename(self.dark_map_filename)
                metadata.add_comment(comment)
                strg = "Dark map ranges from %.2f to %.2f (e)." % \
                    (self.dark_map.min(), self.dark_map.max())
                metadata.add_comment(strg)
                strg = "Dark map mean=%.2f and std=%.2f (e)." % \
                    (np.mean(self.dark_map), np.std(self.dark_map))
                metadata.add_comment(strg)

            if self.flat_map is None:
                strg = "Pixel flat-field: None. All pixels have equal weight."
                metadata.add_comment(strg)
            else:
                comment = "Pixel flat: " + \
                    os.path.basename(self.flat_map_filename)
                metadata.add_comment(comment)
                strg = "Pixel flat-field ranges from %.2f to %.2f." % \
                    (self.flat_map.min(), self.flat_map.max())
                metadata.add_comment(strg)
                strg = "Pixel flat-field mean=%.2f and std=%.2f." % \
                    (np.mean(self.flat_map), np.std(self.flat_map))
                metadata.add_comment(strg)

            if self.gain_map is None:
                strg = "Gain map: None. 1.0 is assumed."
                metadata.add_comment(strg)
            else:
                comment = "Gain map: " + \
                    os.path.basename(self.gain_map_filename)
                metadata.add_comment(comment)
                strg = "Gain map ranges from %.2f to %.2f (e/DN)." % \
                    (self.gain_map.min(), self.gain_map.max())
                metadata.add_comment(strg)
                strg = "Gain map mean=%.2f and std=%.2f (e/DN)." % \
                    (np.mean(self.gain_map), np.std(self.gain_map))
                metadata.add_comment(strg)

        return metadata
     
    def reset(self, nresets=1, new_exposure=False):
        """
        
        Reset the detector
        
        """            
        # Reset the underlying Poisson integrator.
        if self._verbose > 3:
            self.logger.debug( "Resetting the detector." )
        self.pixels.reset(nresets=nresets, new_exposure=new_exposure)
        self.cosmic_ray_count = 0
        self.cosmic_ray_pixel_count = 0
            
    def set_seed(self, seedvalue=None):
        """
        
        Set the seed for the numpy random number generator.
        This function can be used while testing to ensure the set of
        random choices that follow are well defined.
        
        :Parameters:
        
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator.
            If not specified, a value of None will be sent, which
            randomises the seed.
            
        """
        np.random.seed(seedvalue)
 
    def set_readout_mode(self, samplesum, sampleskip=0, refpixsampleskip=3,
                         nframes=1):
        """
         
        Defines the SCA detector readout mode parameters which determine
        the read noise.
         
        :Parameters:
         
        samplesum: int
            The total number of samples when reading a pixel (which
            affects the read noise). Must be at least 1.
        sampleskip: int, optional, default=0
            The number of samples skipped before reading a pixel.
        refpixsampleskip: int, optional, default=3
            The number of samples skipped before reading a reference pixel.
        nframes: int, optional, default=1
            The number of frames per group. Normally 1 for MIRI data, but
            if greater than 1 this parameter further reduces the read
            noise for each group.
             
        """
        if int(samplesum) <= 0:
            strg = "Number of samples per pixel must be at least 1."
            raise ValueError(strg)
        if int(nframes) <= 0:
            strg = "Number of frames per group must be at least 1."
            raise ValueError(strg)
 
        self.samplesum = samplesum
        self.sampleskip = sampleskip
        self.refpixsampleskip = refpixsampleskip
        self.nframes = nframes

    def frame_time(self, samplesum, sampleskip, refpixsampleskip=3,
                   subarray=None, burst_mode=False):
        """
        
        Calculate the frame time for a given subarray and readout mode.
        
        :Parameters:
        
        samplesum: int
            Number of clock cycles to sample a pixel during readout.
            Derived from the readout mode.
        sampleskip: int
            Number of clock cycles to dwell on a pixel before reading it.
            Derived from the readout mode.
        refpixsampleskip: int, optional, default is 3
            Number of clock cycles to dwell on a reference pixel before reading it.
            Derived from the readout mode.
        subarray: tuple of 4 ints, optional, default is None
            If None a full frame readout is assumed. Otherwise this
            parameter should be set to subarray parameters
            (firstrow, firstcol, subrows, subcolumns).
            *NOTE: Rows and columns are numbered from 1 when describing
            a subarray.*
        burst_mode: boolean, optional, default is False
            True if a subarray is being read out in burst mode,
            which skips quickly over unwanted columns.
        
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
            
        :Returns:
        
        ftime: float
            Frame time in seconds.
            
        """
        # Default to the current readout mode
        if samplesum is None:
            samplesum = self.samplesum
        if sampleskip is None:
            sampleskip = self.sampleskip
            
        if int(samplesum) <= 0:
            strg = "Number of samples per readout (samplesum) must be at least 1."
            raise ValueError(strg)

        # First define the physical location of the subarray on the
        # array scanned by a single (normal) amplifier. (Since the
        # amplifiers readout in parallel, it is only necessary to
        # calculate the frame time for one amplifier.)
        # A full frame readout (subarray=None) is a special case.
        if subarray is None:
            rowstart = 1
            rowstop = self.illuminated_shape[0]
            colstart = 1
            colstop = (self.left_columns + self.illuminated_shape[1] + \
                       self.right_columns) / self._normal_amps
        else:
            if len(subarray) == 4:
                firstrow = subarray[0]
                firstcol = subarray[1]
                subrows = subarray[2]
                subcolumns = subarray[3]
        
                rowstart = firstrow
                rowstop = firstrow - 1 + subrows
                colstart = (firstcol-1)//self._normal_amps + 1
                colstop = (firstcol + subcolumns - 2)//self._normal_amps + 1
            else:
                strg = "subarray must be None (=full frame) or " \
                    "(firstrow,firstcol,subrows,subcolumns)"
                raise ValueError(strg)
            
        if self._verbose > 3:
            self.logger.debug( "Frame time calculated from" + \
                " colstart=" + str(colstart) + \
                " colstop=" + str(colstop) + \
                " rowstart=" + str(rowstart) + \
                " rowstop=" + str(rowstop) )

        # The frame RTI calculation is taken from Mike Ressler's
        # "MIRI FPS Exposure Time Calculations (SCE FPGA2)" document,
        # which updates the original description in the MIRI Operations
        # Concept Document.
        detectorrows = self.illuminated_shape[0]
        ampcolumns = 2 + self.illuminated_shape[1]/self._normal_amps
        resetwidth = self._sca['RESET_WIDTH']
        resetoverhead = self._sca['RESET_OVERHEAD']
        frame_clks = frame_rti(detectorrows, ampcolumns, colstart, colstop,
                               rowstart, rowstop, sampleskip, refpixsampleskip,
                               samplesum, resetwidth, resetoverhead, burst_mode) 

        # Multiply the number of clock cycles by the cycle time to get
        # the frame time.
        frame_time =  frame_clks * float(self._sca['CLOCK_TIME'])
        if self._verbose > 3:
            self.logger.debug( "frame_clks=" + str(frame_clks) + \
                               " frame_time=" + str(frame_time) )        
        return frame_time
    
    def clock_time(self):
        """
        
        Return the readout clock time in seconds.
        
        :Parameters:
        
        None
        
        :Returns:
        
        clock_time: float
            The detector clock time in seconds.
        
        """
        return self._sca['CLOCK_TIME']
        
    def exposure_time(self, nints, ngroups, samplesum, sampleskip,
                      refpixsampleskip=3, subarray=None, burst_mode=False,
                      add_initial_resets=False, frame_time=None, nframes=1,
                      groupgap=0 ):
        """
        
        Calculate the exposure time for a given number of
        integrations and groups, subarray and readout mode
        (with the ability to override the frame time).
        
        :Parameters:
        
        nints: int
            The number of integrations per exposure. Must be at least 1.
        ngroups: int
            The number of groups per integration. Must be at least 1.
        samplesum: int
            Number of clock cycles to sample a pixel during readout.
            Derived from the readout mode.
        sampleskip: int
            Number of clock cycles to dwell on a pixel before reading it.
            Derived from the readout mode.
        refpixsampleskip: int, optional, default is 3
            Number of clock cycles to dwell on a reference pixel before reading it.
            Derived from the readout mode.
        subarray: tuple of 4 ints, optional, default is None
            The subarray mode from which to determine the frame time.
            If the frame_time parameter is given explicitly (see below),
            the subarray parameter is ignored.
            If None a full frame readout is assumed. Otherwise this
            parameter should be set to subarray parameters
            (firstrow, firstcol, subrows, subcolumns).
            *NOTE: Rows and columns are numbered from 1.*
        burst_mode: boolean, optional, default is False
            True if a subarray is being read out in burst mode,
            which skips quickly over unwanted columns.
        add_initial_resets: boolean, optional, default is False
            True if this is the first integration and additional reset frames
            are to be added.
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter overrides the readout mode and
            subarray and defines the frame time explicitly. It can be used,
            for example, to calculate the exposure time when simulating
            full frame exposures over a subarray-sized subset.
        nframes: int, optional, default=1
            The number of frames per group. Normally 1 for MIRI data.
        groupgap: int, optional, default=0
            The number of frames dropped between groups. Normally 0 for
            MIRI data.
        
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
            
        :Returns:
        
        etime: tuple of 2 floats
            (exp_time, elapsed_time) where exp_time is the exposure time
            on signal and elapsed_time the total elasped time. The two
            are the same unless extra resets are used between integrations.
            
        """
        if int(nints) <= 0:
            strg = "Number of integrations must be at least 1."
            raise ValueError(strg)
        if int(ngroups) <= 0:
            strg = "Number of groups must be at least 1."
            raise ValueError(strg)
        if int(nframes) <= 0:
            strg = "Number of frames per group must be at least 1."
            raise ValueError(strg)
        if int(groupgap) < 0:
            strg = "Number of dropped frames per group must be positive."
            raise ValueError(strg)

        # Begin by calculating the frame time.
        if frame_time is None:
            ftime = self.frame_time(samplesum, sampleskip,
                                    refpixsampleskip=refpixsampleskip,
                                    subarray=subarray,
                                    burst_mode=burst_mode)
        else:
            ftime = frame_time
        
        # Calculate the group time and the dropped frames time.
        gtime = ftime * int(nframes)
        dtime = ftime * int(groupgap)
        
        # The integration time is the sum of the time taken for the
        # groups plus the time used by the dropped frames between the
        # groups. Multiply by nints to get the exposure time.
        # The MIRI detectors have the option of inserting FRAME_RESETS extra
        # reset frames in between integrations, which can make the elapsed
        # time slightly larger than exp_time. (Normally, FRAME_RESETS=0.)
        int_time = (int(ngroups) * gtime) + ((int(ngroups)-1) * dtime)
        exp_time = int(nints) * int_time
        elapsed_time = exp_time + \
            (int(self._sca['FRAME_RESETS']) * \
             ftime * (nints-1))
        if self._verbose > 3:
            strg = "int_time=" + str(int_time) + \
                " nints=" + str(nints) + \
                " exp_time=" + str(exp_time) + "s" \
                " frm_resets=" + str(self._sca['FRAME_RESETS']) + \
                " frm_time=" + str(ftime) + "s" \
                " elapsed_time=" + str(elapsed_time) + "s"
            self.logger.debug( strg )
            
        return (exp_time, elapsed_time)
        
    def hit_by_cosmic_rays(self, cosmic_ray_list, nframes=1):
        """
        
        Simulate a series of cosmic ray hits on the detector.
        
        :Parameters:
        
        cosmic_ray_list: list of CosmicRay objects
            A list of cosmic ray events hitting the detector.
        nframes: int, optional, default=1
            The number of frames per group. Normally 1 for MIRI data, but
            if greater than 1 this parameter will dilute the effect of a
            cosmic ray hit (since all the frames belonging to a group
            are averaged.
            
        """
        # The effect of each cosmic ray hit is diluted if more than 1
        # frames are averaged to make a group.
        if nframes <= 1:
            energy_mult = 1.0
        else:
            # The net effect is the sum of the probabilities of an event
            # happening within each frame multiplied by the effect of a
            # cosmic ray hit on that frame. (Note that the software doesn't
            # support more than one cosmic ray hit per group.)
            triang = np.array(list(range(1,nframes+1)))
            topsum = triang.sum()
            energy_mult = topsum / (nframes * nframes)

        if self._verbose > 2:
            nhits = len(cosmic_ray_list)
            if nhits > 0:
                strg = "%d cosmic ray hits to the detector" % nhits
                if nframes > 1:
                    strg += " (nframes=%d so energy multiplier=%f)" % \
                        (nframes, energy_mult)
                strg += "."
                self.logger.info( strg )
        
        for cosmic_ray in cosmic_ray_list:
            # The cosmic ray has hit a detector pixel.
            (row, column) = cosmic_ray.get_target_coords()
            throw = np.random.uniform(0.0, 1.0)
            if throw > detector_properties['COSMIC_RAY_LEAKAGE_FRACTION']:
                # A normal cosmic ray hit
                energy_map = energy_mult * cosmic_ray.get_hit_map()
                self.pixels.hit_by_cosmic_ray(energy_map, row, column)
                
                self.cosmic_ray_count += 1
                affected = np.where( energy_map >= 1 )
                if affected:
                    self.cosmic_ray_pixel_count += len(affected[0])
            else:
                # A rare negative cosmic ray event.
                energy = energy_mult * cosmic_ray.get_electrons()
                self.pixels.leak(energy, row, column)

    def integrate(self, photon_flux, time, intnum=0):
        """
        
        Integrate the detector with a certain amount of flux for a certain
        length of time
        
        :Parameters:
        
        photon_flux: array_like
            The photon flux on which to integrate in photons per time unit.
            This must be the same shape and size as the illuminated portion
            of the detector.
        time: float
            The integration time in time units.
        intnum: int
            Integration number, which determines which dark calibration is
            used. Must be 0 or greater.
            
        """
        # The photon flux must be the same size as the illuminated portion
        # of the detector. Other parts get zero illumination.
        if photon_flux.shape == self.illuminated_shape:

            # Map the photon flux onto the central illuminated columns of the
            # detector pixels. The left and right columns receive zero
            # illumination.
            detector_flux = np.zeros(self.detector_shape, dtype=photon_flux.dtype)
            if self.top_rows > 0 and self.right_columns > 0:
                detector_flux[self.bottom_rows:-self.top_rows,
                              self.left_columns:-self.right_columns] \
                   = photon_flux
            else:
                detector_flux[self.bottom_rows:, self.left_columns:] \
                   = photon_flux

            # The detector flux if effectively zero where pixels are bad/dead.
            if self.bad_pixels is not None:
                # Dead zones are assumed not to respond at all.
                dead_zones = np.where((self.bad_pixels & MASK_DEAD) > 0)
                detector_flux[dead_zones] = 0.0
                
            # Multiply the flux by the flat-field
            if self.flat_map is not None:
                detector_flux = detector_flux * self.flat_map

            # The dark current is applied here ONLY if it has been averaged.
            if self.dark_averaged:
                # The dark current affects all parts of the detector
                if self.dark_map is not None:
                    if intnum > 0:
                        dark_flux = self.dark_current * self.dark_map[1,:,:]
                    else:
                        dark_flux = self.dark_current * self.dark_map[0,:,:]
                else:
                    # No dark map, but there is a nominal dark current.
                    dark_flux = self.dark_current * self._sca['DARK_CURRENT']
                electron_flux = detector_flux + dark_flux
            else:
                # DARK simulation is deferred until later.
                electron_flux = detector_flux
                        
            if self._verbose > 3:
                strg = "Integrating on flux array with "
                if self.simulate_dark_current and self.dark_averaged:
                    strg += "darkmin=%.2f e darkmax=%.2f e "   % \
                        (dark_flux.min(), dark_flux.max())
                    strg += "detmin=%.2f e detmax=%.2f e " % \
                        (detector_flux.min(), detector_flux.max())
                strg += "min=%.2f e max=%.2f e for %.2f %s." % \
                    (electron_flux.min(), electron_flux.max(),
                     time, self.time_unit)
                self.logger.debug( strg )
            # >>> Integrate on the flux
            self.pixels.integrate(electron_flux, time)
        else:
            # Faulty flux array given - raise an exception.
            strg = "Photon flux array has the wrong shape: " \
                "(%d, %d) instead of (%d, %d)." % (photon_flux.shape[0],
                                                   photon_flux.shape[1],
                                                   self.illuminated_shape[0],
                                                   self.illuminated_shape[1])
            raise TypeError(strg)

    def wait(self, time, bgflux=0.0):
        """
        
        Rest the detector for a certain length of time.
        
        :Parameters:
        
        time: float
            The elapsed time in time units.
        bgflux: float (optional)
            A uniform background photon flux falling on the detector
            while it is idle. This flux can affect detector persistence
            by filling charge traps. By default there is no background
            flux.
            
        """
        # Wait for the elapsed time
        self.pixels.wait(time, bgflux=bgflux)
        
    def readout(self, subarray=None, total_samples=None, removeneg=True):
        """
        
        Read out the detector (non destructive)
        
        :Parameters:
        
        subarray: tuple of 4 ints, optional, default is None
            If None a full frame readout is assumed. Otherwise this
            parameter should be set to subarray parameters
            (firstrow, firstcol, subrows, subcolumns).
            *NOTE: Rows and columns are numbered from 1.*
        total_samples: int, optional
            The total number of times the pixel is sampled during readout.
            If greater than 1 this reduces the Poisson noise.
            Defaults to the current readout mode.
        removeneg: bool, optional, default=True
            If True, remove negative values from the readout and replace
            them with zero.
        
        :Returns:
        
        read_data: array_like uint32
            The data array read out from the detector.
            
        """
        # Default to the number of samples defined by the readout mode.
        if total_samples is None:
            total_samples = self.samplesum * self.nframes
        # >>> Read out the integrator. The integrator returns an integer
        # array, but this needs to be converted to floating point to
        # prevent the readout noise calculation from wrapping around
        # and generating spurious large values in the reference pixel
        # regions.
        read_data = self.pixels.readout(nsamples=total_samples).astype(np.double)

#         # Apply a random drift to the amplifier electronics
#         random_level_change(amplifier_properties['MAX_REF_DRIFT'])
#         
#         # Apply the readout effects from each of the list of amplifiers
#         # associated with this detector array.
#         for amp in self.amplifier_list:
#             read_data = amp.readout(read_data)
#             # Use for debugging only - tedious output.
#             if self._makeplot and self._verbose > 7:
#                 strg = "after readout from amplifier %s" % amp._count
#                 amp.plot(description=strg)
#                 amp.plot()
#  
#         # Check that all areas of the detector have been consistently
#         # read out. This checks that the rows and columns covered by
#         # the amplifiers do not overlap.
#         if not check_readout():
#             strg = "Inconsistent detector readout. " \
#                 "Check the rows/columns assigned to the amplifiers"
#             raise ValueError(strg)

        # Apply gain and readnoise effects. NOTE: Although the Poisson noise
        # and read noise are not added explicitly in quadrature,
        # combining the two sets of randomly generated offsets has the
        # same effect.
        if self.simulate_read_noise and not (self.readnoise_map is None):
 
            # Take a random sample of numbers from a normal distribution
            # with zero mean and unit variance and distribute them over
            # the detector pixels.
            randArray = np.random.randn(read_data.shape[0],
                                        read_data.shape[1]).astype(np.double)
                 
            # Determine the amount of read noise to by applied to the detector
            # data. The noise is reduced when there is more than one sample.
            if total_samples > 1:
                if self._verbose > 3:
                    self.logger.debug( "Applying readout noise in the range of " + \
                        "%f to %f (e), sampled %d times." % \
                        (np.min(self.readnoise_map), np.max(self.readnoise_map),
                         total_samples) )
                noise = self.readnoise_map / np.sqrt(total_samples)
            else:
                if self._verbose > 3:
                    self.logger.debug( "Applying readout noise in the range of %f to %f (e)." % \
                        (np.min(self.readnoise_map), np.max(self.readnoise_map)) )
                noise = self.readnoise_map
            # Multiply by the noise scaling factor.
            # FIXME: The multiplication is needed, even when x 1.0 because it changes the data type.
            if self.noise_factor != 1.0:
                self.logger.debug( "Scaling noise by %f." % self.noise_factor )
            noise = noise * self.noise_factor

            # Assume that bad pixels flagged as UNRELIABLE_SLOPE are noisy.
            # TODO: The multiplication factor is arbitrary.
            if self.bad_pixels is not None:
                noisy_zones = np.where((self.bad_pixels & MASK_UNRELIABLE_SLOPE) > 0)
                noise[noisy_zones] *= 5.0

            read_data = read_data + (randArray * noise)   
            del randArray
         
        # Apply the gain to convert electrons into DNs
        if self.simulate_gain and not (self.gain_map is None):
            if self._verbose > 3:
                self.logger.debug( "Applying readout gain factor of %f to %f (e/DN)." % \
                    (np.min(self.gain_map), np.max(self.gain_map)) )
            read_data = read_data / self.gain_map
         
        # The smallest DN is zero, so the read noise shouldn't make the
        # reading go negative. If requested, replace negative values
        # with zero.
        if removeneg:
            arenegative = np.where(read_data < 0.0)
            read_data[arenegative] = 0.0
                
        # Extract a subarray from the readout if requested and the data
        # isn't already the correct size.
        if subarray is not None:
            if (subarray[2] != self.illuminated_shape[0]) or \
               (subarray[3] != self.illuminated_shape[1]):
                newdata = self._extract_subarray(read_data, subarray)
                read_data = newdata

        # At this point, the DN values may be translated with the linearity
        # table, to simulate the effect of detector nonlinearity.
        if self.simulate_nonlinearity and NONLINEARITY_BY_TABLE:
            rcolumns = read_data.shape[1]
            rcolmiddle = rcolumns//2
            maxleft = len(self.linearity_table_left)
            maxright = len(self.linearity_table_right)
            if self._verbose > 3:
                self.logger.debug( \
                    "Applying nonlinearity lookup tables of length %d and %d." % \
                    (maxleft,maxright) )
            # TODO : Can this tedious lookup be done more quickly using numpy?
            for rrow in range(0,read_data.shape[0]):
                for rcolumn in range(0,rcolmiddle):
                    oldvalue = int(read_data[rrow][rcolumn])
                    if oldvalue < 0:
                        oldvalue = 0
                    if oldvalue > maxleft:
                        oldvalue = maxleft
                    newvalue = self.linearity_table_left[oldvalue]
                    read_data[rrow][rcolumn] = newvalue
                for rcolumn in range(rcolmiddle,rcolumns):
                    oldvalue = int(read_data[rrow][rcolumn])
                    if oldvalue < 0:
                        oldvalue = 0
                    if oldvalue > maxright:
                        oldvalue = maxright
                    newvalue = self.linearity_table_right[oldvalue]
                    read_data[rrow][rcolumn] = newvalue

        if self._verbose > 4:
            self.logger.debug( "Detector readout: min=" + str(read_data.min()) + \
                " max=" + str(read_data.max()) )
                
        # Return the readout data converted back to unsigned integer
        return read_data.astype(np.uint32)

    def _extract_subarray(self, full_data, subarray):
        """
        
        Extract and return a subarray from full frame data.
        
        :Parameters:
        
        full_data: array_like
            Full frame data.
        subarray: tuple of 4 ints
            Contains subarray parameters
            (firstrow, firstcol, subrows, subcolumns).
            *NOTE: Rows and columns are numbered from 1.*
        
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range, or
            if the subarray does not intersect the full data
            properly.
        TypeError
            Raised if any of the parameters are of the wrong type,
            size or shape.
        
        :Returns:
        
        subarray_data: array_like
            A subarray extracted from the data.
        
        """
        # Subarray must be a tuple or list of 4 numbers.
        if isinstance(subarray,(tuple,list)):
            if len(subarray) != 4:
                strg = "Subarray must be described by a tuple of 4 numbers."
                raise TypeError(strg)
        else:
            strg = "Subarray must be described by a tuple of 4 numbers."
            raise TypeError(strg)
            
        # Determine the area of the illuminated pixels extracted by
        # the subarray mode. Note that 1 is subtracted from the subarray
        # start point because subarray indices start at 1, not 0.
        r1 = subarray[0] - 1
        c1 = subarray[1] - 1
        if subarray[2] <= 0 or subarray[3] <= 0:
            raise ValueError("Zero sized subarrays are not possible.")
        r2 = r1 + subarray[2]
        c2 = c1 + subarray[3]
        subarray_illuminated = full_data[r1:r2, c1:c2]
        if subarray_illuminated.size < 1:
            strg = "Subarray [%d:%d, %d:%d] outside data limits (%d x %d)." % \
                (r1,r2,c1,c2,full_data.shape[0],full_data.shape[1])
            raise ValueError(strg)
            
        # When a subarray is read out the amplifier responsible for
        # reading the reference pixels will insert reference columns
        # to the subarray. When these are deinterlaced they will form
        # extra reference pixels at the top and/or bottom of the array
        # (at the top for MIRI level 1 FITS data). The numbers of such
        # rows included is proportional to the size of the subarray.
        self.subrows_bottom = int( subarray[2] * \
                           self.bottom_rows / self.illuminated_shape[0])
        if self.subrows_bottom > 0:
            x1 = 0
            y1 = subarray[1]
            x2 = x1 + self.subrows_bottom
            y2 = y1 + subarray[3]
            subarray_bottom = full_data[x1:x2, y1:y2]
            
        self.subrows_top = int( subarray[2] * \
                           self.top_rows / self.illuminated_shape[0])
        if self.subrows_top > 0:
            x1 = self.illuminated_shape[0]
            y1 = subarray[1]
            x2 = x1 + self.subrows_top
            y2 = y1 + subarray[3]
            subarray_top = full_data[x1:x2, y1:y2]
        
        if self._verbose > 3:
            strg = "Extracting subarray with illuminated shape (%d %d)" % \
                subarray_illuminated.shape
            if self.subrows_bottom > 0:
                strg += " and a bottom reference zone of shape (%d %d)" % \
                    subarray_top.shape
            if self.subrows_top > 0:
                strg += " and a top reference zone of shape (%d %d)" % \
                    subarray_top.shape
            self.logger.debug( strg )
        
        if self.subrows_bottom > 0 and self.subrows_top > 0:
            array_list = (subarray_bottom, subarray_illuminated, subarray_top)
        elif self.subrows_top > 0:
            array_list = (subarray_illuminated, subarray_top)
        elif self.subrows_bottom > 0:
            array_list = (subarray_bottom, subarray_illuminated)
        else:
            # A single array - no concatenation needed
            return subarray_illuminated
        
        try:
            subarray_combined = np.concatenate(array_list, axis=0)
        except (ValueError, IndexError) as e:
            strg = "Subarray must not extend outside detector area "
            strg += "when there are reference pixels\n (%s)" % e
            raise ValueError(strg)
        return subarray_combined
    
    def __str__(self):
        """
        
        Returns a string describing the DetectorArray object.
        
        """
        strg = "DetectorArray: Detector %s (%d). FPM %s. Detector name %s." % \
            (self.detectorid, self.sca_id, self.fpm_name,  self.detector_name)
        strg += "\n  Detector size: %d rows x %d columns of pixels " % \
            self.detector_shape
        strg += "(%d x %d of which are illuminated). " % \
            self.illuminated_shape
        if self.bad_pixels is not None:
            bad_pixel_count = np.sum(self.bad_pixels > 0)
            strg += "%d bad pixels are applied." % bad_pixel_count
        strg += "\n  Detector temperature is %.2f K." % self.temperature
        strg += "\n  " + self.pixels.__str__()
        return strg

    def plot_readout(self, plotfig=None, plotaxis=None, labels=True,
                     withbar=True, title=''):
        """
        
        Plot the detector readout within the given matplotlib axis.
        This function can can be used to include a plot of this object
        in any figure.
        The plotfig method can be used to create a self-contained
        plot.

        :Parameters:
         
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        withbar: bool, optional
            Set to False to suppress the colour bar. The default is True.
        title: string, optional
            Optional title to be shown above the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
        
        :Side effects:
        
        Note that calling plot has the side effect of calling readout(),
        which can alter simulated effects that change with the number of
        readouts made.
            
        """
        # Display the data as an image with the origin at the bottom.
        # NOTE: plotaxis.imshow() will plot the image but will cause the
        # colorbar() function to fail.
        read_data = self.readout()
        if labels:
            xlabel = 'Columns'
            ylabel = 'Rows'
        else:
            xlabel = ''
            ylabel = ''
        plotaxis = mplt.plot_image(read_data, plotfig=plotfig,
                                   plotaxis=plotaxis, withbar=withbar,
                                   xlabel=xlabel, ylabel=ylabel, title=title)
        return plotaxis

    def plot(self, description=''):
        """
        
        Plot the photon count and detector readout in a self-contained
        matplotlib figure.
        
        :Parameters:
        
        description: string, optional
            Additional description to be shown on the plot, if required.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
        
        :Side effects:
        
        Note that calling plot has the side effect of calling readout(),
        which can alter simulated effects that change with the number of
        readouts made.
            
        """
        # Plot the detector status in two subplots: the expected electron
        # count within the ImperfectIntegrator and the actual DNs obtained
        # from a detector readout.
        # Create a figure to contain the plot (larger than default size).
        fig = mplt.new_figure(1, figsize=(10,10), stitle=description)

        ax1 = mplt.add_subplot(fig, subrows=1, subcols=2, subno=1)
        ax1 = self.pixels.plot(plotfig=fig, plotaxis=ax1, withbar=False,
                               title="Electron count before readout")
        
        ax2 = mplt.add_subplot(fig, subrows=1, subcols=2, subno=2)
        ax2 = self.plot_readout(plotfig=fig, plotaxis=ax2, withbar=False,
                                title="Detector readout (DN)")

        if self._verbose > 1:
            strg2 = "Plotting detector data. "
            strg2 += "Close the plot window to continue."
        else:
            strg2 = ''
        mplt.show_plot(prompt=strg2)
            

#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_poisson_integrator module.
# However, the tests here are made in verbose mode and include plotting,
# so they show what is happening in more detail.
#
if __name__ == '__main__':
    print( "Testing the DetectorArray class" )

    # WHICH TESTS TO RUN?
    TEST_FRAME_RTI = True
    TEST_READINGS = True
    TEST_FRAMETIMES = True

    # MODIFY THESE TWO VARIABLES TO CONTROL THE DEGREE OF INTERACTION
    verbose = 3
    PLOTTING = False        # Set to False to turn off plotting.

    NREADINGS = 12
    NFRAGMENTS = 3
    NINTS = 3

    if TEST_FRAME_RTI:
        print( "\n===Frame RTI tests." )
        print("Col  Col  Row  Row  Smp Smp Rst   Brs    Pred")
        print("Strt Stop Strt Stop Skp Sum Rid   Md     (s)")

        namps = 4
        detectorrows = \
                detector_properties.get('_sca494', 'ILLUMINATED_ROWS')
        ampcolumns = 2 + \
                detector_properties.get('_sca494', 'ILLUMINATED_COLUMNS')/namps
        resetoverhead = \
                detector_properties.get('_sca494', 'RESET_OVERHEAD')
        refpixsampleskip = 3
        clock_time = \
                detector_properties.get('_sca494', 'CLOCK_TIME')
        
        test_parameters = [(  1, 258,   1, 1024, 0, 1, 4, False),
                           (  1, 258,   1, 1024, 0, 1, 6, False),
                           (  1,  64,  65,  320, 0, 1, 4, False),
                           (  1,   5,   1,   16, 0, 1, 4, False),
                           (  1,   5,   1,   16, 0, 1, 4, True),
                           (  2,  65,  65,  320, 0, 1, 4, False),
                           ( 91, 106, 765,  828, 0, 1, 4, False),
                           ( 91, 106, 765,  828, 0, 1, 4, True),
                           ( 96,  99, 787,  802, 0, 1, 4, False),
                           ( 96,  99, 787,  802, 0, 1, 4, True),
                           (181, 244,  65,  320, 0, 1, 4, False),
                           (181, 244,  65,  320, 0, 1, 4, True),
                           
                           (  1, 258,   1, 1024, 1, 8, 4, False),
                           (  1, 258,   1, 1024, 2, 8, 4, False),
                           
                           (  1,  16,   1,   64, 0, 1, 1, False),
                           (  1,  16,   1,   64, 0, 1, 4, False),
                           (  1,  16,   1,   64, 1, 1, 4, False),
                           (  1,  16,   1,   64, 2, 8, 4, False),
                           (  1,  16,   1,   64, 2, 16, 4, False),
                            ]
        for (colstart, colstop, rowstart, rowstop, sampleskip, samplesum,
             resetwidth, burst_mode) in test_parameters:
            frame_clks = frame_rti(detectorrows, ampcolumns,
                                   colstart, colstop, rowstart, rowstop,
                                   sampleskip, refpixsampleskip, samplesum,
                                   resetwidth, resetoverhead, burst_mode)
            frame_time = frame_clks * clock_time
    
            strg = "%4d %4d %4d %4d %3d %3d %3d  %5s %9.5f" % (colstart, colstop,
                                                      rowstart, rowstop,
                                                      sampleskip, samplesum,
                                                      resetwidth, str(burst_mode),
                                                      frame_time)
            print(strg)
                    
    if TEST_READINGS or TEST_FRAMETIMES:
        # Create a DetectorArray object for a detector with a the same
        # area of illuminated pixels as SCA 494 (1024 x 1024).
        nrows = detector_properties.get('_sca494', 'ILLUMINATED_ROWS')
        ncolumns = detector_properties.get('_sca494', 'ILLUMINATED_COLUMNS')
        # Raise the detector temperature so the dark current and artefacts
        # are more apparent.
        detector = DetectorArray('MIRIFULONG', nrows, ncolumns, 7.5,
                                 makeplot=PLOTTING, verbose=verbose)
#         detector = DetectorArray('MIRIFULONG', nrows, ncolumns, 6.7,
#                                  makeplot=PLOTTING, verbose=verbose)
     
        samples = detector_properties.get('READOUT_MODE', 'SLOW')
        print("SLOW mode sample array=%s" % str(samples))
        detector.set_readout_mode(samples[0], samples[1], samples[2])
     
        # Display the contents
        print( "===Status after creation of DetectorArray." )
        print( detector )
     
        darkflux = 0.1 * np.ones([nrows,ncolumns])
        mediumflux = 10.0 * np.ones([nrows,ncolumns])
        brightflux = 100.0 * np.ones([nrows,ncolumns])
 
    if TEST_READINGS:
        # Make 4 normal readings
        NREADINGS = 4
        IPERIODS = 3
        detector.reset()
        for reading in range(1,NREADINGS+1):
     
            # Iterate for 3 integration periods
            for count in range(1,IPERIODS+1):
     
                # Integrate with some flux
                detector.integrate(mediumflux, 30.0)
         
                if verbose > 2:
                    print( "\n===Status after integration period", count )
                    print( detector )
        
            # Read out the integrator
            readout = detector.readout()
            if verbose > 1:
                print( "Reading", reading, "signal:\n", readout )
    #         if PLOTTING:
    #             detector.plot(description='Test plot at reading %d' % reading)
     
        # Make another reading in quick succession.
        # Ensure this reading cannot be less than the previous one.
        # NOTE: For a proper test there need to be some assert statements.
        readout = detector.readout()
        if verbose > 2:
            print( "Another reading", reading, "signal:\n", readout )
        if PLOTTING:
            detector.plot(description='Test plot at end of integration')
    #    detector.hit_by_cosmic_rays(100000.0, 512, 512)
    #    detector.plot(description='Test plot after cosmic ray hit')
         
        # Reset the detector and check the signal has returned to zero
        # (within the expected level of persistence).
        detector.reset()
        if verbose > 1:
            print( "\n===Status after a reset." )
            print( detector )
     
        # Integrate on the detector again and ensure the integration
        # starts again at zero (with the expected level of persistence).
        # NOTE: This should perhaps be checked with an assert statement?
        detector.integrate(mediumflux, 30.0)
        if verbose > 2:
            print( "\n===Status after another integration after the reset." )
            print( detector )
        print( "" )

    if TEST_FRAMETIMES:
        fast_samples = detector_properties.get('READOUT_MODE', 'FAST')
        slow_samples = detector_properties.get('READOUT_MODE', 'SLOW')
    
        subarray_list = detector_properties.get('STANDARD_SUBARRAYS')
        for subarray in subarray_list:
            subarray_properties = detector_properties.get('SUBARRAY', subarray)
            for burst_mode in (False, True):
                fast_time = detector.frame_time(fast_samples[0], fast_samples[1],
                                                refpixsampleskip=fast_samples[2],
                                                subarray=subarray_properties,
                                                burst_mode=burst_mode)
                slow_time = detector.frame_time(slow_samples[0], slow_samples[1],
                                                refpixsampleskip=slow_samples[2],
                                                subarray=subarray_properties,
                                                burst_mode=burst_mode)
                fast_exp = detector.exposure_time(2, 1, fast_samples[0],
                                                  fast_samples[1],
                                                  refpixsampleskip=fast_samples[2],
                                                  subarray=subarray_properties,
                                                  burst_mode=burst_mode,
                                                  nframes=1, groupgap=0)
                slow_exp = detector.exposure_time(2, 1, slow_samples[0],
                                                  slow_samples[1],
                                                  refpixsampleskip=slow_samples[2],
                                                  subarray=subarray_properties,
                                                  burst_mode=burst_mode,
                                                  nframes=1, groupgap=0)
    
                if burst_mode:
                    substr = subarray + " (burst)"
                else:
                    substr = subarray
                strg = "Subarray %18s: Frame time: Fast=%6.3fs; Slow=%6.3fs. " % \
                    (substr, fast_time, slow_time)
                strg += "Elapsed for 2 ints: Fast=%6.3fs; Slow=%6.3fs." % \
                    (fast_exp[1], slow_exp[1])
                print( strg )

    print( "Test finished." )
