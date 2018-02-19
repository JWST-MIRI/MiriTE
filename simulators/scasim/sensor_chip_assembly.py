#!/usr/bin/env python

"""

Module sensor_chip_assembly - Contains the SensorChipAssembly class, which
is the top level class of the MIRI Sensor Chip Assembly (SCA) simulator,
plus associated functions.

The Sensor Chip Assembly class is used to convert a detector illumination
input into a FITSWriter or level 1b FITS file output, simulating the
response of the MIRI detectors.

The class hierarchy is:

    SensorChipAssembly
        CosmicRayEnvironment
            CosmicRay
        MiriIlluminationModel
        MiriBadPixelMaskModel
        MiriExposureModel  or  ExposureData
            MiriMeasuredModel
        DetectorArray
            ImperfectIntegrator
                PoissonIntegrator
            MiriMeasurement
        MiriQuantumEfficiency
            MiriFilter
            
The DetectorArray class simulates the underlying detector array, with the
SensorChipAssembly class managing the application interface: interpreting
the input data, managing the metadata and packaging the output data into
the correct format.

NOTE: There are three variations of SensorChipAssembly class, each of which
is based on the Singleton design pattern. Once an instance is created it
will persist until explicitly destroyed. This feature is used to implement
the detector latency simulation, but it does produce side effects. The
constructor is called only once, and a setup() method is used to change the
detector properties once the a SensorChipAssembly object has been created.
    
The following configuration files are use to describe the properties of the
MIRI SCA and its environment:

    cosmic_ray_properties
    detector_properties
    
The simulation is also based on the calibration data contained in the MIRI
bad pixel mask, gain, nonlinearity, dark current and pixel flat-field
Calibration Data Products (CDPs).

:History:

18 Jun 2010: Created.
28 Jun 2010: Reworked to recommended Python syntax standards.
05 Jul 2010: Corrected a typo in the input file variable name in function
             simulate_sca.
07 Jul 2010: Amplifier class added.
15 Jul 2010: Detector and amplifier properties imported from parameter
             files and read noise and dark current determined by
             detector temperature.
20 Jul 2010: Forward the readout mode to the DetectorArray object.
             Plot the exposure data when in verbose mode.
26 Jul 2010: cosmic_ray package included. Detector frame time calculated.
             Readout mode parameters obtained from detector properties.
01 Aug 2010  Major changes to sort out integration time definition.
             More header information written to the output file.
02 Aug 2010: Added ability to average groups and integrations.
             Added detector subarray mode.
06 Aug 2010: Documentation formatting problems corrected.
11 Aug 2010: FITS header updated.
16 Aug 2010: Modified to work with different focal plane modules.
24 Aug 2010: Read detector configuration data from FITS files.
01 Sep 2010: Improved attribute checking for public functions.
13 Sep 2010: Added modulation of the input flux by a fringe map. Flux
             array saved in between integrations to improve efficiency.
             Software version string passed in rather than defined here.
14 Sep 2010: Added flags to turn Poisson noise or read noise on and off.
01 Oct 2010: Dropped groups parameter added to readout mode.
07 Oct 2010: FITS header keywords updated after reading "Definition
             of the MIRI FM ILT Keywords" by T.W.Grundy.
11 Oct 2010: Added flags to turn dark current or amplifier effects on
             and off. Added makeplot flag.
12 Oct 2010: Number of frames and dropped frames included in detector
             readout mode and used in exposure time calculation.
             FITSWriter option added for saving exposure data.
18 Oct 2010: Check for no QE filter. Added qe_adjust flag.
05 Nov 2010: Take left and right columns, bottom and top rows and well
             depth from detector properties. Pass number of frames per
             group to cosmic ray functions. Subarray mode read from
             illumination data header.
16 Nov 2010: Added function to set the random number generator seed, for
             testing. Added simulate_bad_pixels flag.
23 Nov 2010: Some bugs corrected in the processing of 3-D illumination
             data. Input of subarray data now supported.
26 Nov 2010: Handling of subarray modes considerably simplified. Any
             combination of input and output modes is now supported.
15 Dec 2010: Integration and group commentary added.
07 Jan 2011: Detector and amplifier properties classified and looked up
             by SCA_ID rather than by FPM_ID.
26 Jan 2011: MiriMetadata object introduced to manage the FITS header
             keywords in the exposure function.
18 Feb 2011: Updated to match change of group name from DQ to DQL in
             mirikeyword.py.
04 Mar 2011: Documentation tweaks to resolve bad formatting.
10 Mar 2011: FLT keyword group imported.
15 Mar 2011: set_data method added so that scasim can be called directly
             from other Python software. QuantumEfficiency class moved to
             miritools.
22 Mar 2011: Plotting functions modified to use matplotlib figures
             and axes more flexibly (unfortunately the colorbar
             function doesn't work properly). QuantumEfficiency class
             plotting functions changed.
30 Mar 2011: Very large detector flux levels rejected before they
             cause a problem in PoissonIntegrator. simulate_sca_fromdata
             function provided as an alternative to simulate_sca.
             Bugs corrected in the plot_ramp and set_data functions.
             Added simulate_sca_fromdata function. Code duplication
             removed. Default values for parameters will now be taken
             from the header of the input file, if provided.
01 Apr 2011: Documentation formatting problems corrected.
20 Apr 2011: Intensity scaling factor added - for debugging.
21 Apr 2011: Set overwrite=True in to_fits_header function to ensure
             existing header keywords in the input data are
             overwritten with new values. TSAMPLE should be the clock
             time in microseconds, not sample time in seconds. COLSTART,
             COLSTOP, ROWSTART and ROWSTOP keywords added to FITS header.
             Integrations executed before assembling metadata.
03 May 2011: Added the ability to parse a non-standard input subarray
             mode.
05 May 2011: Additional subarray checks. Subarray default adjusted so
             the output mode is never non-standard.
14 Jul 2011: Use of exceptions documented.
14 Sep 2011: Exposure data written in uint16 format rather than float32.
14 Sep 2011: Modified to keep up with the separation of
             miritools.miridata into miritools.metadata and
             miritools.miricombdata
29 Sep 2011: _prepare_for_new_data and _new_detector functions modified
             to allow persistence effects to be simulated.
             simulate_sca_list function added so that multiple exposures
             can be simulated in one go, with persistence effects passed
             from one exposure to the next.
             Extra resets are inserted between integrations when the
             FRAME_RESETS parameter is greater than zero.
20 Sep 2011: Destructor added to remove large objects. Integration data
             explicitly deleted when finished with (rather than relying
             on Python garbage collector). Reporting of the reading of
             the illumination map now happens in data_maps.
25 Oct 2011: cosmic_ray_properties, detector_properties and
             amplifier_properties imported using ParameterFileManager.
             This allows new parameter files to be substituted by putting
             new versions into the current working directory, without the
             need to rebuild the source code.
02 Nov 2011: Deprecated use of the .has_key() dictionary method removed.
11 Jan 2012: Renamed symbols to avoid potential name clash: file-->thisfile
             and format-->fileformat.
04 Mar 2012: add_to_header replaced with to_fits_header
22 Mar 2012: Metadata group import now specifies a wildcard.
29 Mar 2012: QuantumEfficiency class now imported from miritools.filters
             module.
13 Apr 2012: Major changes to implement the Filter and QuantumEfficiency
             classes as a JwstDataProduct.
26 Apr 2012: Brought up to date with changes in Metadata class.
21 May 2012: Plotting modified to use the miriplot utilities.
13 Nov 2012: Major restructuring of the package folder. Import statements
             updated.
14 Jan 2013: olddataproduct subpackage split from dataproduct.
05 Jun 2013: Tidied up constructor and moved bad pixel mask and dark map
             access to separate functions. Rationalised attribute names.
10 Jun 2013: Modified to use the new MIRI exposure model as an alternative.
             Replaced the old MiriMetadata model with a simpler one in
             data_maps. Metadata keywords and FITS header access functions
             changed to match the new data model.
18 Jun 2013: Added ability to specify the output file format and switch
             between old and new models.
04 Jul 2013: Old commented out code removed.
19 Aug 2013: Changed subarray keywords to match the change made in
             jwst_lib (SUBXSTRT --> SUBSTRT1, etc...).
25 Oct 2013: Make the primary keyword for determining the detector ID
             'DETECTOR' rather than 'SCA_ID'. Pass detectorid as a
             parameter instead of scaid.
23 Apr 2014: Some global attributes renamed.
10 Jun 2014: Allow simulation of reference pixels to be turned off.
             Corrected problem passing "verbose" parameter to illumination
             model. Explicitly remove the exposure object at the beginning
             of each new exposure.
17 Jun 2014: Improved memory management. Do not delete the detector object
             whenever there is a change of data unless the shape has changed.
             Changed verbosity definitions and some log messages. Added the
             ability to override detector properties and define the detector
             frame time explicitly. Removed the old "preflash" latency
             implementation.
20 Jun 2014: Replaced 'SUBMODE' keyword with 'SUBARRAY' and 'READOUT'
             keyword with 'READPATT'.
             Detector data converted to uint32 format rather than uint16
             format, to prevent data being truncated when
             simulate_gain=False.
25 Jun 2014: fileformat and simulate_ref_pixels typos corrected in the
             global sca_ functions.
11 Jul 2014: Switch over to using the new MiriQuantumEfficiency class,
             based on the jwst_lib data models. Modified documentation.
18 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
02 Mar 2015: Added subarray burst mode as an option.
05 Mar 2015: Trap a ValueError as well as a KeyError when attempting to set
             a metadata keyword. Default output data shape changed to
             hypercube.
06 Mar 2015: Modified the selection of output data format and shape.
20 May 2015: Reimplemented SensorChipAssembly class as a singleton.
21 May 2015: Amplifier level removed and gain and read noise defined by
             calibration data products.
10 Jun 2015: Obsolete calls to fits_history and fits_comment functions
             replaced by calls to generic functions.
05 Aug 2015: Missing metadata added to data created by the
             simulate_sca_pipeline function. More realistic exposure start
             and end times included in the metadata.
14 Aug 2015: Simulation flag setting functions combined into one.
             Inconsistency in name of bad pixel flag corrected.
19 Aug 2015: Removed MiriCubeModel.
02 Sep 2015: "slope" option added to output data shape.
             Added tests for subarray modes. Corrected bug which caused
             the previous PIXELDQ array to be remembered when there was
             no bad pixel mask.
04 Sep 2015: Corrected a bug in the sizing of the PIXELDQ array for subarrays.
08 Sep 2015: Legacy illumination map and fringe map models removed.
             Made compatible with Python 3.
17 Sep 2015: Declare three separate singleton classes, so the caller can
             simulate 1, 2 or 3 independent detectors.
22 Sep 2015: Improved the management of the clock time during an exposure.
             When there are subsequent exposures, the start time of each
             exposure is, by default, set of to end time of the previous
             exposure. The start time can also be defined explicitly.
             Fixed a problem where the "setup" method erased the detector
             object and spoiled the singleton behaviour.
25 Sep 2015: simulate_sca and simulate_sca_list global functions replaced
             by simulate_files method. simulate_sca_pipeline global function
             replaced by simulate_pipe method. simulate_sca_fromdata method
             removed. Cosmic ray library only reloaded when the cosmic ray
             is changed.
06 Oct 2015: Attribute names modified to be compatible with documentation.
             Typo corrected in get_exposure_times function.
07 Oct 2015: Made exception catching Python 3 compatible.
03 Nov 2015: Display the exposure times when testing.
18 Nov 2015: Over-sized illumination maps truncated to match maximum
             detector illumination size.
04 Dec 2015: Added a logger.
23 Feb 2016: Corrected bug in FITSWriter output log string.
25 Feb 2016: Cosmic ray simulation suffix parameter added, to choose between
             different variations of cosmic ray simulation.
10 Mar 2016: Bad pixel simulation moved to Detector object. Do not copy
             'units' metadata from illumination model to exposure model
             data object.
23 Mar 2016: Obtain MASK, PIXELFLAT and GAIN calibrations from CDP files.
             NOTE: This generates more error messages, because CDP files
             are not available for all combinations of readout modes.
30 Mar 2016: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
31 Mar 2016: SUBSTRT1 and SUBSTRT2 metadata parameters changed to start at 1.
             SUBSIZE1 metadata parameter changed to include reference pixels.
03 May 2016: Added the REFOUT extension needed by the STScI pipeline
             to the level 1b format exposure data.
04 May 2016: Set the EXP_TYPE keyword if not already defined. Made the
             creation of the GROUPDQ array in the output data optional.
05 May 2016: Log strings shortened and subarray passed to DetectorArray object.
             simulate_amp_effects split into simulate_gain and
             simulate_nonlinearity. Shortened metadata comments.
06 May 2016: Use the pixel flat-field appropriate for each subarray mode
             (and filter or band, if appropriate).
20 May 2016: Do not overwrite the data type metadata in the output ramp data.
             Increased the number of groups in the output test data.
06 Jun 2016: Added some missing simulate_ref_pixels and simulate_nonlinearity
             flags. Added flags to allow the version numbers of imported CDPs
             to be specified. Subarray burst mode made True by default.
             Improved documentation.
05 Jul 2016: CDP-6 dark maps added. Integration number added to dark map
             selection. Dark maps are now used unnormalised, and the dark
             current vs temperature measurement is instead normalised at 6.7K.
14 Jul 2016: Added simulate_drifts and simulate_latency flags. Transfer the
             NON-SCIENCE data quality flag to the PIXELDQ array.
19 Jul 2016: Shift World Coordinates metadata to account for reference
             columns and subarray mode.
03 Aug 2016: Corrected a bug in simulate_pipe() where the input subarray mode
             was mixed up with the output subarray mode. Renamed the
             set_illumination() function parameters to make this clearer.
05 Aug 2016: Reset the reference output array size when the subarray mode is
             FULL. Shuffle the list of subarrays to improve test coverage.
05 Sep 2016: Test data for subarrays should specify detector MIRIMAGE.
28 Sep 2016: miri.miritools.dataproduct renamed miri.datamodels and
             miri.miritools renamed miri.tools.
02 Nov 2016: Bug corrected in simulate_files: only save the simulation to
             an output file when the file name is not empty. Check for
             uppercase or lowercase data shape.
10 Jan 2017: If the subarray mode changes, only read a new flat-field CDP
             if flat-field simulation is switched on. Filter out NaN
             values before using illumination data. Corrected bug in flux
             warning messages. Added more subarray tests.
07 Feb 2017: include_dq flag split into include_pixeldq and include_groupdq,
             so that the pixeldq and groupdq arrays can be separately
             excluded from the output exposure data.
14 Feb 2017: Corrected the subarray to FULL simulation tests.
13 Jun 2017: Added cdp_ftp_path parameter.
16 Jun 2017: Modified the extreme saturation, "junk data" check so that
             extreme data are truncated, rather than throwing an exception.
             Corrected mistake in default setting of ftp_path.
07 Jul 2017: Added the ability to add a full dark calibration (but do not
             activate it).
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Jul 2017: Removed the 'IPC' variant of the cosmic ray library.
18 Jul 2017: Ability to add the full dark calibration debugged and activated.
             Bug in the application of the pixel flat-field to subarray data
             corrected. Fixed "False," bug in the transfer of simulator
             flags to detector.add_calibration_data() function.
27 Jul 2017: Report the number of cosmic ray events in the metadata.
             Added an option to the wait() method to simulate a non-zero
             background flux while idling.
05 Sep 2017: Added clear_exposure_data() method.
11 Sep 2017: Start and end times must be communicated to the exposure
             data model as MJD rather than seconds.
12 Sep 2017: Correct the WCS info for ramp data. Look for WCS information
             in the intensity data of the illumination model and write
             WCS information to 'SCI' HDU of the exposure data.
             Corrected the WCS correction for the reference columns, including
             for subarray data which does not touch the left edge.
21 Sep 2017: Additional WCS and WCS-like keywords copied from the INTENSITY
             metadata of the input illumination model.
03 Oct 2017: Added test simulations of all known MRS bands.
13 Oct 2017: New frame time calculation from Mike Ressler.
             SLOW mode now uses 8 out of 9 samples. READOUT_MODE now defines
             samplesum, sampleskip and refpixsampleskip parameters separately.
01 Nov 2017: Do not simulate detector drifts and latents in SLOW mode, since
             the measured parameters are only valid for FAST mode.
08 Jan 2018: Import version number from miri package.
14 Feb 2018: Added version number for nonlinearity CDP.
19 Feb 2018: Non-linearity correction moved here.

@author: Steven Beard

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function
from astropy.extern import six

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

#import warnings
import time, os
# The simulator uses the mathematical and array processing functions of
# numpy.
import numpy as np

from miri import __version__

# General purpose quantum efficiency class
from miri.datamodels.miri_filters import MiriQuantumEfficiency

# MIRI data models
from miri.datamodels import MiriMeasuredModel

# Old and new MIRI exposure data models. The old model is compatible with DHAS
# and the new model is compatible with the STScI pipeline.
from miri.simulators.scasim.exposure_data import ExposureData, Metadata, \
    get_file_header
from miri.datamodels.sim import MiriExposureModel, \
    MiriIlluminationModel

# The following modules form part of the sca simulator.
from miri.simulators.scasim.cosmic_ray import CosmicRayEnvironment, CosmicRay, \
    load_cosmic_ray_library, load_cosmic_ray_random
from miri.simulators.scasim.detector import DetectorArray, SIM_CDP_FTP_PATH, \
    NONLINEARITY_BY_TABLE

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Search for the cosmic ray and detector parameters files and parse them
# into properties dictionaries. The files are searched for in 3 places:
# (a) The current directory
# (b) The directory where this Python file is being executed from
# (c) The miri.simulators.scasim installation directory.
from miri.tools.filesearching import ParameterFileManager, make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
cosmic_ray_properties = ParameterFileManager(
                            "cosmic_ray_properties.py",
                            search_path=search_path,
                            description="cosmic ray properties",
                            logger=LOGGER)
detector_properties = ParameterFileManager(
                            "detector_properties.py",
                            search_path=search_path,
                            description="detector properties",
                            logger=LOGGER)

#
# Global helper function
#
def _report_sca_parameters(inputfile, outputfile, detectorid, readout_mode,
                           fringemap, fileformat, datashape, qe_adjust,
                           simulate_poisson_noise, simulate_read_noise,
                           simulate_ref_pixels, simulate_bad_pixels,
                           simulate_dark_current, simulate_flat_field,
                           simulate_gain, simulate_nonlinearity,
                           simulate_drifts, simulate_latency,
                           cdp_ftp_path,
                           readnoise_version, bad_pixels_version,
                           flat_field_version, linearity_version, gain_version,
                           logger=LOGGER):
    """
    
    Helper function to report the SCA simulation parameters.
    (Saves duplication of code between simulate_sca and
    simulate_sca_fromdata.)
    
    """
    if inputfile:
        if isinstance(inputfile, (tuple,list)):
            strg = "Simulating " + str(readout_mode) + \
                " detector readout for " + str(detectorid) + \
                " from list of illumination files:"
            for thisfile in inputfile:
                strg += "   %s" % thisfile
            logger.info( strg )
        elif isinstance(inputfile, MiriIlluminationModel):
            detectorid = ''
            if hasattr(inputfile, 'meta'):
                if hasattr(inputfile.meta, 'instrument'):
                    if hasattr(inputfile.meta.instrument, 'detector'):
                        detectorid = inputfile.meta.instrument.detector
            detectorid = inputfile.meta.instrument.detector
            readout_mode = ''
            if hasattr(inputfile, 'meta'):
                if hasattr(inputfile.meta, 'exposure'):
                    if hasattr(inputfile.meta.instrument, 'readpatt'):
                        readout_mode = inputfile.meta.exposure.readpatt
            logger.info( "Simulating " + str(readout_mode) + \
                " detector readout for " + str(detectorid) + \
                " from illumination data model of shape " + \
                str(inputfile.intensity.shape) + ".")
            logger.debug( inputfile.get_meta_str() )
        else:
            logger.info( "Simulating " + str(readout_mode) + \
                " detector readout for " + str(detectorid) + \
                " from illumination file:\n   %s" % inputfile )
    else:
        logger.info( "Simulating " + str(readout_mode) + \
            " detector readout for " + str(detectorid) + \
            " from illumination data array" )

    if fringemap is not None and fringemap:
        logger.info( "Simulation modulated by a fringe map contained in:\n   %s" % fringemap )

    if isinstance(outputfile,(tuple,list)):
        strg = "Results will be written to a list of " + str(fileformat) + \
            " FITS files with " + str(datashape) + " data shape:"
        for thisfile in outputfile:
            if thisfile:
                strg += "\n   %s" % thisfile
            else:
                strg += "\n   (no output)"
        logger.info( strg )
    elif outputfile:
        logger.info( "Results will be written to a " + str(fileformat) + \
                     " FITS file with " + \
                     str(datashape) + "data shape:\n   %s" % outputfile )
    else:
        logger.info( "Results will be returned in an exposure data model." )
        
    # NOTE: This duplicates "Simulation control flags" message
    switched_off = []
    if not qe_adjust:
        switched_off.append("Quantum efficiency")
    if not simulate_poisson_noise:
        switched_off.append("Poisson noise")
    if not simulate_read_noise:
        switched_off.append("Read noise")
    if not simulate_ref_pixels:
        switched_off.append("Reference pixels")
    if not simulate_bad_pixels:
        switched_off.append("Bad pixels")
    if not simulate_dark_current:
        switched_off.append("Dark current")
    if not simulate_flat_field:
        switched_off.append("Pixel flat-field")
    if not simulate_gain:
        switched_off.append("Bias and Gain")
    if not simulate_nonlinearity:
        switched_off.append("Non-linearity")
    if ('SLOW' in readout_mode) or (not simulate_drifts):
        switched_off.append("Detector drifts")
    if ('SLOW' in readout_mode) or (not simulate_latency):
        switched_off.append("Detector latency")
    if switched_off:
        logger.debug( "NOTE: The following effects are switched off: %s" % \
            "\'" + "\', \'".join(switched_off) + "\'" )

    if cdp_ftp_path:
        logger.debug("CDP search path folder: " + \
                     str(cdp_ftp_path))

    if readnoise_version:
        logger.debug("Read noise CDP restricted to version number: " + \
                     str(readnoise_version))
    if bad_pixels_version:
        logger.debug("Bad pixels CDP restricted to version number: " + \
                     str(bad_pixels_version))
    if flat_field_version:
        logger.debug("Pixel flat-field CDP restricted to version number: " + \
                     str(flat_field_version))
    if linearity_version:
        logger.debug("Nonlinearity CDP restricted to version number: " + \
                     str(linearity_version))
    if gain_version:
        logger.debug("Gain CDP restricted to version number: " + \
                     str(gain_version))

def _get_parameter_defaults(fpm, metadata, readout_mode, subarray, frame_time,
                            temperature, cosmic_ray_mode, verbose=2,
                            logger=LOGGER):
    """
    
    Helper function to obtain appropriate defaults for parameters
    that have not been explicitly set.
    (Saves duplication of code between simulate_sca and
    simulate_sca_fromdata.)
    
    """
    # If the readout mode is not specified, obtained a default from
    # the FITS metadata of the input file. Failing that, obtain a
    # default value from the detector properties.
    if readout_mode is None:
        if metadata is not None and 'READPATT' in metadata:
            readout_mode = metadata['READPATT']
            if verbose > 2:
                logger.info( "Readout mode %s obtained from FITS metadata." % \
                    readout_mode )
        else:
            readout_mode = detector_properties['DEFAULT_READOUT_MODE']
            if verbose > 2:
                logger.info( "Readout mode defaulted to " + \
                    "%s from detector properties." % readout_mode )
    else:
        if verbose > 2:
            logger.info( "Readout mode explicitly set to %s." % readout_mode )

    # If the output subarray mode is not specified, obtained a default
    # from the FITS metadata of the input file, as long as this is a
    # known subarray mode. Failing that, obtain a default value from
    # the detector properties.
    if subarray is None:
        if metadata is not None and 'SUBARRAY' in metadata:
            subarray = metadata['SUBARRAY']
            if subarray in detector_properties['SUBARRAY']:
                if verbose > 2:
                    logger.info( "Subarray mode %s obtained from FITS metadata." % \
                        subarray )
            else:
                nonstandard = subarray
                subarray = detector_properties['DEFAULT_SUBARRAY']
                if verbose > 2:
                    strg = "Subarray mode %s obtained from FITS metadata " % \
                        nonstandard
                    strg += "is non-standard, so output subarray mode "
                    strg += "defaulted to %s from detector properties." % \
                        subarray
                    logger.info( strg )
        else:
            subarray = detector_properties['DEFAULT_SUBARRAY']
            if verbose > 2:
                logger.info( "Subarray mode defaulted to " + \
                    "%s from detector properties." % subarray )
    else:
        if verbose > 2:
            logger.info( "Subarray mode explicitly set to %s." % subarray )

    if frame_time is None:
        if metadata is not None and 'TFRAME' in metadata:
            frame_time = metadata['TFRAME']
            if verbose > 2:
                strg = "Frame time of %f seconds obtained " % \
                    frame_time
                strg += "from FITS metadata "
                strg += "(overriding the readout mode and subarray)."
                logger.info( strg )
    else:
        strg = "Frame time of %f seconds specified explicitly " % frame_time
        strg += "(overriding the readout mode and subarray)."
        logger.info( strg )

    # If the detector temperature is not specified, use the target
    # temperature from the detector properties.
    if temperature is None:
        temperature = fpm['TARGET_TEMPERATURE']
        if verbose > 3:
            logger.debug( "Temperature defaulted to %fK from detector properties." % \
                temperature )
    else:
        if verbose > 2:
            logger.info( "Temperature explicitly set to %fK." % temperature )

    # If the cosmic ray mode is not specified, obtained a default from
    # the FITS metadata of the input file. Failing that, obtain a default
    # value from the cosmic ray properties.
    if cosmic_ray_mode is None:
        if metadata is not None and 'CRMODE' in metadata:
            cosmic_ray_mode = metadata['CRMODE']
            if verbose > 3:
                logger.debug( "Cosmic ray mode %s obtained from FITS metadata." % \
                    cosmic_ray_mode )
        else:
            cosmic_ray_mode = cosmic_ray_properties['DEFAULT_CR_MODE']
            if verbose > 3:
                logger.debug( "Cosmic ray mode defaulted to " + \
                    "%s from cosmic ray properties." % cosmic_ray_mode )
    else:
        if verbose > 2:
            logger.info( "Cosmic ray mode explicitly set to %s." % cosmic_ray_mode )
            
    return (readout_mode, subarray, frame_time, temperature, cosmic_ray_mode)

def _clock_to_mjd( clock_seconds ):
    """
    
    Helper function which converts clock time in seconds to MJD in days
    
    """
    # Modified Julian date of the "zero epoch" of the system clock (1/1/70)
    MJD_ZEROPOINT = 40587.0
    # Number of seconds per day.
    SECONDS_PER_DAY = 86400.0
    mjd_days = MJD_ZEROPOINT + (clock_seconds/SECONDS_PER_DAY)
    return mjd_days

def _mjd_to_clock( mjd_days ):
    """
    
    Helper function which converts MJD in days to clock time in seconds.
    
    """
    # Modified Julian date of the "zero epoch" of the system clock (1/1/70)
    MJD_ZEROPOINT = 40587.0
    # Number of seconds per day.
    SECONDS_PER_DAY = 86400.0
    clock_seconds = (mjd_days - MJD_ZEROPOINT) * SECONDS_PER_DAY
    return clock_seconds
 

# A metaclass which only allows one instance of a class to exist at a time.
class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = \
                super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class SensorChipAssembly(object):
    """

    Class Sensor Chip Assembly - Simulates the behaviour of the MIRI
    detectors. The class converts detector illumination data (typically
    provided in a FITS file generated by another MIRI simulator and read
    by the read_data method) and generates simulated detector ramp data,
    which may be written to a FITSWriter or level 1b FITS file by the
    write_data method.

    :Parameters:
    
    See setup method.   
                
    """
    
    def __init__(self, logger=LOGGER):
        """
        
        Constructor for class SensorChipAssembly.
        
        :Parameters:
        
        logger: Logger object (optional)
            A Python logger to handle the I/O. This parameter can be used
            by a caller to direct the output to a different logger, if
            the default defined by this module is not suitable.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("scasim")
        self._setup = False
        self._makeplot = False
        self._verbose = 0
        self.cosmic_ray_mode = None
        self.cosmic_ray_env = None
        self.detectorid = None
        self._sca = None
        self.detector = None
        self.shape = ()
        self.clock_time = time.time()
        
        # Initial values for the simulation flags
        self.simulate_poisson_noise = True
        self.simulate_read_noise = True
        self.simulate_ref_pixels = True
        self.simulate_bad_pixels = True
        self.simulate_dark_current = True
        self.simulate_flat_field = True
        self.simulate_gain = True
        self.simulate_nonlinearity = True
        self.simulate_drifts = True
        self.simulate_latency = True
      
    def setup(self, detectorid, readout_mode='FAST', subarray='FULL',
              burst_mode=True, inttime=None, ngroups=None, nints=None,
              temperature=6.7, cosmic_ray_mode='SOLAR_MIN',
              fileformat='STScI', include_pixeldq=True, include_groupdq=False,
              qe_adjust=True, simulate_poisson_noise=True,
              simulate_read_noise=True, simulate_ref_pixels=True,
              simulate_bad_pixels=True, simulate_dark_current=True,
              simulate_flat_field=True, simulate_gain=True,
              simulate_nonlinearity=True,
              simulate_drifts=True, simulate_latency=True,
              cdp_ftp_path=SIM_CDP_FTP_PATH,
              readnoise_version='', bad_pixels_version='',
              flat_field_version='', linearity_version='', gain_version='',
              makeplot=False, verbose=2):
        """
        
        Set up a SensorChipAssembly simulation.
        
        :Parameters:
        
        detectorid: string
            The detector ID, identifying a particular detector.
            The MIRI instrument has three detectors: 'MIRIMAGE',
            'MIRIFULONG' and 'MIRIFUSHORT'. (These correspond to the imager
            and the long and short wavelength arms of the MRS respectively.)
            readout_mode: string, optional, default='FAST'
        Readout mode, which can be one of the following common options:
        
            * 'SLOW' - 10 samples per readout and defaults of ngroups=10
              and nints=1  
            * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
              nints=10
            * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
              averaged to reduce data volume.
            * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
              averaged to reduce data volume.
          
            The following unusual options are also available for testing:
        
            * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
              4 frames per group and a gap of 8 frames between each group.
            * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
              averaged and default ngroups=1.
            * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
              averaged and default ngroups=4.
            * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
              4 frames per group and a gap of 8 frames between each group.
        
        subarray: string, optional, default='FULL'
            Detector subarray mode for output. This can be one 'FULL',
            'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
            'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
            full-frame and the other modes read out portions of the detector
            as described in the MIRI Operational Concept Document.
            The simulator can also accept the other subarray modes defined
            in detector_properties.py for test purposes.
        burst_mode: boolean, optional, default is True
            True if a subarray is being read out in burst mode,
            which skips quickly over unwanted columns.
        inttime: float, optional
            The integration time in seconds to be simulated. This parameter
            will be ignored if ngroups is provided as well.
            *NOTE: Any requested integration time will be rounded up to
            the time resulting from the nearest whole number of groups.*
        ngroups: int, optional
            The number of readout groups. If None, this is derived from
            the integration time and readout mode. Otherwise the parameter
            overrides the integration time and sets the number of groups
            directly. There must be at least one group.
        nints: int, optional
            The number of integrations per exposure. It must be at least 1.
            If None, the default set by the readout mode is used.
            The total exposure time is nints x inttime.
        temperature: float, optional, default=6.7
            The temperature of the detector, which will determine the dark
            current and read noise.
        cosmic_ray_mode: string, optional, default='SOLAR_MIN'
            The cosmic ray environment mode to be simulated. Available modes
            are:
        
            * 'NONE' - No cosmic rays.
            * 'SOLAR_MIN' - Solar minimum
            * 'SOLAR_MAX' - Solar maximum
            * 'SOLAR_FLARE' - Solar flare (worst case scenario)
            
        fileformat: string, optional, default='STScI'
            The kind of file format to be written.
            
            * 'STScI' - use the STScI level 1 model for the JWST
              DMS pipeline.
            * 'FITSWriter' - emulate the format written by the FITSWriter
              during MIRI VM, FM and CV tests and read by the DHAS.
        
        include_pixeldq: boolean, optional, default=True
            A flag that may be used to switch on the inclusion of the
            PIXELDQ data array in the output exposure data.
            By default, this array is included (and contains the bad pixel
            mask used to generate the simulated data).
        include_groupdq: boolean, optional, default=False
            A flag that may be used to switch on the inclusion of the
            GROUPDQ data array in the output exposure data.
            By default, this array is not included.
        qe_adjust: boolean, optional, default=True
            A flag that may be used to switch off quantum efficiency
            adjustment (for example to observe what effects in a simulation
            are caused by QE). *NOTE: When QE adjustment is turned off,
            the input illumination is assumed to be in electrons/second.*
        simulate_poisson_noise: boolean, optional, default=True
            A flag that may be used to switch off Poisson noise (for example
            to observe what effects in a simulation are caused by Poisson
            noise).
        simulate_read_noise: boolean, optional, default=True
            A flag that may be used to switch off read noise (for example
            to observe what effects in a simulation are caused by read
            noise).
        simulate_ref_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of reference
            pixels in the data. 
        simulate_bad_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of bad
            pixels in the data, even if a bad pixel map containing bad pixels
            is specified in the detector properties.
        simulate_dark_current: boolean, optional, default=True
            A flag that may be used to switch off the addition of dark
            current (for example to observe what effects in a simulation are
            caused by dark current).
        simulate_flat_field: boolean, optional, default=True
            A flag that may be used to switch off the simulation of the pixel
            flat-field (for example to observe the effect of quantum
            efficiency in isolation).
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
        cdp_ftp_path: str, optional, default=None
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
        # Note that all parameters (which don't have None as an option)
        # are explicitly converted into the expected data type, since
        # the values extracted from the properties dictionary can
        # sometimes be converted by Python into strings, which then
        # upsets the formatted I/O statements.
        self._makeplot = bool(makeplot)
        self._verbose = int(verbose)
        if verbose > 3:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("+++SensorChipAssembly object setup for" + \
                " detector=" + str(detectorid) + \
                " readout mode=" + str(readout_mode) + \
                " subarray mode=" + str(subarray) + \
                " burst mode=" +  str(burst_mode) + \
                " inttime=" + str(inttime) + \
                " ngroups=" + str(ngroups) + \
                " nints=" + str(nints) + \
                " temperature=" + str(temperature) + "K" + \
                " cosmic ray mode=" + str(cosmic_ray_mode) )
                
        # Extract the properties of the particular sensor chip assembly
        # being simulated.
        if self.detectorid != detectorid:
            self.detectorid = detectorid
            # Erase any existing detector object when the detector ID changes
            if self.detector is not None:
                del self.detector
                self.detector = None
        try:
            self._sca = detector_properties.get('DETECTORS_DICT',
                                                self.detectorid)
        except KeyError:
            strg = "%s is not a known detector ID." % self.detectorid
            raise KeyError(strg)

        # Set the starting readout mode. This function also checks the
        # validity of the inttime, ngroups and nints attributes and initialises
        # the collection of attributes describing the readout mode.
        self.set_readout_mode(readout_mode, inttime=inttime, ngroups=ngroups,
                              nints=nints)
        self.readout_mode_previous = ''
        
        # Initialise the output and input subarray modes.
        self.subarray_str = subarray
        self.subarray_previous = ''
        try:
            self.subarray = detector_properties.get('SUBARRAY', subarray)
        except KeyError:
            strg = "Unrecognised or badly defined output subarray mode: %s" % \
                subarray
            raise ValueError(strg)
        # The input subarray mode is not known until an illumination map
        # is provided.
        self.subarray_input_str = 'FULL'
        self.subarray_input = None
        self.subarray_burst_mode = burst_mode
                    
        if verbose > 1:
            self.logger.info( self.readout_mode_str(prefix='  ') )
        
        # There is no illumination map, flux data, fringe map or fringe data
        # until they are read in and calculated.
        self.illumination_map = None
        self.flux = None
        self.fringe_map = None
        self.fringe_map_data = None
        
        # There is no metadata until the illumination map is read.
        self.metadata = None
        
        # There is no integration data or exposure_data object until an
        # integration.
        self.integration_data = None
        self.exposure_data = None
        self.include_pixeldq = include_pixeldq
        self.include_groupdq = include_groupdq
        self.fileformat = fileformat
        
        try:
            self.temperature = float(temperature)
        except (TypeError, ValueError) as e:
            strg = "SCA temperature must be a valid number."
            strg += "\n %s" % e
            raise ValueError(strg)
        if verbose > 1:
            self.logger.info( self.temperature_str() )
        
        if qe_adjust:
            if verbose > 2:
                self.logger.info( "Loading QE from \`%s\`" % self._sca['QE_FILE'] )
            # Load the QE information. Plot it if it isn't null.
            #self.qe = QuantumEfficiency(self._sca['QE_FILE'])
            self.qe = MiriQuantumEfficiency(self._sca['QE_FILE'])
            if makeplot and self.qe is not None:
                self.qe.plot()
        else:
            # No QE adjustment
            self.qe = None
        
        self.define_cosmic_ray_env(cosmic_ray_mode)
            
        self.set_simulation_flags(readout_mode, qe_adjust=qe_adjust,
                                  simulate_poisson_noise=simulate_poisson_noise,
                                  simulate_read_noise=simulate_read_noise,
                                  simulate_ref_pixels=simulate_ref_pixels,
                                  simulate_bad_pixels=simulate_bad_pixels,
                                  simulate_dark_current=simulate_dark_current,
                                  simulate_flat_field=simulate_flat_field,
                                  simulate_gain=simulate_gain,
                                  simulate_nonlinearity=simulate_nonlinearity,
                                  simulate_drifts=simulate_drifts,
                                  simulate_latency=simulate_latency)

        self.cdp_ftp_path       = cdp_ftp_path
        self.readnoise_version  = readnoise_version
        self.bad_pixels_version = bad_pixels_version
        self.flat_field_version = flat_field_version
        self.linearity_version  = linearity_version
        self.gain_version       = gain_version

        self._setup = True
        
    def __del__(self):
        """
        
        Destructor for class SensorChipAssembly.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        if self._setup:
            try:
                # Objects created during setup
                if self.qe is not None:
                    del self.qe
                if self.cosmic_ray_env is not None:
                    del self.cosmic_ray_env
            
                # Objects created by class methods
            
                if self.detector is not None:
                    del self.detector
                if self.illumination_map is not None:
                    del self.illumination_map
                if self.metadata is not None:
                    del self.metadata
                if self.fringe_map_data is not None:
                    del self.fringe_map_data
                if self.fringe_map is not None:
                    del self.fringe_map
                if self.exposure_data is not None:
                    del self.exposure_data
            except Exception:
                pass

    def define_cosmic_ray_env(self, cosmic_ray_mode):
        """
        
        Define the cosmic ray environment for the simulation.
        
        Cosmic ray simulations for the MIRI sensors are assumed to be
        contained in files named CRs_SiAs_35_<mode>_<number>_<suffix>.fits
        where
        
        * <mode> is the cosmic ray mode.
        * <number> is a simulation run number (selected at random by
          this function).

        :Parameters:
        
        cosmic_ray_mode: string, optional, default='SOLAR_MIN'
            The cosmic ray environment mode to be simulated. Available modes
            are:
        
            * 'NONE' - No cosmic rays.
            * 'SOLAR_MIN' - Solar minimum
            * 'SOLAR_MAX' - Solar maximum
            * 'SOLAR_FLARE' - Solar flare (worst case scenario)
        
        """
        # Only set up a new cosmic ray environment when the mode has changed,
        # or if the environment doesn't exist.
        if (cosmic_ray_mode == self.cosmic_ray_mode) and \
           (self.cosmic_ray_env is not None):
            if self._verbose > 1:
                self.logger.info( self.cosmic_ray_str() + " (No change.)" )
            return
        
        # Set up the cosmic ray environment, based on the required cosmic
        # ray mode. Choose one of the 10 library files available for that
        # mode at random.
        self.cosmic_ray_mode = cosmic_ray_mode
        if self._verbose > 1:
            self.logger.info( self.cosmic_ray_str() )
        wordz = self.cosmic_ray_mode.split('+')
        cr_mode = wordz[0]
        if len(wordz) > 1:
            suffix = wordz[1]
        else:
            suffix = ''
        if cr_mode != 'NONE':
            filenum = int(np.random.uniform(0, 9))
            try:
                crname = cosmic_ray_properties.get('CR_LIBRARY_FILES', \
                                                   cr_mode)
                if crname is not None:
                    if suffix:
                        filename = '%s0%d_%s.fits' % (crname, filenum, suffix)
                    else:
                        filename = '%s0%d.fits' % (crname, filenum)
                else:
                    filename = None
            except KeyError:
                strg = "Unrecognised cosmic ray mode: %s" % cosmic_ray_mode
                raise KeyError(strg)
                        
            # Delete any previous cosmic ray environment
            if self.cosmic_ray_env is not None:
                del self.cosmic_ray_env
            try:
                self.cosmic_ray_env = load_cosmic_ray_library(filename,
                                                    cosmic_ray_mode=cr_mode,
                                                    convolve_ipc=True,
                                                    verbose=self._verbose,
                                                    logger=self.toplogger)
            except Exception:
                # If the cosmic ray library file cannot be read, fall back
                # to using a random distribution.
                if filename is not None:
                    strg = "***Cosmic ray library '%s' not readable." % \
                        filename
                else:
                    strg = "***Cosmic ray library is not defined."
                strg += " Falling back to RANDOM mode."
                self.logger.warn(strg)
                self.cosmic_ray_env = load_cosmic_ray_random(
                                                    cosmic_ray_mode=cr_mode,
                                                    verbose=self._verbose,
                                                    logger=self.toplogger)
        else:
            self.cosmic_ray_env = load_cosmic_ray_random(
                                                    cosmic_ray_mode=cr_mode,
                                                    verbose=self._verbose)
     
    def set_seed(self, seedvalue=None):
        """
        
        Set the seed for the numpy random number generator.
        This function can be used while testing to ensure the
        randomised effects that follow are well defined.
        
        :Parameters:
        
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator.
            If not specified, a value of None will be sent, which
            randomises the seed.
            
        """
        np.random.seed(seedvalue)

    def set_simulation_flags(self, readout_mode, qe_adjust=True,
                             simulate_poisson_noise=True,
                             simulate_read_noise=True,
                             simulate_ref_pixels=True,
                             simulate_bad_pixels=True,
                             simulate_dark_current=True,
                             simulate_flat_field=True,
                             simulate_gain=True,
                             simulate_nonlinearity=True,
                             simulate_drifts=True,
                             simulate_latency=True):
        """
        
        Set the simulation control flags.
        
        """
        if self._verbose > 1:
            strg = "Simulation control flags:"
            strg += "\n\tQuantum efficiency simulation turned "
            if qe_adjust:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tPoisson noise simulation turned "
            if simulate_poisson_noise:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tRead noise simulation turned "
            if simulate_read_noise:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tReference pixels simulation turned "
            if simulate_ref_pixels:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tBad pixels simulation turned "
            if simulate_bad_pixels:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tDark current simulation turned "
            if simulate_dark_current:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tFlat-field simulation turned "
            if simulate_flat_field:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tAmplifier bias and gain turned "
            if simulate_gain:
                strg += "ON."
            else:
                strg += "OFF."
            strg += "\n\tDetector non-linearity effects turned "
            if simulate_nonlinearity:
                strg += "ON."
            else:
                strg += "OFF."
            if 'FAST' in readout_mode:
                strg += "\n\tDetector drift effects turned "
                if simulate_drifts:
                    strg += "ON."
                else:
                    strg += "OFF."
                strg += "\n\tDetector latency effects turned "
                if simulate_latency:
                    strg += "ON."
                else:
                    strg += "OFF."
            else:
                strg += "\n\tDetector drift and latency effects not "
                strg += "simulated for SLOW readout mode."
                simulate_drifts = False
                simulate_latency = False
            self.logger.info( strg )
        # If any of the flags have changed, delete the detector and warn
        # the user
        if self.simulate_poisson_noise != simulate_poisson_noise or \
           self.simulate_read_noise != simulate_read_noise or \
           self.simulate_ref_pixels != simulate_ref_pixels or \
           self.simulate_bad_pixels != simulate_bad_pixels or \
           self.simulate_dark_current != simulate_dark_current or \
           self.simulate_flat_field != simulate_flat_field or \
           self.simulate_gain != simulate_gain or \
           self.simulate_nonlinearity != simulate_nonlinearity or \
           self.simulate_drifts != simulate_drifts or \
           self.simulate_latency != simulate_latency:
            if self.detector is not None:
                self.logger.warn("Simulation flags have changed. " + \
                        "Restarting simulation with a new detector. " + \
                        "Latent effects up to this point will be erased.")
                del self.detector
                self.detector = None
        # Now change the flags
        self.simulate_poisson_noise = simulate_poisson_noise
        self.simulate_read_noise = simulate_read_noise
        self.simulate_ref_pixels = simulate_ref_pixels
        self.simulate_bad_pixels = simulate_bad_pixels
        self.simulate_dark_current = simulate_dark_current
        self.simulate_flat_field = simulate_flat_field
        self.simulate_gain = simulate_gain
        self.simulate_nonlinearity = simulate_nonlinearity
        self.simulate_drifts = simulate_drifts
        self.simulate_latency = simulate_latency
        
    def read_data(self, filename, scale=1.0 ):
        """
        
        Read the detector illumination data from a file.
                
        :Parameters:
     
        filename: string
            The name of the file to be opened.
        scale: float, optional, default=1.0
            An optional scale factor to apply to the intensity data.
            This is for debugging and testing only - the input data
            should already be scaled to photons/second/pixel.
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
        TypeError
            Raised if any of the parameters are of the wrong type,
            size or shape.
                      
        """
        if not self._setup:
            strg = "No simulation parameters defined. Please use the setup() "
            strg += " method before defining input data."
            raise AttributeError(strg)
        
        # Tidy up and prepare for new data.
        self._prepare_for_new_data()

        # Create the illumination map object and read the data from the file.        
        self.illumination_map = MiriIlluminationModel(filename)
        
        # The intensity data must be at least 2-D.
        if self.illumination_map.intensity.ndim < 2:
            strg = "Illumination map must be at least 2-D "
            strg += "(%d-D data provided)" % self.illumination_map.intensity.ndim
            raise TypeError(strg)
        
        # The illumination map cannot be larger than the detector.
        # Give a warning and truncate an over-sized map.
        illum_shape = self.illumination_map.get_illumination_shape()
        if illum_shape[0] > self._sca['ILLUMINATED_COLUMNS'] or \
           illum_shape[1] > self._sca['ILLUMINATED_ROWS']:
            strg = "***Illumination map of size %d x %d is too large! " % \
                illum_shape
            strg += "Truncating to detector size of %d x %d pixels." % \
                (self._sca['ILLUMINATED_COLUMNS'], self._sca['ILLUMINATED_ROWS'])
            self.logger.warn(strg)
            self.illumination_map.truncate([self._sca['ILLUMINATED_COLUMNS'],
                                            self._sca['ILLUMINATED_ROWS']] )

        if scale is not None and scale != 1.0:
            self.illumination_map.apply_scale(scale)
        if self.illumination_map.meta.subarray.name and \
           self.illumination_map.meta.subarray.name != "GENERIC":
            self.subarray_input_str = self.illumination_map.meta.subarray.name
        else:
            # If the input subarray mode is undefined or GENERIC but the
            # output subarray mode matches the input data shape, assume the
            # input and output subarray modes are the same. Otherwise assume
            # full frame.
            if self.subarray is not None:
                ishape = self.illumination_map.get_illumination_shape()
                if (self.subarray[2] == ishape[0]) and \
                   (self.subarray[3] == ishape[1]):
                    self.subarray_input_str = self.subarray_str
                    if self._verbose > 2:
                        self.logger.info( "Input subarray mode assumed to be %s." % \
                            self.subarray_input_str )
                else:
                    self.subarray_input_str = 'FULL'
            else:
                self.subarray_input_str = 'FULL'

        # Get the subarray location and dimensions from the detector
        # properties, or failing that attempt to parse the input subarray
        # string into a list of 4 integers. Any blank or comma separated
        # list is allowed, with or without surrounding brackets.
        if self.subarray_input_str in detector_properties['SUBARRAY']:
            self.subarray_input = \
                detector_properties.get('SUBARRAY',self.subarray_input_str)
        else:
            try:
                self.subarray_input = eval(self.subarray_input_str)
            except Exception:
                words = self.subarray_input_str.split()
                try:
                    firstrow = int(words[0].strip(','))
                    firstcol = int(words[1].strip(','))
                    subrows = int(words[2].strip(','))
                    subcols = int(words[3].strip(','))
                    self.subarray_input = (firstrow, firstcol,
                                            subrows, subcols)
                except Exception:
                    strg = "Unrecognised subarray mode specified "
                    strg += "in input file: %s" % self.subarray_input_str
                    raise ValueError(strg)

        # A subarray definition must contain exactly 4 values.
        if isinstance(self.subarray_input, (list,tuple)) and \
           (len(self.subarray_input) != 4):
            strg = "Subarray mode %s declared in input file " % \
                self.subarray_input_str
            strg += "translates to subarray list of wrong length: %s" % \
                self.subarray_input.__str__()
            raise TypeError(strg)

        if self._verbose > 2:
            if self.subarray_input is None:
                strg = "Input subarray mode assumed %s." % \
                    self.subarray_input_str
            else:
                strg = "Input subarray mode is %s " % self.subarray_input_str
                strg += "(%d %d %d %d)." % self.subarray_input
            self.logger.info( strg )
        
        # Check the integrity of the subarray declared in the input data.
        ishape = self.illumination_map.get_illumination_shape()
        if self.subarray_input is not None:
            # The location must be within the expected detector limits.
            drows = self._sca['ILLUMINATED_ROWS']
            dcols = self._sca['ILLUMINATED_COLUMNS'] + \
                        self._sca['LEFT_COLUMNS'] + self._sca['RIGHT_COLUMNS']
            if (self.subarray_input[0] < 1) or \
               (self.subarray_input[0] + self.subarray_input[2] > dcols) or \
               (self.subarray_input[1] < 1) or \
               (self.subarray_input[1] + self.subarray_input[3] > drows):
                strg = "Subarray mode %s declared in input file " % \
                    self.subarray_input_str
                strg += "translates to invalid location: %s" % \
                    self.subarray_input.__str__()
                raise ValueError(strg)
            
            # The data must have the same shape and size as the subarray.
            if (ishape[0] != self.subarray_input[2]) or \
               (ishape[1] != self.subarray_input[3]):
                strg = "Subarray mode %s declared in input file " % \
                    self.subarray_input_str
                strg += "is incompatible with the actual data shape "
                strg += "(%d x %d)" % ishape
                raise ValueError(strg)
        else:
            # The input file is assumed full frame. Warn if the data
            # array size might exceed the size of the calibration frames.
            fullcolumns = self._sca['LEFT_COLUMNS'] + \
                            self._sca['ILLUMINATED_COLUMNS'] + \
                            self._sca['RIGHT_COLUMNS']
            fullframe = (self._sca['ILLUMINATED_ROWS'], fullcolumns)
            if (ishape[0] > fullframe[0]) or (ishape[1] > fullframe[1]):
                if self.simulate_bad_pixels or self.simulate_dark_current:
                    strg = "WARNING: Input file data shape (%d x %d) " % ishape
                    strg += "may be too large for the bad pixel or dark "
                    strg += "frames, which expect up to %d x %d. " % fullframe
                    strg += "Try --nobadpixels --nodark."
                    self.logger.warn( strg )
                    
        # Give a warning if the input subarray mode is smaller than the
        # output subarray mode.
        if self.subarray is not None:
            oshape = (self.subarray[2], self.subarray[3])
            if (oshape[0] > ishape[0]) or (oshape[1] > ishape[1]):
                if self._verbose > 0:
                    strg = "NOTE: Output subarray (%d x %d) " % oshape
                    strg += "is larger than input array (%d x %d)." % ishape
                    strg += " There will be some gaps in the output data."
                    self.logger.warn( strg )
       
        # The previous flux data is invalid.
        self.flux = None

        # Create a new detector object matching the size of the
        # data contained in the illumination map.
        self._new_detector()

    def set_illumination(self, illumination_map, scale=None,
                         subarray_input=None ):
        """
        
        Define the detector illumination by providing an MiriIlluminationModel
        object directly. This function is useful when calling scasim
        directly from Python functions that implement a pipeline.
                
        :Parameters:
    
        illumination_map: MiriIlluminationModel object
            A MIRI illumination map object which describes the input
            illumination.
        scale: float, optional, default=None
            An optional scale factor to apply to the intensity data.
            If not provided, no scaling is applied.
            This is for debugging and testing only - the input data
            should already be scaled to photons/second/pixel.
        subarray_input: string, optional
            The subarray mode of the INPUT data. If not specified, the
            subarray mode will be taken from the illumination map
            metadata, and if neither are present it will be assumed
            full frame.
            
        :Raises:
    
        AttributeError
            Raised if simulation parameters undefined.
        TypeError
            Raised if the illumination map is of the wrong type,
            size or shape.

        """
        if not self._setup:
            strg = "No simulation parameters defined. Please use the setup() "
            strg += " method before defining input data."
            raise AttributeError(strg)

        assert isinstance(illumination_map, MiriIlluminationModel)
        intensity = illumination_map.intensity
        wavelength = illumination_map.wavelength
        
        # The intensity data must be at least 2-D.
        if intensity.ndim < 2:
            strg = "Illumination map must be at least 2-D "
            strg += "(%d-D data provided)" % intensity.ndim
            raise TypeError(strg)
        
        # If the wavelength array is 1-D and the intensity array is 3-D
        # the wavelength array needs to be reshaped into a 3-D array.
        if illumination_map._isvalid(intensity) and \
           illumination_map._isvalid(wavelength):
            if len(wavelength.shape) == 1 and len(intensity.shape) > 2:
                if self._verbose > 2:
                    self.logger.info( "Reshaping wavelength array from 1-D to 3-D." )
                sz1 = illumination_map.wavelength.shape[0]
                illumination_map.wavelength.shape = [sz1, 1, 1]

        # The illumination data size must be a multiple of 8 so the
        # reference output can be reshaped successfully.
        # The illumination map cannot be larger than the detector.
        # Give a warning and truncate an over-sized map.
        illum_shape = illumination_map.get_illumination_shape()
        if (illum_shape[0] % 8 != 0) or (illum_shape[0] % 8 != 0):
            strg = "Illumination map rows and columns must be a multiple of 8 "
            strg += "(%d x %d data provided)" % illum_shape
            raise TypeError(strg)
            
        if illum_shape[0] > self._sca['ILLUMINATED_COLUMNS'] or \
           illum_shape[1] > self._sca['ILLUMINATED_ROWS']:
            strg = "***Illumination map of size %d x %d is too large! " % \
                illum_shape
            strg += "Truncating to detector size of %d x %d pixels." % \
                (self._sca['ILLUMINATED_COLUMNS'], self._sca['ILLUMINATED_ROWS'])
            self.logger.warn(strg)
            illumination_map.truncate([self._sca['ILLUMINATED_COLUMNS'],
                                       self._sca['ILLUMINATED_ROWS']] )
            
        # Scale the intensity data if needed.
        if scale is not None and scale != 1.0:
            illumination_map.apply_scale(scale)

        # Use the subarray mode supplied or get the subarray mode from
        # the supplied metadata.
        if subarray_input is None:
            if hasattr(illumination_map, 'meta'):
                if hasattr(illumination_map.meta, 'subarray'):
                    if hasattr(illumination_map.meta.subarray, 'name'):
                        subarray_input = illumination_map.meta.subarray.name
        if subarray_input is not None:
            self.subarray_input_str = str(subarray_input)
        else:
            # If the input subarray mode is undefined but the output subarray
            # mode matches the input data shape, assume the input and output
            # subarray modes are the same. Otherwise assume full frame.
            if self.subarray is not None:
                ishape = illumination_map.get_illumination_shape()
                if (self.subarray[2] == ishape[0]) and \
                   (self.subarray[3] == ishape[1]):
                    self.subarray_input_str = self.subarray_str
                    if self._verbose > 2:
                        self.logger.info( "Input subarray mode assumed to be %s." % \
                            self.subarray_input_str )
                else:
                    self.subarray_input_str = 'FULL'
            else:
                self.subarray_input_str = 'FULL'
       
        # Define the illumination map and subarray mode.
        self._set_illumination_map(illumination_map, subarray_input)

    def _set_illumination_map(self, illumination_map, subarray_input=None ):
        """
        
        Test function to set the illumination map directly from a
        pre-defined illumination_map object, as an alternative to
        reading a file.
                
        :Parameters:
     
        illumination_map: IlluminationMap object
            An illumination map to be used.
        subarray_input: string, optional
            The subarray mode of the INPUT data. If not specified, full
            frame will be assumed.
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
         
        """
        # Verify the input subarray mode provided
        if subarray_input is not None:
            self.subarray_input_str = subarray_input
            try:
                self.subarray_input = \
                    detector_properties.get('SUBARRAY', subarray_input)
            except KeyError:
                # N.B. This function does not attempt to parse an
                # unrecognised input subarray string. Just report an error.
                strg = "Unrecognised input subarray mode specified: %s" % \
                    self.subarray_input_str
                raise ValueError(strg)
        else:
            self.subarray_input_str = 'FULL'
            self.subarray_input = None
        
        # Tidy up and prepare for new data.
        self._prepare_for_new_data()
        
        # Set the illumination map object.        
        self.illumination_map = illumination_map
        
        # The previous flux data is invalid.
        self.flux = None

        # Create a new detector object matching the size of the
        # data contained in the illumination map.
        self._new_detector()

    def _prepare_for_new_data(self):
        """
        
        Helper function to tidy up in preparation for new data.
        
        """            
        # Delete any previous illumination map object
        if self.illumination_map is not None:
            del self.illumination_map
            self.illumination_map = None
            
        # Delete any previous MIRI metadata object and begin with a fresh one.
        if self.metadata is not None:
            del self.metadata
        self.metadata = Metadata(description='Metadata created by SCA simulator')

    def _new_detector(self):
        """
        
        Helper function to set up a detector object for new data.
        
        """
        # All subarray combinations are supported.
        #
        # If the input data is full frame, create a new detector object
        # matching the shape and size of the data contained in the
        # illumination map. Otherwise create full-sized detector object
        # into which the input subarray will be inserted.
        if self.subarray_input is None:
            ishape = self.illumination_map.get_illumination_shape()
        else:
            ishape = (self._sca['ILLUMINATED_ROWS'], \
                      self._sca['ILLUMINATED_COLUMNS'])

        # Query the filter name or band from the illumination model.
        mirifilter = self.illumination_map.get_fits_keyword('FILTER')
        miriband = self.illumination_map.get_fits_keyword('BAND')
            
        # A new detector object is only needed if the detector ID
        # or the shape of the illumination map has changed.
        if (self.detector is None) or (ishape != self.shape):
            # The shape has changed.
            self.shape = ishape

            if self._verbose > 1:
                self.logger.info( "Creating a new detector object for " + \
                    str(self.shape[0]) + " rows x " + \
                    str(self.shape[1]) + " columns." )

            if self.simulate_ref_pixels:
                left_columns  = self._sca['LEFT_COLUMNS']
                right_columns = self._sca['RIGHT_COLUMNS']
                bottom_rows   = self._sca['BOTTOM_ROWS']
                top_rows      = self._sca['TOP_ROWS']
            else:
                left_columns = 0
                right_columns = 0
                bottom_rows = 0
                top_rows = 0
            well_depth = self._sca['WELL_DEPTH']
        
            # Clear the old detector object and create a new one.
            if self.detector is not None:
                del self.detector
            if 'SLOW' in self.readout_mode:
                readpatt = 'SLOW'
            else:
                readpatt = 'FAST'
            self.detector = DetectorArray(
                                self.detectorid,
                                self.shape[0], self.shape[1],
                                self.temperature,
                                left_columns=left_columns,
                                right_columns=right_columns,
                                bottom_rows=bottom_rows,
                                top_rows=top_rows,
                                well_depth=well_depth,
                                readpatt=readpatt, subarray=self.subarray_str,
                                mirifilter=mirifilter, miriband=miriband,
                                simulate_poisson_noise=self.simulate_poisson_noise,
                                simulate_read_noise=self.simulate_read_noise,
                                simulate_bad_pixels=self.simulate_bad_pixels,
                                simulate_dark_current=self.simulate_dark_current,
                                simulate_flat_field=self.simulate_flat_field,
                                simulate_gain=self.simulate_gain,
                                simulate_nonlinearity=self.simulate_nonlinearity,
                                simulate_drifts=self.simulate_drifts,
                                simulate_latency=self.simulate_latency,
                                cdp_ftp_path=self.cdp_ftp_path,
                                readnoise_version=self.readnoise_version,
                                bad_pixels_version=self.bad_pixels_version,
                                flat_field_version=self.flat_field_version,
                                linearity_version=self.linearity_version,
                                gain_version=self.gain_version,
                                makeplot=self._makeplot,
                                verbose=self._verbose,
                                logger=self.toplogger)
            
        # If a new detector is not needed, but the readout mode or the
        # output subarray mode have changed, then the calibration data
        # needs to be updated.
        elif (self.subarray_str != self.subarray_previous) or \
             (self.readout_mode != self.readout_mode_previous):
            if 'SLOW' in self.readout_mode:
                readpatt = 'SLOW'
            else:
                readpatt = 'FAST'
            self.detector.simulate_poisson_noise=self.simulate_poisson_noise
            self.detector.simulate_read_noise=self.simulate_read_noise
            self.detector.simulate_bad_pixels=self.simulate_bad_pixels
            self.detector.simulate_dark_current=self.simulate_dark_current
            self.detector.simulate_flat_field=self.simulate_flat_field
            self.detector.simulate_gain=self.simulate_gain
            self.detector.simulate_nonlinearity=self.simulate_nonlinearity
            self.detector.simulate_drifts=self.simulate_drifts
            self.detector.simulate_latency=self.simulate_latency
            self.detector.add_calibration_data(self.detectorid,
                    readpatt=readpatt, subarray=self.subarray_str,
                    mirifilter=mirifilter, miriband=miriband,
                    cdp_ftp_path=self.cdp_ftp_path,
                    bad_pixels_version=self.bad_pixels_version,
                    flat_field_version=self.flat_field_version,
                    linearity_version=self.linearity_version,
                    readnoise_version=self.readnoise_version,
                    gain_version=self.gain_version)

        self.subarray_previous = self.subarray_str
        self.readout_mode_previous = self.readout_mode
 
    def read_fringe_map(self, filename):
        """
        
        Read fringe map data from a file.
                
        :Parameters:
     
        filename: string
            The name of the file to be opened.
             
        :Prerequisites:
        
        Can only be used after the detector illumination data have been
        imported by the read_data method.
        
        :Raises:
        
        AttributeError
            Raised if the read_data method has not been executed first.
        
        """
        if self._verbose > 1:
            self.logger.info( "Reading fringe map file: %s" % filename )
        if self.illumination_map is None:
            strg = "IlluminationMap object not defined - use read_data first."
            raise AttributeError(strg)

        # The fringe map must be the same shape as the illumination map
        if self.subarray_input is None:
            ishape = self.illumination_map.get_illumination_shape()
        else:
            ishape = (self._sca['ILLUMINATED_ROWS'], \
                      self._sca['ILLUMINATED_COLUMNS'])

        self.fringe_map = filename
        fringe_model = MiriMeasuredModel( filename )
        
        # Get a windowed fringe map scaled so the mean value is 1.0
        fringe_map = fringe_model.data[0:ishape[0], 0:ishape[1]]
        fringe_mean = np.mean(fringe_map)
        if fringe_mean > 0.0:
            self.fringe_map_data = fringe_map / fringe_mean
        else:
            self.fringe_map_data = fringe_map / fringe_mean
            
        # The flux must be recalculated if the fringe map changes.
        self.flux = None
        del fringe_model
        
    def set_readout_mode(self, mode, inttime=None, ngroups=None, nints=None,
                         samplesum=None):
        """
        
        Defines the SCA detector readout mode and integration parameters.
        
        :Parameters:
        
        mode: string
            Readout mode, which can be one of the following options.
            
            * 'SLOW' - 9 samples per readout and defaults of ngroups=10
              and nints=1
            * 'FAST' - 1 sample per readout and defaults of ngroups=1
              and nints=10
            * 'FASTINTAVG' - same as FAST but with groups of 4
              integrations averaged to reduce data volume.
            * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
              averaged to reduce data volume.
              
        inttime: float, optional
            The integration time in seconds. This parameter also defines
            the number of readout groups, unless ngroups is specified.
            The integration time must be positive, as there must be at
            least one group.
            *NOTE: Any requested integration time will be rounded up to
            the time resulting from the nearest whole number of groups.*
        ngroups: int, optional
            The number of readout groups. If None, this is derived from
            the integration time and readout mode. If specified, this
            parameter defines the number of groups and integration time
            and overrides the inttime parameter. There must be at least
            one group.
        nints: int, optional
            The number of integrations required. If None, the default
            for the readout mode will be used. There must be at least
            one integration.
            *NOTE: nints can be safely changed without affecting the
            readout mode.*
        samplesum: int, optional
            Normally None, but if specified allows samplesum to be defined
            directly - overriding the values defined in the detector
            properties file. This is only for testing weird readout
            modes, since MIRI users cannot change samplesum explicitly.
            There must be at least one sample.
            *NOTE: Altering this value from its default will change
            the readout mode to 'TESTnn'.*
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
            
        """
        # Define default parameters based on the readout mode.
        self.readout_mode = mode
        try:
            mode_tuple = detector_properties.get('READOUT_MODE', mode)
            self.samplesum_def = int(mode_tuple[0])
            self.sampleskip_def = int(mode_tuple[1])
            self.refpixsampleskip_def = int(mode_tuple[2])
            self.nframes = int(mode_tuple[3])
            self.groupgap = int(mode_tuple[4])
            self.ngroups_def = int(mode_tuple[5])
            self.nints_def = int(mode_tuple[6])
            self.avggrps = int(mode_tuple[7])
            self.avgints = int(mode_tuple[8])
        except (KeyError, IndexError, TypeError, ValueError):
            strg = "Unrecognised or badly defined readout mode: %s" % mode
            raise ValueError(strg)

        # Determine how the actual detector readout parameters are
        # defined.
        if ngroups is not None:
            if int(ngroups) <= 0:
                raise ValueError("Number of groups must be > 0")
            # Integration time is determined by ngroups.
            # The inttime parameter is ignored.
            self.time_mode = 'groups_to_time'
            self.inttime = None
            self.inttime_req = inttime
            # If avggrps is not 1 round the specified number of groups
            # up to a whole multiple of avggrps.
            if self.avggrps > 1:
                self.ngroups = \
                    int(np.ceil(ngroups/self.avggrps) * self.avggrps)
            else:
                self.ngroups = int(ngroups)
        elif inttime is not None:
            # ngroups is determined by the integration time.
            self.time_mode = 'time_to_groups'
            self.inttime = float(inttime)
            if self.inttime <= 0.0:
                raise ValueError("Integration time must be positive.")
            self.inttime_req = None
            self.ngroups = None
        else:
            # If neither inttime nor ngroups are specified, the defaults
            # set by the readout mode are used.
            self.time_mode = 'default'
            self.inttime = None
            self.inttime_req = None
            self.ngroups = self.ngroups_def

        # If a number of integrations has been provided use it.
        # Otherwise use the read mode default.
        if nints is not None:
            if int(nints) <= 0:
                raise ValueError("Number of integrations must be > 0")
            # If avgints is not 1 round the specified number of integrations
            # up to a whole multiple of avgints.
            if self.avgints > 1:
                self.nints = int(np.ceil(nints/self.avgints) * self.avgints)
            else:
                self.nints = int(nints)
        else:
            self.nints = self.nints_def

        # Overriding nsample is for theoretical testing only
        # and will cause the readout mode to become 'TESTnn'
        if samplesum is None:
            self.samplesum = self.samplesum_def
            self.sampleskip = self.sampleskip_def
            self.refpixsampleskip = self.refpixsampleskip_def
        else:
            if int(samplesum) <= 0:
                raise ValueError("Number of samples must be > 0")
            self.samplesum = int(samplesum)
            if self.samplesum > 2:
                self.sampleskip = 1
            else:
                self.sampleskip = 0
            self.refpixsampleskip = self.refpixsampleskip_def
            self.readout_mode = 'TEST%d' % self.samplesum
            
        # This flag will trigger a recalculation of the integration time
        # during the next integration after a DetectorArray object has been
        # created. (The final integration time depends on the detector
        # frame time determined by the DetectorArray object).
        self.inttime_calculated = False
        
        if self._verbose > 2:
            strg = "Setting readout mode to %s (sampleskip=%d, samplesum=%d): " % \
                (mode, self.sampleskip, self.samplesum)
            if self.time_mode == 'groups_to_time':
                strg += "\n  integration time determined from ngroups=%d" % \
                    self.ngroups
            elif self.time_mode == 'time_to_groups':
                strg += "\n  ngroups determined from integration time=%.2f" % \
                    self.inttime
            else:
                strg += "\n  ngroups=%d from mode default" % self.ngroups
                
            if nints is not None:
                strg += " and nints specified to be %d." % self.nints
            else:
                strg += " and nints defaulting to %d." % self.nints
            self.logger.info( strg )

    def wait(self, elapsed_time, bgflux=0.0):
        """
        
        Simulate waiting or a certain elapsed time.
                
        :Parameters:
        
        elapsed_time: float
            The elapsed time, in seconds.
        bgflux: float (optional)
            A uniform background photon flux falling on the detector
            while it is idle. This flux can affect detector persistence
            by filling charge traps. By default there is no background
            flux.
            
        :Returns:
        
        None.
            
        :Prerequisites:
        
        The detector readout mode and integration time should have been
        defined with the set_readout_mode method.
        
        :Raises:
        
        AttributeError
            Raised if the read_data method has not been executed first.
        ValueError
            Raised if any of the parameters are out of range.
        
        """
        # A wait is not possible until a DetectorArray object
        # has been successfully defined.
        if self.detector is None:
            strg = "DetectorArray object not defined - " + \
                "use read_data or set_illumination first."
            raise AttributeError(strg)

        # Elapsed must be positive.
        assert elapsed_time > 0.0

        if self._verbose > 2:
            self.logger.info( "Waiting for %.2f seconds." % elapsed_time )
        
        # Generate some cosmic ray events during this elapsed time.
        rows = self.shape[0]
        columns = self.shape[1]
        pixsize = self._sca['PIXEL_SIZE']
            
        cosmic_ray_list = \
            self.cosmic_ray_env.generate_events(rows, columns,
                                                elapsed_time, pixsize)
        # Hit the detector with the cosmic rays.
        self.detector.hit_by_cosmic_rays(cosmic_ray_list, self.nframes)
        del cosmic_ray_list    
            
        # Wait for the elapsed time.
        self.detector.wait(elapsed_time, bgflux=bgflux)
        self.clock_time += elapsed_time
          
    def integration(self, intnum, frame_time=None):
        """
        
        Simulate a detector integration consisting of one or more groups
        of readouts.
                
        :Parameters:
        
        intnum: int
            Integration number, which determines where data are
            stored in the simulated science data structure.
            Must be 0 or greater.
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter overrides the frame time
            obtained from the current detector readout mode and
            subarray. This can be used, for example, to simulate
            the timing of full frame exposures with subarray-sized
            data (to save memory). If None, the frame time obtained
            from the current readout mode and subarray is used.
            
        :Returns:
        
        integration_data: array_like uint32
            An array of integration data.
            
        :Prerequisites:
        
        Can only be used after the detector illumination data have been
        imported by the read_data method.
        The detector readout mode and integration time should have been
        defined with the set_readout_mode method.
        
        :Raises:
        
        AttributeError
            Raised if the read_data method has not been executed first.
        ValueError
            Raised if any of the parameters are out of range.
        
        """
        # An integration is not possible until a DetectorArray object
        # and illumination map have been successfully defined.
        if self.illumination_map is None:
            strg = "IlluminationMap object not defined - use read_data first."
            raise AttributeError(strg)
        if self.detector is None:
            strg = "DetectorArray object not defined - use read_data first."
            raise AttributeError(strg)
        
        # The integration number cannot be negative.
        if int(intnum) < 0:
            strg = "The integration number must be a positive integer."
            raise ValueError(strg)

        # Frame time must be positive.
        if frame_time is not None:
            assert frame_time > 0.0
        
        # Define the number of groups from the readout mode and integration
        # time provided, unless specified explicitly. There must be at least
        # one group.
        if not self.inttime_calculated:
            if self.time_mode == 'time_to_groups':
                self.ngroups = self._time_to_groups(
                                        self.inttime,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode,
                                        frame_time=frame_time)
            else:
                self.inttime = self._groups_to_time(
                                        self.ngroups,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode,
                                        frame_time=frame_time)
                self.inttime_req = self.inttime
            self.inttime_calculated = True

#         # Forward the readout mode to the detector.
#         self.detector.set_readout_mode(self.samplesum, self.sampleskip,
#                                        self.refpixsampleskip, self.nframes)
       
        # Get the combined illumination from the illumination map.
        # This calculation only needs to be done once unless the
        # illumination map or the QE measurement changes.
        if self.flux is None:
            if self.subarray_input is None:
                flux = self.illumination_map.get_illumination( usefilter=self.qe )
            else:
                dshape = self.detector.illuminated_shape
                location = (self.subarray_input[0], self.subarray_input[1])
                dleft  = self._sca['LEFT_COLUMNS']
                dright = self._sca['RIGHT_COLUMNS']                
                flux = self.illumination_map.get_illumination_enlarged(
                            dshape, usefilter=self.qe, location=location,
                            leftcrop=dleft, rightcrop=dright)

            # Check the amount of flux is within a sensible range.
            # Give a warning if the flux is likely to saturate the
            # detector. Reject any truly silly flux levels.
            wherenan = np.where( flux == np.nan )
            if len(wherenan[0]) > 0:
                strg = "***Input flux array contains "
                strg += "%d NaN values." % len(wherenan[0])
                strg += " These will be replaced by 0.0."
                self.logger.warn(strg)
                flux[ wherenan ] = 0.0 
            
            if flux.min() < 0.0:
                whereneg = np.where( flux < 0.0 )
                strg = "***Input flux array contains "
                strg += "%d negative values." % len(whereneg[0])
                self.logger.warn(strg)

            # Calculate an approximate flux that would saturate the detectors
            # after only one frame time.
            if frame_time is None:
                ftime = self.detector.frame_time(
                                        self.samplesum,
                                        self.sampleskip,
                                        refpixsampleskip=self.refpixsampleskip,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode)
            else:
                ftime = frame_time
            if ftime <= 0.0:
                ftime = 1.0
            saturation = self._sca['WELL_DEPTH'] / float(ftime)
            if flux.min() > saturation:
                # Junk data detected. All the pixels will saturate.
                strg = "Detector illumination flux is too high "
                strg += "(%g to %g photons/s/pixel)." % (flux.min(), flux.max())
                strg += "\n   It will saturate ALL the detector pixels. "
                strg += "It should be well below the first frame "
                strg += "saturation level of %g photons/s/pixel." % \
                    float(saturation)
                raise ValueError(strg)
            if flux.max() > (100.0 * saturation):
                # Seriously saturated data detected. It needs to be
                # truncated to prevent an overflow.
                whereover = np.where(flux > 100.0 * saturation)
                strg = "Detector illumination flux is too high "
                strg += "(up to %g photons/s/pixel)." % flux.max()
                strg += "\n   Values that would cause an overflow "
                strg += "will be truncated at %g photons/s/pixel." % \
                    float(saturation * 100.0)
                self.logger.error(strg)
                flux[whereover] = saturation * 100.0
            if flux.max() > saturation:
                wheresat = np.where(flux > saturation)
                strg = "***Input flux array contains values that "
                strg += "could saturate at least %d pixels." % \
                    len(wheresat[0])
                strg += "\n   Maximum flux of %g photons/s/pixel " % flux.max()
                strg += "is greater than first frame "
                strg += "saturation level of %g photons/s/pixel." % \
                    float(saturation)
                self.logger.warn(strg)
            
            if self.fringe_map_data is None:
                self.flux = flux
            else:
                self.flux = flux * self.fringe_map_data
            del flux
        
        if self._verbose > 2:
            strg = "Flux data obtained from input file: shape=" + \
                str(self.flux.shape) + \
                " min=" + str(self.flux.min()) + " max=" + str(self.flux.max())
            if self.fringe_map_data is not None:
                strg +=" and multiplied by fringe map: shape=" + \
                str(self.fringe_map_data.shape) + \
                " min=" + str(self.fringe_map_data.min()) + \
                " max=" + str(self.fringe_map_data.max())
            self.logger.info(strg)
        
        # >>>
        # >>> Reset the detector.
        # >>> Insert an extra FRAME_RESETS resets if needed.
        # >>>
        extra_resets = self._sca['FRAME_RESETS']
        if self._verbose > 4:
            self.logger.debug( "Applying %d extra frame resets." % extra_resets )
        nresets = 1 + extra_resets
        if intnum == 0:
            self.detector.reset(nresets=nresets, new_exposure=True)
        else:
            self.detector.reset(nresets=nresets)            
                
        # >>>
        # >>> Step through the groups.        
        # >>>
        if self._verbose > 1:
            if self.ngroups > 1:
                self.logger.info( "Simulating %d groups for integration %d." % \
                                  (self.ngroups, intnum+1) )
            else:
                self.logger.info( "Simulating %d group for integration %d." % \
                                  (self.ngroups, intnum+1) )
        for group in range(0, self.ngroups):
            if self._verbose > 2:
                self.logger.info( "Integration %d group %d:" % (intnum+1, group+1) )
                
            # The group time is the frame time x number of frames per group.
            if frame_time is None:
                ftime = self.detector.frame_time(
                                        self.samplesum,
                                        self.sampleskip,
                                        refpixsampleskip=self.refpixsampleskip,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode)
            else:
                ftime = frame_time
            time = self.nframes * ftime
            # For all but the first group, the integration includes
            # the group gap.
            if group > 0:
                time += self.groupgap * ftime
            
            # Generate some cosmic ray events during this group of frames
            # interval.
            # TODO: The next few lines could be taken outside the group loop.
            rows = self.shape[0]
            columns = self.shape[1]
            pixsize = self._sca['PIXEL_SIZE']
            
            cosmic_ray_list = \
                self.cosmic_ray_env.generate_events(rows, columns,
                                                    time, pixsize)
            # Hit the detector with the cosmic rays.
            self.detector.hit_by_cosmic_rays(cosmic_ray_list, self.nframes)
            del cosmic_ray_list    
            
            # >>>
            # >>> Integrate on the flux, read the detector and update the
            # >>> exposure data.
            # >>>
            self.detector.integrate(self.flux, time, intnum=intnum)
            total_samples = self.nframes * self.samplesum
            if total_samples < 1:
                total_samples = 1
            # Assist the garbage collector by discarding the previous
            # integration data.
            if self.integration_data is not None:
                del self.integration_data
            self.integration_data = \
                self.detector.readout(subarray=self.subarray,
                                      total_samples=total_samples)
            if self._makeplot and self._verbose > 7:
                mplt.plot_image2D(self.integration_data,
                    xlabel='Columns', ylabel='Rows', withbar=True,
                    title='Readout for integration %d, group %d' % \
                        (intnum,group))
            
            self.exposure_data.set_group(self.integration_data, group, intnum)

        # Return the last set of integration data.
        return self.integration_data

    def exposure(self, nints=None, frame_time=None, start_time=None):
        """
        
        Simulate a detector exposure consisting of one or more
        integrations, each of which consists of one of more groups of
        readouts.
                
        :Parameters:
        
        nints: int, optional
            The number of integrations in this exposure. It must be at
            least 1. If None, the default set by the readout mode is
            used. The total exposure time is nints x integration time.
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter overrides the frame time
            obtained from the current detector readout mode and
            subarray. This can be used, for example, to simulate
            the timing of full frame exposures with subarray-sized
            data (to save memory). If None, the frame time obtained
            from the current readout mode and subarray is used.
        start_time: string or float, optional
            The required clock time of the start of the exposure (MJD days).
            Strings other than 'NOW' are converted to floating point.
            If set to 'NOW', the current date-time is obtained from
            the system clock. By default, the internally stored clock
            time is used.

        :Returns:
        
        integration_data: array_like uint32
            An array of data from the final integration.
        
        :Prerequisites:
        
        Can only be used after the detector illumination data have been
        imported by the read_data method.
        The detector readout mode and integration time should have been
        defined with the set_readout_mode method.
        
        :Raises:
        
        AttributeError
            Raised if the read_data method has not been executed first.
        ValueError
            Raised if any of the parameters are out of range.
        
        """
        # An exposure is not possible until a DetectorArray object
        # and illumination map have been successfully defined.
        if self.illumination_map is None:
            strg = "IlluminationMap object not defined - use read_data first."
            raise AttributeError(strg)
        if self.detector is None:
            strg = "DetectorArray object not defined - use read_data first."
            raise AttributeError(strg)

        # The number of integrations must be at least 1.
        if nints is not None:
            if nints >= 1:
                self.nints = int(nints)
            else:
                self.nints = 1
                
        # Frame time must be positive.
        if frame_time is not None:
            assert frame_time > 0.0

        # Define the number of groups from the readout mode and integration
        # time provided, unless specified explicitly. There must be at least
        # one group.
        if not self.inttime_calculated:
            if self.time_mode == 'time_to_groups':
                self.ngroups = self._time_to_groups(
                                        self.inttime,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode,
                                        frame_time=frame_time)
            else:
                self.inttime = self._groups_to_time(
                                        self.ngroups,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode,
                                        frame_time=frame_time)
                self.inttime_req = self.inttime
            self.inttime_calculated = True

#         # Forward the readout mode to the detector.
#         self.detector.set_readout_mode(self.samplesum, self.sampleskip,
#                                        self.refpixsampleskip, self.nframes)

        # If needed, create a new exposure data object.
        self._new_exposure_data()

        # >>>
        # >>> Execute each integration in turn.
        # >>>
        if self._verbose > 1:
            if self.nints > 1:
                self.logger.info( "Simulating %d integrations." % self.nints )
            else:
                self.logger.info( "Simulating %d integration." % self.nints )
        for intnum in range(0, self.nints):
            self.integration(intnum, frame_time=frame_time)
            
        # If the DARK calibration is not averaged, it is added here.
        if self.detector.simulate_dark_current and not self.detector.dark_averaged:
            if self.detector.dark_map is not None:
                self.logger.info("Adding the DARK calibration from %s" % \
                                 self.detector.dark_map_filename)
                try:
                    self.exposure_data.add_dark( self.detector.dark_map,
                                                 clipvalue=None )
                except ValueError, e:
                    # Catch an exception and instead issue a log message.
                    strg = str(e)
                    strg += ": DARK addition skipped."
                    self.logger.error(strg)
                    
        # If the nonlinearity correction is done by translation table
        # it is applied to the exposure data here.
        if self.detector.simulate_nonlinearity and NONLINEARITY_BY_TABLE:
            rcolumns = self.exposure_data.data.shape[3]
            rcolmiddle = rcolumns//2
            if self.detector.linearity_table_left is not None:
                self.logger.info("Correcting nonlinearity from %s" % \
                                 self.detector.linearity_filename)
                self.exposure_data.apply_translation( \
                    self.detector.linearity_table_left,
                    fromcolumn=0, tocolumn=rcolmiddle )
            if self.detector.linearity_table_right is not None:
                self.exposure_data.apply_translation( \
                    self.detector.linearity_table_right,
                    fromcolumn=rcolmiddle, tocolumn=rcolumns )        
        
        # Copy the primary metadata from the illumination map to the
        # exposure data and append some additional information.
        metadata = Metadata("Metadata of simulated exposure")
        metadata.from_data_object( self.illumination_map )
        
        # Import the primary metadata (but ignore the TYPE keyword)
        for key in metadata.keys():
            if key == 'TYPE':
                continue # Output data will have a different TYPE
            self.metadata[key] = metadata[key]
        comments = metadata.get_comments()
        for comment in comments:
            self.metadata.add_comment(comment)
        histories = metadata.get_history()
        for history in histories:
            self.metadata.add_history(history)

        # Also import the INTENSITY HDU metadata (which contains
        # World Coordinates keywords).
        intensity_metadata = Metadata('Intensity metadata')
        intensity_metadata.from_data_object(self.illumination_map,
                                            hdu_name='INTENSITY')
         
        # Set up the primary metadata according to MIRI-UM-00004-RAL
        simname = "MIRI SCA simulator (MiriTE %s)" % __version__
        strg = "Processed by %s" % simname
        if hasattr(self.metadata, 'add_history'):
            self.metadata.add_history(strg)

        # Housekeeping metadata
        self.metadata["ORIGIN"] = "MIRI European Consortium"
        self.metadata["TELESCOP"] = "JWST"
        self.metadata["INSTRUME"] = "MIRI"
        self.metadata["CREATOR"] = simname

        # Cosmic ray mode
        wordz = self.cosmic_ray_mode.split("+")
        if len(wordz) > 1:
            self.metadata["CRMODE"] = wordz[0]
            self.metadata["CRVAR"] = wordz[1]
        else:
            self.metadata["CRMODE"] = self.cosmic_ray_mode
        if self.cosmic_ray_mode != 'NONE':
            self.metadata = self.cosmic_ray_env.set_metadata(self.metadata)
        self.metadata['CRNUM'] = self.detector.cosmic_ray_count
        self.metadata['CRNUMPIX'] = self.detector.cosmic_ray_pixel_count

        # Detector properties and readout mode
        self.metadata["READPATT"] = self.readout_mode
        self.metadata["NSAMPLES"] = self.samplesum
        self.metadata["SMPSKIP"] = self.sampleskip
        self.metadata["DETROWS"] = self.shape[0]
        self.metadata["DETCOLS"] = self.shape[1]
        #self.metadata["ZROFRAME"] = False

        # Subarray mode
        self.metadata["SUBARRAY"] = self.subarray_str
        (subrows, subcols) = self.detector.get_subarray_shape(self.subarray)
        if self.subarray is not None:
            #(namps, nref) = self._amplifier_count()
            self.metadata["SUBSTRT1"] = self.subarray[1]
            self.metadata["SUBSTRT2"] = self.subarray[0]
            self.metadata["SUBSIZE1"] = subcols
            self.metadata["SUBSIZE2"] = self.subarray[2]
        else:
            self.metadata["SUBSTRT1"] = 1 # 1 indexed
            self.metadata["SUBSTRT2"] = 1 # 1 indexed
            self.metadata["SUBSIZE1"] = subcols
            self.metadata["SUBSIZE2"] = self.shape[0]
        self.metadata["SUBBURST"] = self.subarray_burst_mode

        # Adjust the World Coordinates metadata contained in the INTENSITY
        # metadata, shifting the reference point to account for the
        # reference columns and subarray mode.
        # The first two axes of the exposure data are pixels, and the
        # rest are groups and integrations
        naxes = 2                        
        if self.exposure_data is not None:
            naxes = len(self.exposure_data.shape)
        intensity_metadata['WCSAXES'] = naxes 
                   
        if 'CRPIX1' in self.metadata:
            crpix1 = intensity_metadata['CRPIX1']
        else:
            crpix1 = 1
        if self.subarray_input_str == 'FULL':
            substrt1 = self.metadata["SUBSTRT1"]
        else:
            substrt1 = 1
        if substrt1 > 1:
            # Reference columns not included in subarrays that do not
            # touch the left hand edge.
            intensity_metadata['CRPIX1'] = \
                        crpix1 + 1 - substrt1
        else:
            intensity_metadata['CRPIX1'] = \
                        crpix1 + 1 - substrt1 + self.detector.left_columns
                        
        if 'CRPIX2' in self.metadata:
            crpix2 = intensity_metadata['CRPIX2']
        else:
            crpix2 = 1
        if self.subarray_input_str == 'FULL':
            substrt2 = self.metadata["SUBSTRT2"]
        else:
            substrt2 = 1
        intensity_metadata['CRPIX2'] = crpix2 + 1 - substrt2

        # Add the extra WCS metadata which describes ramp data
        if naxes > 2:
            intensity_metadata['CTYPE3'] = ''
            intensity_metadata['CUNIT3'] = 'groups'
            intensity_metadata['CRVAL3'] = 0.0
            intensity_metadata['CRPIX3'] = 0.0
            intensity_metadata['CDELT3'] = 1.0
        if naxes > 3:
            intensity_metadata['CTYPE4'] = ''
            intensity_metadata['CUNIT4'] = 'integrations'
            intensity_metadata['CRVAL4'] = 0.0
            intensity_metadata['CRPIX4'] = 0.0
            intensity_metadata['CDELT4'] = 1.0

        # Copy the detector metadata to the primary metadata.         
        # NOTE: This function call should only happen AFTER
        # extracting a subarray.
        self.metadata = self.detector.set_metadata( self.metadata )
        
        # Define detailed exposure parameters.
        clktime =  self.detector.clock_time()
        clktime *= 1.0E6   # Convert from seconds to microseconds
        if frame_time is None:
            ftime = self.detector.frame_time(
                                    self.samplesum,
                                    self.sampleskip,
                                    refpixsampleskip=self.refpixsampleskip,
                                    subarray=self.subarray,
                                    burst_mode=self.subarray_burst_mode)
        else:
            ftime = frame_time
        gtime = ftime * self.nframes
        self.metadata["TSAMPLE"] = clktime
        self.metadata["TFRAME"] = ftime
        self.metadata["TGROUP"] = gtime

        self.metadata["ROWRSETS"] = self.refpixsampleskip
        self.metadata["FRMRSETS"] = self._sca['FRAME_RESETS'] # Old keyword
        self.metadata["NRESETS"]  = self._sca['FRAME_RESETS']  # New keyword
                
        # Add integration and exposure metadata keywords.
        self.metadata["EFFINTTM"] = self.inttime
        exptime = self.inttime * self.nints
        reqtime = self.inttime_req * self.nints
        self.metadata["EXPTIME"]  = exptime  # Old keyword
        self.metadata["EFFEXPTM"] = exptime  # New keyword
        self.metadata["REQTIME"] = reqtime
        self.metadata["DETTEMP"] = self.temperature
        
        # Estimate the duration of the exposure in elapsed time.
        # Account for the extra frame resets in between integrations
        # and add an estimated readout time in between exposures.
        # TODO: Improve the readout time overhead estimate.
        extra_time = (int(self._sca['FRAME_RESETS']) * ftime * (self.nints-1))
        extra_time += 1.0e-6 * self.metadata["SUBSIZE1"] * self.metadata["SUBSIZE2"]
        duration = exptime + extra_time
        self.metadata['DURATION'] = duration
        
        # Add the remaining exposure data keywords
        self.metadata["NINTS"] = self.exposure_data.nints
        self.metadata["INTAVG"] = self.exposure_data.intavg
        self.metadata["NGROUPS"] = self.exposure_data.ngroups
        self.metadata["NFRAMES"] = self.exposure_data.nframes
        self.metadata["GROUPGAP"] = self.exposure_data.groupgap
        self.metadata["GRPAVG"] = self.exposure_data.grpavg
        self.metadata['READPATT'] = self.exposure_data.readpatt
            
        # Add a description of the source of the QE measurement
        if hasattr(self.metadata, 'add_comment'):
            if self.qe is not None:
                comment = "QE: From " + os.path.basename(self._sca['QE_FILE'])
                self.metadata.add_comment(comment)
                qecm = self.qe.get_comments()
                self.metadata.add_comment(qecm)
            else:
                self.metadata.add_comment("QE: No QE simulation.")

        # In verbose mode, display the metadata.
        if self._verbose > 3:
            self.logger.debug( str(self.metadata) )

        # Copy primary metadata to the exposure data FITS header.
        for keyw in self.metadata.keys():
            value = self.metadata[keyw]
            try:
                self.exposure_data.set_fits_keyword(keyw, value,
                                                    hdu_name='PRIMARY')
            except KeyError as e:
                strg = "\nKeyword %s cannot be set to %s." % (keyw, str(value))
                strg += " " + str(e)
                strg += " Ignored."
                self.logger.warn(strg)
            except ValueError as e:
                strg = "\nKeyword %s cannot be set to %s." % (keyw, str(value))
                strg += " " + str(e)
                strg += " Ignored."
                self.logger.warn(strg)
             
        # Copy comment and history records to the exposure data FITS header.
        comments = self.metadata.get_comments()
        if isinstance(comments, (tuple,list)):
            for comment in comments:
                self.exposure_data.add_comment(str(comment))
        else:
            self.exposure_data.add_comment(str(comments))
        histories = self.metadata.get_history()
        if isinstance(histories, (tuple,list)):
            for history in histories:
                self.exposure_data.add_history(str(history))
        else:
            self.exposure_data.add_history(str(histories))

        # Copy the World Coordinates metadata from the INTENSITY metadata
        # to the exposure data FITS header.
        # NOTE: Only the first 2 axes of the illumination model WCS
        # parameters are relevant.
        if 'WCSAXES' in intensity_metadata.keys():
            wcsaxes = intensity_metadata['WCSAXES']
        else:
            wcsaxes = 2
            
        crpix = []
        crval = []
        ctype = []
        cunit = []
        cdelt = []
        if wcsaxes > 0:
            # Only the first 2 axes of the PC array are relevant
            pc = np.zeros( [wcsaxes, 2] )
            pc_defined = False
        else:
            pc = None
            pc_defined = True
        
        for axis in range(0, wcsaxes):
            axis1 = axis + 1
            crpix_kw = 'CRPIX%d' % axis1
            crval_kw = 'CRVAL%d' % axis1
            ctype_kw = 'CTYPE%d' % axis1
            cunit_kw = 'CUNIT%d' % axis1
            cdelt_kw = 'CDELT%d' % axis1
            if crpix_kw in intensity_metadata.keys():
                crpix.append( intensity_metadata[crpix_kw] )
            else:
                crpix.append( 0.0 )
            if crval_kw in intensity_metadata.keys():
                crval.append( intensity_metadata[crval_kw] )
            else:
                crval.append( 0.0 )
            if ctype_kw in intensity_metadata.keys():
                ctype.append( intensity_metadata[ctype_kw] )
            else:
                ctype.append( '' )
            if cunit_kw in intensity_metadata.keys():
                cunit.append( intensity_metadata[cunit_kw] )
            else:
                cunit.append( '' )
            if cdelt_kw in intensity_metadata.keys():
                cdelt.append( intensity_metadata[cdelt_kw] )
            else:
                cdelt.append( 1.0 )
                
            for other in range(0, 2):
                other1 = other + 1
                pc_kw = 'PC%d_%d' % (axis1, other1)
                if pc_kw in intensity_metadata.keys():
                    pc_defined = True
                    pc[axis,other] = intensity_metadata[pc_kw]

        # If no PC values have been copied from the INTENSITY metadata
        # leave the PC array blank.
        if not pc_defined:
            pc = None

        # Additional, WCS-related information.
        if 'S_REGION' in intensity_metadata.keys():
            s_region = intensity_metadata['S_REGION']
        else:
            s_region = None
        if 'WAVSTART' in intensity_metadata.keys():
            waverange_start = intensity_metadata['WAVSTART']
        else:
            waverange_start = None
        if 'WAVEND' in intensity_metadata.keys():
            waverange_end = intensity_metadata['WAVEND']
        else:
            waverange_end = None
        if 'SPORDER' in intensity_metadata.keys():
            spectral_order = intensity_metadata['SPORDER']
        else:
            spectral_order = None
        if 'V2_REF' in intensity_metadata.keys():
            v2_ref = intensity_metadata['V2_REF']
        else:
            v2_ref = None
        if 'V3_REF' in intensity_metadata.keys():
            v3_ref = intensity_metadata['V3_REF']
        else:
            v3_ref = None
        if 'VPARITY' in intensity_metadata.keys():
            vparity = intensity_metadata['VPARITY']
        else:
            vparity = None
        if 'V3I_YANG' in intensity_metadata.keys():
            v3yangle = intensity_metadata['V3I_YANG']
        else:
            v3yangle = None
        if 'RA_REF' in intensity_metadata.keys():
            ra_ref = intensity_metadata['RA_REF']
        else:
            ra_ref = None
        if 'S_REGION' in intensity_metadata.keys():
            s_region = intensity_metadata['S_REGION']
        else:
            s_region = None
        if 'DEC_REF' in intensity_metadata.keys():
            dec_ref = intensity_metadata['DEC_REF']
        else:
            dec_ref = None
        if 'ROLL_REF' in intensity_metadata.keys():
            roll_ref = intensity_metadata['ROLL_REF']
        else:
            roll_ref = None
                    
        self.exposure_data.set_wcs_metadata(wcsaxes=wcsaxes, crpix=crpix,
                crval=crval, ctype=ctype, cunit=cunit, cdelt=cdelt, pc=pc,
                s_region=s_region, waverange_start=waverange_start,
                waverange_end=waverange_end, spectral_order=spectral_order,
                v2_ref=v2_ref, v3_ref=v3_ref, vparity=vparity,
                v3yangle=v3yangle, ra_ref=ra_ref, dec_ref=dec_ref,
                roll_ref=roll_ref)

        # Set the exposure start and end times and update the stored clock time.
        if start_time is None:
            start_time = _clock_to_mjd( self.clock_time )
        self.exposure_data.set_exposure_times( start_time=start_time)

        (exposure_time, duration, start_time, mid_time, end_time) = \
            self.exposure_data.get_exposure_times()
        self.clock_time = _mjd_to_clock( end_time )
        
        # In verbose mode, display statistics for the last integration
        # and group.
        if self._verbose > 3:
            self.logger.debug( self.exposure_data.statistics() )
 
        # Return the data from the last integration.
        return self.integration_data

    def _new_exposure_data(self):
        """
        
        Helper function to set up a exposure data object for a new
        exposure or integration.
        
        """     
        # Force a new exposure data object for each exposure.
        if self.exposure_data is not None:
            del self.exposure_data
            self.exposure_data = None
        # If needed, create the object to hold the exposure data.
        # NOTE: The output subarray mode is assumed to be fixed once
        # the exposure data object has been created. Changing the
        # output subarray mode will invalidate the exposure data
        # object.
        if self.exposure_data is None:
            # Get exposure data size and reference output size according
            # to the subarray mode.
            data_shape = self.detector.get_subarray_shape(self.subarray)
            refout_shape = self.detector.refout_shape

            if self._verbose > 1:
                self.logger.info( "Creating exposure_data with " + \
                    str(data_shape[0]) + " rows x " + \
                    str(data_shape[1]) + " columns plus " + \
                    str(self.ngroups) + " groups and " + \
                    str(self.nints) + " ints." )

            # The exposure data model used depends on the choice of
            # output file format.
            if self.fileformat.upper() == 'STSCI':
                # Convert the bad pixel mask into a PIXELDQ array.
                if self.detector.bad_pixels is not None and self.include_pixeldq:
                    new_dq = np.zeros_like( self.detector.bad_pixels )
                    # TODO: This is a fixed mapping which might change.
                    # See "https://confluence.stsci.edu/display/JWSTPWG/JWST+Calibration+Reference+Files%3A+File+Formats+for+the+Build+6+Pipeline"
                    # DO_NOT_USE, DEAD, HOT, UNRELIABLE_SLOPE, RC, NON_SCIENCE
                    oldflags = (1,    2,    4,        8,    16, 512)
                    newflags = (1, 1024, 2048, 16777216, 16384, 512)
                    for oldfl, newfl in zip(oldflags, newflags):
                        testmask = self.detector.bad_pixels & oldfl
                        where = np.where(testmask != 0)
                        new_dq[where] = new_dq[where] | newfl           

                    # If needed, extract a subarray from the bad pixel mask.
                    # Note that subarray rows and columns start at 1 but
                    # numpy array indices start at 0.
                    if self.subarray is not None:
                        r1 = int(self.subarray[0]) - 1
                        c1 = int(self.subarray[1]) - 1
                        r2 = r1 + data_shape[0]
                        c2 = c1 + data_shape[1]
                        window_dq = new_dq[r1:r2, c1:c2]
                    else:
                        window_dq = new_dq
                else:
                    # Either there is no bad pixel mask or the PIXELDQ array
                    # is not to be saved.
                    new_dq = None
                    window_dq = None 
 
                self.exposure_data = MiriExposureModel(
                                        data_shape[0],
                                        data_shape[1],
                                        self.ngroups, self.nints,
                                        self.readout_mode,
                                        refout_shape[0],
                                        refout_shape[1],
                                        pixeldq=window_dq,
                                        grpavg=self.avggrps,
                                        intavg=self.avgints,
                                        nframes=self.nframes,
                                        groupgap=self.groupgap,
                                        include_pixeldq=self.include_pixeldq,
                                        include_groupdq=self.include_groupdq)
                # Initialise the metadata from the illumination model,
                # but do not copy the data units or data type information.
                if self.illumination_map is not None:
                    self.exposure_data.copy_metadata(self.illumination_map,
                                                     ignore=['units', 'type'])
                # Define the correct data units
                if self.simulate_gain:
                    self.exposure_data.meta.data.units = 'DN'
                    self.exposure_data.meta.refout.units = 'DN'
                else:
                    self.exposure_data.meta.data.units = 'electrons'
                    self.exposure_data.meta.refout.units = 'electrons'
# UNNECESSARY. ALREADY DEFINED IN CONSTRUCTOR AND COPIED TO METADATA IN EXPOSURE()
#                 # Ensure the exposure parameters are included in the metadata
#                 self.exposure_data.meta.exposure.nints = self.nints
#                 self.exposure_data.meta.exposure.ngroups = self.ngroups
#                 self.exposure_data.meta.exposure.nframes = self.nframes
#                 self.exposure_data.meta.exposure.groupgap = self.groupgap
#                 self.exposure_data.meta.exposure.groups_averaged = self.avggrps
#                 self.exposure_data.meta.exposure.integrations_averaged = self.avgints
                
                # Define the exposure type (if not already given)
                if not self.exposure_data.meta.exposure.type:
                    self.exposure_data.set_exposure_type(self.detectorid,
                            filter=None, subarray=self.subarray_str)
            else:
                self.exposure_data = ExposureData(
                                        data_shape[0], data_shape[1],
                                        self.ngroups, self.nints,
                                        self.readout_mode,
                                        grpavg=self.avggrps,
                                        intavg=self.avgints,
                                        nframes=self.nframes,
                                        groupgap=self.groupgap)
                # Define the correct data units
                if self.simulate_gain:
                    self.exposure_data.set_fits_keyword('BUNIT', 'DN')
                else:
                    self.exposure_data.set_fits_keyword('BUNIT', 'electrons')

    def clear_exposure_data(self):
        """
        
        Explicitly clear the exposure data contained within the object.
        Use this function to reduce memory usage when the data object
        returned by simulate_pipe() is no longer required.
        
        :Parameters:
        
        None
        
        :Returns:
        
        None
        
        """     
        # Delete the exposure data if it exists.
        if self.exposure_data is not None:
            del self.exposure_data
            self.exposure_data = None
    
    def _groups_to_time(self, ngroups, subarray=None, burst_mode=True,
                        frame_time=None):
        """
        
        Convert a number of groups to an integration time.
        
        """
        # Define the integration time from the frame time, frames per group
        # and number of groups provided. There must be at least one group.
        # If there are 2 or more groups and a non-zero group gap, the
        # time taken for the dropped frames must be added to the
        # integration time.
        if frame_time is None:
            ftime = self.detector.frame_time(
                                    self.samplesum,
                                    self.sampleskip,
                                    refpixsampleskip=self.refpixsampleskip,
                                    subarray=self.subarray,
                                    burst_mode=self.subarray_burst_mode)
        else:
            ftime = frame_time
        if ngroups > 1:
            time = ftime * int(ngroups) * self.nframes
            if self.groupgap > 0:
                time += ftime * self.groupgap * (int(ngroups)-1)
        else:
            time = ftime * self.nframes

        if self._verbose > 2:
            if subarray is None:
                substr = 'FULL'
            else:
                substr = subarray
            if burst_mode:
                burststr = " (burst mode)"
            else:
                burststr = ""
            strg = "samplesum=%d, sampleskip=%d and subarray=%s%s gives frame time=%.3fs; " % \
                (self.samplesum, self.sampleskip, substr, burststr, ftime)
            strg += "ngroups=%d nframes=%d groupgap=%d " % \
                (ngroups, self.nframes, self.groupgap)
            strg += "gives an integration time of %.2fs." % time
            self.logger.info( strg )
                
        return time

    def _time_to_groups(self, time, subarray=None, burst_mode=False,
                        frame_time=None):
        """
        
        Convert an integration time to number of groups.
        
        """
        # Preserve the requested integration time.
        self.inttime_req = float(time)
        # Define the number of groups from the readout mode and integration
        # time provided, rounding the number of groups up so the integration
        # time is at least the amount specified. If avggrps is not 1 the
        # number is also rounded up to the nearest whole multiple of avggrps.
        # There must be at least one group.
        if frame_time is None:
            ftime = self.detector.frame_time(
                                    self.samplesum,
                                    self.sampleskip,
                                    refpixsampleskip=self.refpixsampleskip,
                                    subarray=self.subarray,
                                    burst_mode=self.subarray_burst_mode)
        else:
            ftime = frame_time
        gtime = ftime * self.nframes
        dtime = ftime * self.groupgap
        ngroups = int(np.ceil((self.inttime_req+dtime) / (gtime+dtime)))
        if ngroups < 1:
            ngroups = 1
        if self.avggrps > 1:
            ngroups = int(np.ceil(ngroups/self.avggrps) * self.avggrps)
        if ngroups < 2:
            self.groupgap = 0
            dtime = 0.0
            
        # Record the actual integration time.
        self.inttime = (ngroups * gtime) + ((ngroups-1) * dtime)
            
        if self._verbose > 2:
            if subarray is None:
                substr = 'FULL'
            else:
                substr = subarray
            if burst_mode:
                burststr = " (burst mode)"
            else:
                burststr = ""
            strg = "samplesum=%d, sampleskip=%d and subarray=%s%s gives frame time=%.3fs; " % \
                (self.samplesum, self.sampleskip, substr, burststr, ftime)
            strg +=" Integration time of %.2fs with nframes=%d groupgap=%d " % \
                (time, self.nframes, self.groupgap)
            if ngroups > 1:
                strg += "converts to %d groups." % ngroups
            else:
                strg += "converts to %d group." % ngroups
            self.logger.info( strg )
                
        return ngroups

    def get_exposure_times(self):
        """
        
        Return the exposure time metadata to the caller
        
        :Returns:
        
        (exposure_time, duration, start_time, mid_time, end_time): tuple of 5 floats
            The exposure time in seconds, clock duration in seconds, exposure
            start time, mid time and end time in MJD days. Undefined values are
            returned None.
       
        """
        # There must be some exposure data containing the end time.
        if self.exposure_data is None:
            strg = "Exposure data not defined - " + \
                "use integration or exposure first."
            raise AttributeError(strg)
        
        return self.exposure_data.get_exposure_times()
          
    def write_data(self, filename, datashape='hypercube', overwrite=False):
        """
        
        Write the simulated data to a given output FITS file.
        
        :Parameters:
        
        filename: string
            The name of the file to be created.
        datashape: string, optional, default='hypercube'
            The SCI data shape to be written.
            
            * 'hypercube' - write the SCI data to a 4 dimensional FITS
              image with separate columns x rows x groups x integrations
              dimensions.
            * 'cube' - append the groups and integrations to make a
              3-dimensional FITS image with columns x rows x (groups and
              integrations) dimensions.
            * 'slope' - combine the groups and integrations with a crude
              straight line fit to make slope data, contained in a
              2-dimensional FITS image with columns x rows dimensions.

        overwrite: bool, optional, default=False
            Parameter passed to pyfits.HDUlist.writeto
            
        :Prerequisites:
        
        Can only be used after the exposure data have been generated by
        the integration or exposure methods.
        
        :Raises:
        
        AttributeError
            Raised if no exposure data exists.
        
        """
        # There must be some exposure data to write out.
        if self.exposure_data is None:
            strg = "Exposure data not defined - " + \
                "use integration or exposure first."
            raise AttributeError(strg)

        # Write the exposure data to the given output file with the given
        # file format and shape.
        if self.fileformat.upper() == 'STSCI':
            if datashape == 'cube' or datashape == 'CUBE':
                # The exposure data must be crunched down into a cube.
                cube_data = self.exposure_data.cube_data()
                if self.include_pixeldq:
                    cube_model = MiriMeasuredModel(data=cube_data,
                                            dq=self.exposure_data.pixeldq)
                else:
                    cube_model = MiriMeasuredModel(data=cube_data)
                # Copy the metadata.
                cube_model.copy_metadata(self.exposure_data)
                # The new exposure data has only 3 WCS aces
                cube_model.meta.wcsinfo.wcsaxes = 3
                # Write the cube data to the given output file.
                cube_model.save(filename, overwrite=overwrite)
                # Tidy up the temporary cube model.
                del cube_data, cube_model
            elif datashape == 'slope' or datashape == 'SLOPE':
                # The exposure data must be straight-line fitted to make slope data.
                slope_data = self.exposure_data.slope_data(diff_only=True)
                if self.include_pixeldq:
                    slope_model = MiriMeasuredModel(data=slope_data,
                                            dq=self.exposure_data.pixeldq)
                else:
                    slope_model = MiriMeasuredModel(data=slope_data)

                # Copy the metadata.
                # TODO: This does not copy the simulation metadata because it
                # is not present the the MiriMeasuredModel schema.
                slope_model.copy_metadata(self.exposure_data)
                # The new slope data has only 2 WCS aces
                cube_model.meta.wcsinfo.wcsaxes = 2
                 # Write the cube data to the given output file.
                slope_model.save(filename, overwrite=overwrite)
                # Tidy up the temporary cube model.
                del slope_data, slope_model
            else:
                # Leave the data in hypercube format.         
                # Write the exposure data to the given output file.
                self.exposure_data.save(filename, overwrite=overwrite)
        else:
            self.exposure_data.save(filename, fileformat=self.fileformat,
                                    datashape=datashape, overwrite=overwrite)
        
        # Delete the exposure data object so a new one can be created
        # next time an exposure is started.
        del self.exposure_data
        self.exposure_data = None
        
    def __str__(self):
        """
        
        Returns a string describing the SensorChipAssembly object.
        
        """
        strg = "SensorChipAssembly: Detector ID %s, name %s." % \
            (self._sca["SCA_ID"], self._sca["DETECTOR"])
        strg += "\n" + self.readout_mode_str(prefix='  ')
        strg += "\n" + self.temperature_str(prefix='  ')
        strg += "\n" + self.cosmic_ray_str(prefix='  ')
        if self.detector is not None and self._verbose > 2:
            strg += "\n" + self.detector.__str__()
        return strg
    
    def readout_mode_str(self, prefix=''):
        """
        
        Returns a string describing the readout mode.
        
        """
        strg = prefix + "Detector readout mode is " + \
            "%s (samplesum=%d, sampleskip=%d, nframe=%d, groupgap=%s) " % \
            (self.readout_mode, self.samplesum, self.sampleskip,
             self.nframes, self.groupgap)
        strg += "\n%swith %d integrations " % (prefix, self.nints)
        
        if self.time_mode == 'groups_to_time':
            strg += "and ngroups=%d defined explicitly." % self.ngroups
        elif self.time_mode == 'time_to_groups':
            if self.ngroups is None:
                strg += "and ngroups to be defined from the " + \
                    "integration time of %.2f seconds." % self.inttime
            else:
                strg += "ngroups=%d defined from the " % self.ngroups
                strg += "integration time of %.2f seconds." % self.inttime        
        else:
            strg += "and ngroups=%d defined by default." % self.ngroups

        if self.subarray is None:
            strg += "\n%sDetector subarray mode is %s." % \
                (prefix, self.subarray_str)
        else:
            strg += "\n%sDetector subarray mode is %s " % \
                (prefix, self.subarray_str)
            strg += "(%d %d %d %d)." % self.subarray
        return strg

    def simulate_files(self, inputfile, outputfile, detectorid, scale=1.0,
            fringemap=None, readout_mode=None, subarray=None, burst_mode=True,
            frame_time=None, inttime=None, ngroups=None, nints=None,
            start_time=None, wait_time=0.0, temperature=None,
            cosmic_ray_mode=None, fileformat='STScI',
            datashape='hypercube', include_pixeldq=True, include_groupdq=False,
            overwrite=False, qe_adjust=True, simulate_poisson_noise=True,
            simulate_read_noise=True, simulate_ref_pixels=True,
            simulate_bad_pixels=True, simulate_dark_current=True,
            simulate_flat_field=True, simulate_gain=True,
            simulate_nonlinearity=True, simulate_drifts=True,
            simulate_latency=True, cdp_ftp_path=SIM_CDP_FTP_PATH,
            readnoise_version='', bad_pixels_version='',
            flat_field_version='', linearity_version='', gain_version='',
            makeplot=False, seedvalue=None, verbose=2):
        """
    
        Runs the MIRI SCA simulator on the given input file, generating the
        given output file.
                 
        :Parameters:
    
        inputfile: string or list of strings
            If a string, the name of the FITS file containing detector
            illumination data (normally created by another MIRI simulator).
            If a list of strings, a list of files to be processed.
        outputfile: string or list of strings
            If a string, the name of the output file to contain the
            simulated SCA data.
            If a list, the 
        detectorid: string
            Detector ID, identifying a particular detector.
            The MIRI instrument has three detectors: 'MIRIMAGE',
            'MIRIFULONG' and 'MIRIFUSHORT'. (These correspond to the imager
            and the long and short wavelength arms of the MRS respectively.)
        scale: float, optional, default=1.0
            An optional scale factor to apply to the intensity data.
            This is for debugging and testing only - the input data
            should already be scaled to photons/second/pixel.
        fringemap: string, optional
            The name of a file containing a fringe map to be used to
            modulate the input illumination. If not specified, no fringe
            map will be used. 
        readout_mode: string, optional
            Readout mode, which can be one of the following common options:
        
            * 'SLOW' - 10 samples per readout and defaults of ngroups=10
               nints=1.
            * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
              nints=10.
            * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
              averaged.
            * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
              averaged.
          
            The following unusual options are also available for testing:
        
            * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
              4 frames per group and a gap of 8 frames between each group.
            * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
              averaged and default ngroups=1.
            * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
              averaged and default ngroups=4.
            * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
              4 frames per group and a gap of 8 frames between each group.
          
            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file. Failing that, it will
            be obtained from the detector properties default.
        subarray: string, optional
            Detector subarray mode for output. This can be one 'FULL',
            'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
            'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
            full-frame and the other modes read out portions of the detector
            as described in the MIRI Operational Concept Document.
            The simulator can also accept the other subarray modes defined
            in detector_properties.py for test purposes.
            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file (unless the input data
            specifies a non-standard subarray). Failing that, it will be
            obtained from the detector properties default (FULL).
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter overrides the frame time obtained
            from the detector readout mode and subarray. This can be
            used, for example, to simulate the timing of full frame
            exposures with subarray-sized data (to save memory). If None,
            the frame time obtained from the current readout mode and
            subarray is used.
        inttime: float, optional
            The integration time in seconds to be simulated. This parameter
            will be ignored if ngroups is provided as well.
            *NOTE: Any requested integration time will be rounded up to
            the time resulting from the nearest whole number of groups.*
        ngroups: int, optional
            The number of readout groups. If None, this is derived from
            the integration time and readout mode. Otherwise the parameter
            overrides the integration time and sets the number of groups
            directly. There must be at least one group.
        nints: int, optional
            The number of integrations per exposure. It must be at least 1.
            If None, the default set by the readout mode is used.
            The total exposure time is nints x inttime.
        start_time: string or float, optional
            The required clock time of the start of the exposure (MJD days).
            Strings other than 'NOW' are converted to floating point.
            If set to 'NOW', the current date-time is obtained from
            the system clock. By default, the internally stored clock
            time is used.
        wait_time: float, optional
            The wait time in seconds that has elapsed since the previous
            exposure. The default is 0.0 seconds.
        temperature: float, optional, default=detector target temperature
            The temperature of the detector, which will determine the dark
            current and readout noise.
        cosmic_ray_mode: string, optional, default=cosmic ray properties
            The cosmic ray environment mode to be simulated. Available
            modes are:
        
            * 'NONE' - No cosmic rays.
            * 'SOLAR_MIN' - Solar minimum
            * 'SOLAR_MAX' - Solar maximum
            * 'SOLAR_FLARE' - Solar flare (worst case scenario)

            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file. Failing that, it will
            be obtained from the cosmic ray properties default.
        fileformat: string, optional, default='STScI'
            The kind of file format to be written.
            
            * 'STScI' - use the STScI level 1 model for the JWST
              DMS pipeline.
            * 'FITSWriter' - emulate the format written by the FITSWriter
              during MIRI VM, FM and CV tests and read by the DHAS.
              
        datashape: string, optional, default='hypercube'
            The SCI data shape to be written.
            
            * 'hypercube' - write the SCI data to a 4 dimensional FITS image
              with separate columns x rows x groups x integrations
              dimensions.
            * 'cube' - append the groups and integrations to make a 3
              dimensional FITS image with columns x rows x (groups and
              integrations) dimensions.
            * 'slope' - combine the groups and integrations with a crude
              straight line fit to make slope data, contained in a
              2-dimensional FITS image with columns x rows dimensions.
          
        include_pixeldq: boolean, optional, default=True
            A flag that may be used to switch on the inclusion of the
            PIXELDQ data array in the output exposure data.
            By default, this array is included (and contains the bad pixel
            mask used to generate the simulated data).
        include_groupdq: boolean, optional, default=False
            A flag that may be used to switch on the inclusion of the
            GROUPDQ data array in the output exposure data.
            By default, this array is not included.
        overwrite: bool, optional, default=False
            Parameter passed to pyfits.HDUlist.writeto
        qe_adjust: boolean, optional, default=True
            A flag that may be used to switch off quantum efficiency
            adjustment (for example to observe what effects in a simulation
            are caused by QE). *NOTE: When QE adjustment is turned off,
            the input illumination is assumed to be in electrons/second.*
        simulate_poisson_noise: boolean, optional, default=True
            A flag that may be used to switch off Poisson noise (for example
            to observe what effects in a simulation are caused by Poisson
            noise).
        simulate_read_noise: boolean, optional, default=True
            A flag that may be used to switch off read noise (for example
            to observe what effects in a simulation are caused by read
            noise).
        simulate_ref_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of reference
            pixels in the data. 
        simulate_bad_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of bad
            pixels in the data, even if a bad pixel map containing bad pixels
            is specified in the detector properties.
        simulate_dark_current: boolean, optional, default=True
            A flag that may be used to switch off the addition of dark
            current (for example to observe what effects in a simulation are
            caused by dark current).
        simulate_flat_field: boolean, optional, default=True
            A flag that may be used to switch off the simulation of the pixel
            flat-field (for example to observe the effect of quantum
            efficiency in isolation).
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
        cdp_ftp_path: str, optional, default=None
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
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator before
            generating the test data.
            If not specified, a value of None will be sent, which
            randomises the seed.
        verbose: int, optional, default=2
            Verbosity level. Activates print statements when non-zero.
        
            * 0 - no output at all except for error messages
            * 1 - warnings and errors only
            * 2 - normal output
            * 3, 4, 5 - additional commentary
            * 6, 7, 8 - extra debugging information
        
        :Raises:
    
        ValueError
            Raised if any of the simulation parameters are out of range.
            Also raised if the value of a parameter is invalid or
            inappropriate.
        TypeError
            Raised if any of the simulation parameters are of the
            wrong type, size or shape.
        KeyError
            Raised if any of the simulation parameters do not contain
            a recognised keyword.
        IndexError
            Raised if an attempt is made to access beyond the bounds of
            the data.
        IOError
            Raised if a file cannot be read, interpreted or written.
        ImportError
            A delayed ImportError is raised if there is an attempt to use
            an optional library (such as matplotlib) which is not
            available. The software fails immediately if a compulsory
            library is not available.
        AttributeError
            Raised if simulation methods are executed in the wrong order,
            meaning that a necessary attribute is missing. A programming
            error.
        
        """
        if verbose > 1:
            _report_sca_parameters(inputfile, outputfile, detectorid,
                                   readout_mode, fringemap,
                                   fileformat, datashape, qe_adjust,
                                   simulate_poisson_noise, simulate_read_noise,
                                   simulate_ref_pixels, simulate_bad_pixels,
                                   simulate_dark_current, simulate_flat_field,
                                   simulate_gain, simulate_nonlinearity,
                                   simulate_drifts, simulate_latency,
                                   cdp_ftp_path,
                                   readnoise_version, bad_pixels_version,
                                   flat_field_version, linearity_version,
                                   gain_version,
                                   logger=self.logger)

        # Set the seed for the np.random function.
        np.random.seed(seedvalue)

        # Extract the properties of the particular sensor chip assembly
        # being simulated.
        try:
            fpm = detector_properties.get('DETECTORS_DICT', str(detectorid))
        except KeyError:
            strg = "%s is not a known detector ID." % str(detectorid)
            raise KeyError(strg)

        # The simulation may be applied either to a single file or a list of
        # files. If a single file name has been specified, convert it to a list.
        if not isinstance(inputfile, (tuple,list)):
            inputfile = [inputfile]
            outputfile = [outputfile]
        
        # If any of these arguments are explicitly set to None they need
        # to be defaulted. The default values are either extracted from
        # the FITS header of the (first) input file or, if not found, set
        # to a fixed default.
        firstfile = inputfile[0]
        header = get_file_header(firstfile)
        (readout_mode, subarray, frame_time, temperature, cosmic_ray_mode) = \
            _get_parameter_defaults(fpm, header, readout_mode, subarray,
                                    frame_time, temperature, cosmic_ray_mode,
                                    verbose=verbose, logger=self.logger)
    
        # Set up the required detector parameters.
        self.setup(detectorid, readout_mode=readout_mode, subarray=subarray,
              burst_mode=burst_mode, inttime=inttime, ngroups=ngroups,
              nints=nints, temperature=temperature,
              cosmic_ray_mode=cosmic_ray_mode, fileformat=fileformat,
              include_pixeldq=include_pixeldq, include_groupdq=include_groupdq,
              qe_adjust=qe_adjust,
              simulate_poisson_noise=simulate_poisson_noise,
              simulate_read_noise=simulate_read_noise,
              simulate_ref_pixels=simulate_ref_pixels,
              simulate_bad_pixels=simulate_bad_pixels,
              simulate_dark_current=simulate_dark_current,
              simulate_flat_field=simulate_flat_field,
              simulate_gain=simulate_gain,
              simulate_nonlinearity=simulate_nonlinearity,
              simulate_drifts=simulate_drifts,
              simulate_latency=simulate_latency,
              cdp_ftp_path=cdp_ftp_path,
              readnoise_version=readnoise_version,
              bad_pixels_version=bad_pixels_version,
              flat_field_version=flat_field_version, 
              linearity_version=linearity_version, 
              gain_version=gain_version,
              makeplot=makeplot, verbose=verbose)
        for ii in range(0, len(inputfile)):
            infile = inputfile[ii]
            outfile = outputfile[ii]
            if verbose > 1:
                separator = "Exposure %d:" % (ii+1)
                if outfile:
                    self.logger.info( "\n%s %s --> %s" % (separator, infile, outfile) )
                else:
                    self.logger.info( "\n%s %s --> (no output file)" % (separator, infile) )

            # Read the detector illumination from the given FITS file,
            # scaling if it required.
            self.read_data(infile, scale=scale)
    
            # Read the fringe map if a name has been provided.
            if fringemap is not None and fringemap:
                self.read_fringe_map(fringemap, ftype='FITS')

            # Wait for the given elapsed time.
            if wait_time > 0.0:
                self.wait(wait_time)

            # Simulate an exposure and return some simulated data.
            simulated_data = self.exposure(frame_time=frame_time,
                                           start_time=start_time)
    
            # Report the exposure times.
            if verbose > 1:
                (exposure_time, duration, start_time, mid_time, end_time) = \
                    self.get_exposure_times()
                strg = "Exposure time %.2fs (duration %.2fs) " % \
                    (exposure_time, duration)
                if verbose > 2:
                    strg += "started at %.2f and finished at %.2f" % \
                        (start_time, end_time)
                self.logger.info( strg )
    
            # Plot the simulated data if requested.
            if makeplot:
                strg = "Simulated exposure data from %s" % inputfile
                self.plot(description=strg)
                if verbose > 3:
                    # Development and testing - plot a ramp for the central pixel
                    strg2 = "\n%s" % inputfile
                    self.plot_ramp(simulated_data.shape[0]/2,
                                  simulated_data.shape[1]/2, strg2)

            if outfile:
                # Write the results to a data cube or level 1 FITS file.
                if verbose > 1:
                    if fileformat == 'FITSWriter':
                        self.logger.info(
                            "Writing FITSWriter fileformat FITS file: " + \
                            outfile )
                    elif fileformat == 'STSCI' or fileformat == 'STScI':
                        strg = "Writing level 1 FITS file matching STScI ramp data model"
                        if self.include_pixeldq:
                            strg += " +PIXELDQ"
                        if self.include_groupdq:
                            strg += " +GROUPDQ"
                        strg += ": "
                        strg += outfile
                        self.logger.info(strg)
                    else:
                        self.logger.info(
                            "Writing OLD FORMAT level 1 FITS file: " + outfile )
                self.write_data(outfile, datashape=datashape, overwrite=overwrite)

    def simulate_pipe(self, illumination_map, scale=1.0,
                      fringemap=None, readout_mode=None, subarray=None,
                      burst_mode=True, frame_time=None,
                      inttime=None, ngroups=None, nints=None,
                      start_time=None, wait_time=0.0, temperature=None,
                      cosmic_ray_mode=None, include_pixeldq=True,
                      include_groupdq=False, overwrite=False, qe_adjust=True,
                      simulate_poisson_noise=True, simulate_read_noise=True,
                      simulate_ref_pixels=True, simulate_bad_pixels=True,
                      simulate_dark_current=True, simulate_flat_field=True,
                      simulate_gain=True, simulate_nonlinearity=True,
                      simulate_drifts=True, simulate_latency=True,
                      cdp_ftp_path=SIM_CDP_FTP_PATH,
                      readnoise_version='', bad_pixels_version='',
                      flat_field_version='', linearity_version='',
                      gain_version='',
                      makeplot=False, seedvalue=None, verbose=2):
        """
    
        Runs the MIRI SCA simulator on the given MiriIluminationModel
        and returns a MiriExposureModel
        This is an alternative to the file-based method, simulate_files(),
        which can be used by pipeline code which needs to convert one
        data model into another.
        
        NOTE: The clear_exposure_data() method can be used to delete
        the internal reference to the exposure data and save memory in
        between simulations.
        
        :Caveats:
        
        It is possible to create a MiriIlluminationModel object of any
        size and shape and give it to this function, and it is common to
        test the simulation using small datasets. Please note the following
        caveats when giving SCASim arbitrary-sized illumination arrays or
        arbitrary parameter settings:
        
           * If you give SCASim an illumination array of arbitrary size and
             a subarray mode of 'FULL', SCASim will create a detector with
             the same number of pixels as the illumination array. If you
             specify a particular subarray mode, the input array must either
             be full-sized (1024x1024) or be exactly the same size as the
             named subarray.
           * The frame time of the detector readout depends on the number
             of pixels simulated - the smaller the detector the faster the
             frame time. The frame time also depends on whether reference
             pixels are included (by setting the simulate_ref_pixels flag).
             Therefore, to make accurate timing tests or flux tests, you
             need to provide a full-sized array (1024x1024) and use
             simulate_ref_pixels=True. Be aware that the setting of the
             qe_adjust flag will affect the measured flux level.
           * Some simulation steps (bad pixels, dark current, pixel
             flat-field, and gain) use MIRI calibration data products (CDPs).
             Results are unpredictable if the CDP has been defined for a
             different-sized detector than the one being simulated. If in
             doubt, ensure the input array is full-size or sized according
             to a named subarray mode.
                 
        :Parameters:
    
        illumination_map: MiriIlluminationModel object
            A MIRI illumination map object which describes the input
            illumination.
        scale: float, optional, default=1.0
            An optional scale factor to apply to the intensity data.
            This is for debugging and testing only - the input data
            should already be scaled to photons/second/pixel.
        fringemap: string, optional
            The name of a file containing a fringe map to be used to
            modulate the input illumination. If not specified, no fringe
            map will be used. 
        readout_mode: string, optional
            Readout mode, which can be one of the following common options:

            * 'SLOW' - 10 samples per readout and defaults of ngroups=10
              and nints=1.
            * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
              nints=10.
            * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
              averaged.
            * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
              averaged.
          
            The following unusual options are also available for testing:
        
            * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
              4 frames per group and a gap of 8 frames between each group.
            * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
              averaged and default ngroups=1.
            * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
              averaged and default ngroups=4.
            * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
              4 frames per group and a gap of 8 frames between each group.
        
            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file. Failing that, it will
            be obtained from the detector properties default.
        subarray: string, optional
            Detector subarray mode for output. This can be one 'FULL',
            'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
            'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
            full-frame and the other modes read out portions of the detector
            as described in the MIRI Operational Concept Document.
            The simulator can also accept the other subarray modes defined
            in detector_properties.py for test purposes.
            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file (unless the input data
            specifies a non-standard subarray). Failing that, it will be
            obtained from the detector properties default (FULL).
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter overrides the frame time obtained
            from the detector readout mode and subarray. This can be
            used, for example, to simulate the timing of full frame
            exposures with subarray-sized data (to save memory). If None,
            the frame time obtained from the current readout mode and
            subarray is used.
        inttime: float, optional
            The integration time in seconds to be simulated. This parameter
            will be ignored if ngroups is provided as well.
            *NOTE: Any requested integration time will be rounded up to
            the time resulting from the nearest whole number of groups.*
        ngroups: int, optional
            The number of readout groups. If None, this is derived from
            the integration time and readout mode. Otherwise the parameter
            overrides the integration time and sets the number of groups
            directly. There must be at least one group.
        nints: int, optional
            The number of integrations per exposure. It must be at least 1.
            If None, the default set by the readout mode is used.
            The total exposure time is nints x inttime.
        start_time: string or float, optional
            The required clock time of the start of the exposure (MJD days).
            Strings other than 'NOW' are converted to floating point.
            If set to 'NOW', the current date-time is obtained from
            the system clock. By default, the internally stored clock
            time is used.
        wait_time: float, optional
            The wait time in seconds that has elapsed since the previous
            exposure. The default is 0.0 seconds.
        temperature: float, optional, default=detector target temperature
            The temperature of the detector, which will determine the dark
            current and readout noise.
        cosmic_ray_mode: string, optional
            The cosmic ray environment mode to be simulated. Available
            modes are:
        
            * 'NONE' - No cosmic rays.
            * 'SOLAR_MIN' - Solar minimum
            * 'SOLAR_MAX' - Solar maximum
            * 'SOLAR_FLARE' - Solar flare (worst case scenario)

            If this parameter is not explicitly given it will be obtained
            from the FITS header of the input file. Failing that, it will
            be obtained from the cosmic ray properties default.
        include_pixeldq: boolean, optional, default=True
            A flag that may be used to switch on the inclusion of the
            PIXELDQ data array in the output exposure data.
            By default, this array is included (and contains the bad pixel
            mask used to generate the simulated data).
        include_groupdq: boolean, optional, default=False
            A flag that may be used to switch on the inclusion of the
            GROUPDQ data array in the output exposure data.
            By default, this array is not included.
        qe_adjust: boolean, optional, default=True
            A flag that may be used to switch off quantum efficiency
            adjustment (for example to observe what effects in a simulation
            are caused by QE). *NOTE: When QE adjustment is turned off,
            the input illumination is assumed to be in electrons/second.*
        simulate_poisson_noise: boolean, optional, default=True
            A flag that may be used to switch off Poisson noise (for example
            to observe what effects in a simulation are caused by Poisson
            noise).
        simulate_read_noise: boolean, optional, default=True
            A flag that may be used to switch off read noise (for example
            to observe what effects in a simulation are caused by read
            noise).
        simulate_ref_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of reference
            pixels in the data. 
        simulate_bad_pixels: boolean, optional, default=True
            A flag that may be used to switch off the inclusion of bad
            pixels in the data, even if a bad pixel map containing bad pixels
            is specified in the detector properties.
        simulate_dark_current: boolean, optional, default=True
            A flag that may be used to switch off the addition of dark
            current (for example to observe what effects in a simulation are
            caused by dark current).
        simulate_flat_field: boolean, optional, default=True
            A flag that may be used to switch off the simulation of the pixel
            flat-field (for example to observe the effect of quantum
            efficiency in isolation).
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
        cdp_ftp_path: str, optional, default=None
            If specified, a list of folders (or folders) on the SFTP host
            where the MIRI CDPs are held to be searched, consisting of a
            list of folder names separated by a ":" delimiter.
            Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'
            If not specified, the default CDP repository at Leuven is used.
        cdp_ftp_path: str, optional, default=None
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
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator before
            generating the test data.
            If not specified, a value of None will be sent, which
            randomises the seed.
        verbose: int, optional, default=2
            Verbosity level. Activates print statements when non-zero.
        
            * 0 - no output at all except for error messages
            * 1 - warnings and errors only
            * 2 - normal output
            * 3, 4, 5 - additional commentary
            * 6, 7, 8 - extra debugging information
    
        :Returns:
    
        exposure_map: MiriExposureModel
            A MIRI exposure data model containing the simulated data.
    
        :Raises:
    
        ValueError
            Raised if any of the simulation parameters are out of range.
            Also raised if the value of a parameter is invalid or
            inappropriate.
        TypeError
            Raised if any of the simulation parameters are of the
            wrong type, size or shape.
        KeyError
            Raised if any of the simulation parameters do not contain
            a recognised keyword.
        IndexError
            Raised if an attempt is made to access beyond the bounds of
            the data.
        IOError
            Raised if a file cannot be read, interpreted or written.
        ImportError
            A delayed ImportError is raised if there is an attempt to use
            an optional library (such as matplotlib) which is not
            available. The software fails immediately if a compulsory
            library is not available.
        AttributeError
            Raised if simulation methods are executed in the wrong order,
            meaning that a necessary attribute is missing. A programming
            error.
        
        """
        if verbose > 1:
            _report_sca_parameters(illumination_map, '', '', '', fringemap,
                               '', '', qe_adjust,
                               simulate_poisson_noise, simulate_read_noise,
                               simulate_ref_pixels, simulate_bad_pixels,
                               simulate_dark_current, simulate_flat_field,
                               simulate_gain, simulate_nonlinearity,
                               simulate_drifts, simulate_latency,
                               cdp_ftp_path,
                               readnoise_version, bad_pixels_version,
                               flat_field_version, linearity_version,
                               gain_version,
                               logger=self.logger)

        # Set the seed for the np.random function.
        np.random.seed(seedvalue)

        # Extract the properties of the particular sensor chip assembly
        # being simulated.
        detectorid = illumination_map.meta.instrument.detector
        try:
            fpm = detector_properties.get('DETECTORS_DICT', str(detectorid))
        except KeyError:
            strg = "%s is not a known detector ID." % str(detectorid)
            raise KeyError(strg)

        # Get the metadata from the given data model
        description = 'Metadata extracted from %s' % illumination_map.get_title()
        metadata = Metadata(description)
        #intensity_metadata = Metadata('Intensity metadata')
        #wavelength_metadata = Metadata('Wavelength metadata')
        metadata.from_data_object(illumination_map, hdu_name='PRIMARY')
        #intensity_metadata.from_data_object(illumination_map,
        #                                    hdu_name='INTENSITY')
        #wavelength_metadata.from_data_object(illumination_map,
        #                                     hdu_name='WAVELENGTH')

        # If any of these arguments are explicitly set to None they need
        # to be defaulted. The default values are either extracted from
        # the metadata contained in the data object or, if not found, set
        # to a fixed default.
        (readout_mode, subarray, frame_time, temperature, cosmic_ray_mode) = \
            _get_parameter_defaults(fpm, metadata, readout_mode, subarray,
                                    frame_time, temperature, cosmic_ray_mode,
                                    verbose=verbose, logger=self.logger)

        # Set up the required detector parameters.
        self.setup(detectorid, readout_mode=readout_mode, subarray=subarray,
            burst_mode=burst_mode, inttime=inttime, ngroups=ngroups,
            nints=nints, temperature=temperature,
            cosmic_ray_mode=cosmic_ray_mode,
            include_pixeldq=include_pixeldq, include_groupdq=include_groupdq,
            qe_adjust=qe_adjust,
            simulate_poisson_noise=simulate_poisson_noise,
            simulate_read_noise=simulate_read_noise,
            simulate_ref_pixels=simulate_ref_pixels,
            simulate_bad_pixels=simulate_bad_pixels,
            simulate_dark_current=simulate_dark_current,
            simulate_flat_field=simulate_flat_field,
            simulate_gain=simulate_gain,
            simulate_nonlinearity=simulate_nonlinearity,
            simulate_drifts=simulate_drifts,
            simulate_latency=simulate_latency,
            cdp_ftp_path=cdp_ftp_path,
            readnoise_version=readnoise_version,
            bad_pixels_version=bad_pixels_version,
            flat_field_version=flat_field_version, 
            linearity_version=linearity_version, 
            gain_version=gain_version,
            makeplot=makeplot, verbose=verbose)
   
        # Set up the illumination data from illumination map object.
        self.set_illumination(illumination_map, scale=scale)    
        # Read the fringe map if a name has been provided.
        if fringemap is not None and fringemap:
            self.read_fringe_map(fringemap, ftype='FITS')

        # Wait for the given elapsed time.
        if wait_time > 0.0:
            self.wait(wait_time)
    
        # Simulate an exposure and return some simulated data.
        simulated_data = self.exposure(frame_time=frame_time,
                                       start_time=start_time)

        # Report the exposure times.
        if verbose > 1:
            (exposure_time, duration, start_time, mid_time, end_time) = \
                self.get_exposure_times()
            strg = "Exposure time %.2fs (duration %.2fs) " % (exposure_time,
                                                          duration)
            if verbose > 2:
                strg += "started at %.2f and finished at %.2f" % \
                    (start_time, end_time)
            self.logger.info( strg )
    
        # Plot the simulated data if requested.
        if makeplot:
            strg = "Simulated exposure data"
            self.plot(description=strg)
            if verbose > 3:
                # Development and testing - plot a ramp for the central pixel
                self.plot_ramp(simulated_data.shape[0]/2,
                               simulated_data.shape[1]/2)

        # Return the exposure data model.
        exposure_data = self.exposure_data
        return exposure_data
   
    def temperature_str(self, prefix=''):
        """
        
        Returns a string describing the detector temperature
        
        """
        strg = "%sDetector temperature = %.2f K " % (prefix, self.temperature)
        strg += "(which affects dark current and read noise)."
        return strg
    
    def cosmic_ray_str(self, prefix=''):
        """
        
        Returns a string describing the cosmic ray environment.
        
        """
        wordz = self.cosmic_ray_mode.split('+')
        if len(wordz) > 1:
            strg = "%sCosmic ray environment is %s (%s variant)." % \
                (prefix, wordz[0], wordz[1])
        else:
            strg = "%sCosmic ray environment is %s." % \
                (prefix, self.cosmic_ray_mode)
        return strg

    def plot(self, description=''):
        """
        
        Plot the simulated science data. Each set of integration data
        is plotted in a separate pyplot figure and the image read out
        from each group is plotted in a separate subplot.
        
        :Parameters:
        
        description: string, optional
            Additional description to be shown on the plot, if required.
            
        :Requires:
        
        matplotlib.pyplot
    
        ImportError
            Raised if the matplotlib plotting library is not available.
            
        """
        # There must be some exposure data to write out.
        if self.exposure_data is None:
            strg = "Exposure data not defined - " + \
                "use integration or exposure first."
            raise AttributeError(strg)
        
        # Plot the exposure data.
        self.exposure_data.plot(description=description)

    def plot_ramp(self, row, column, description='', frame_time=None,
                  show_ints=False):
        """
        
        Plot a ramp of signal vs integration and group at a particular
        row,column of the simulated data.
        
        :Parameters:
        
        row: int
            The row at which the ramp is to be plotted.
        column: int
            The column at which the ramp is to be plotted.
        description: string, optional
            Additional description to be shown on the plot, if required.
        frame_time: float, optional
            The detector frame time, in seconds.
            If specified, this parameter allows the frame time obtained
            from the detector properties to be overridden. This can be
            used, for example, to simulate the timing of full frame
            data with a subarray. If None, the frame time obtained from
            the detector properties for the current subarray and readout
            mode is used.
        show_ints: boolean, optional, default=False
            Set to True to distinguish the ramps belonging to different
            integrations. Only works when a single row and column is
            specified.
            
        """
        # There must be some exposure data to write out.
        if self.exposure_data is None:
            strg = "Exposure data not defined "
            strg += "- use integration or exposure first."
            raise AttributeError(strg)
        
        # Plot the exposure ramp.
        (exp, elapsed) = self.detector.exposure_time(
                                        self.nints,
                                        self.ngroups,
                                        self.samplesum,
                                        self.sampleskip,
                                        refpixsampleskip=self.refpixsampleskip,
                                        nframes=self.nframes,
                                        groupgap=self.groupgap,
                                        subarray=self.subarray,
                                        burst_mode=self.subarray_burst_mode,
                                        frame_time=frame_time)
        stime = elapsed / (self.nints * self.ngroups)
        self.exposure_data.plot_ramp(row, column, stime=stime,
                                     tunit='seconds', show_ints=show_ints,
                                     description=description)

@six.add_metaclass(Singleton)
class SensorChipAssembly1(SensorChipAssembly):
    # Only one instance of this class is allowed to exist.
    """

    Class Sensor Chip Assembly - Simulates the behaviour of the MIRI
    detectors. The class converts detector illumination data (typically
    provided in a FITS file generated by another MIRI simulator and read
    by the read_data method) and generates simulated detector data,
    which may be written to a FITSWriter or level 1 FITS file by the
    write_data method.
    
    Singleton instance 1

    :Parameters:
    
    See setup method.   
                
    """
    def __init__(self, logger=LOGGER):
        """
         
        Constructor for class SensorChipAssembly1.
         
        :Parameters:
         
        None. The class is a singleton and the constructor is only called once.
         
        """
        super(SensorChipAssembly1, self).__init__(logger=logger)

# Only one instance of this class is allowed to exist.
@six.add_metaclass(Singleton)
class SensorChipAssembly2(SensorChipAssembly):
    """

    Class Sensor Chip Assembly - Simulates the behaviour of the MIRI
    detectors. The class converts detector illumination data (typically
    provided in a FITS file generated by another MIRI simulator and read
    by the read_data method) and generates simulated detector data,
    which may be written to a FITSWriter or level 1 FITS file by the
    write_data method.
    
    Singleton instance 2

    :Parameters:
    
    See setup method.   
                
    """
    def __init__(self, logger=LOGGER):
        """
         
        Constructor for class SensorChipAssembly2.
         
        :Parameters:
         
        None. The class is a singleton and the constructor is only called once.
         
        """
        super(SensorChipAssembly2, self).__init__(logger=logger)


# Only one instance of this class is allowed to exist.
@six.add_metaclass(Singleton)
class SensorChipAssembly3(SensorChipAssembly):
    """

    Class Sensor Chip Assembly - Simulates the behaviour of the MIRI
    detectors. The class converts detector illumination data (typically
    provided in a FITS file generated by another MIRI simulator and read
    by the read_data method) and generates simulated detector data,
    which may be written to a FITSWriter or level 1 FITS file by the
    write_data method.
    
    Singleton instance 3

    :Parameters:
    
    See setup method.   
                
    """
    def __init__(self, logger=LOGGER):
        """
         
        Constructor for class SensorChipAssembly3.
         
        :Parameters:
         
        None. The class is a singleton and the constructor is only called once.
         
        """
        super(SensorChipAssembly3, self).__init__(logger=logger)

#
# Public simulator functions
#
def simulate_sca(inputfile, outputfile, detectorid, scale=1.0, fringemap=None,
                 readout_mode=None, subarray=None, burst_mode=True,
                 frame_time=None, inttime=None, ngroups=None, nints=None,
                 start_time=None, wait_time=0.0, temperature=None,
                 cosmic_ray_mode=None, fileformat='STScI',
                 datashape='hypercube', include_pixeldq=True,
                 include_groupdq=False, overwrite=False, qe_adjust=True,
                 simulate_poisson_noise=True, simulate_read_noise=True,
                 simulate_ref_pixels=True, simulate_bad_pixels=True,
                 simulate_dark_current=True, simulate_flat_field=True,
                 simulate_gain=True, simulate_nonlinearity=True,
                 simulate_drifts=True, simulate_latency=True,
                 cdp_ftp_path=SIM_CDP_FTP_PATH,
                 readnoise_version='', bad_pixels_version='',
                 flat_field_version='', linearity_version='', gain_version='',
                 makeplot=False, seedvalue=None, verbose=2, logger=LOGGER):
    """
    
    Runs the MIRI SCA simulator on the given input file, generating the
    given output file.
                 
    :Parameters:
    
    inputfile: string
        The name of the FITS file containing detector illumination data
        (normally created by another MIRI simulator).
    outputfile: string
        The name of the FITSWriter or level 1 FITS file to contain the
        simulated SCA data.
    detectorid: string
        Detector ID, identifying a particular detector.
        The MIRI instrument has three detectors: 'MIRIMAGE',
        'MIRIFULONG' and 'MIRIFUSHORT'. (These correspond to the imager
        and the long and short wavelength arms of the MRS respectively.)
    scale: float, optional, default=1.0
        An optional scale factor to apply to the intensity data.
        This is for debugging and testing only - the input data
        should already be scaled to photons/second/pixel.
    fringemap: string, optional
        The name of a file containing a fringe map to be used to
        modulate the input illumination. If not specified, no fringe
        map will be used. 
    readout_mode: string, optional
        Readout mode, which can be one of the following common options:
        
        * 'SLOW' - 10 samples per readout and defaults of ngroups=10
          and nints=1.
        * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
          nints=10.
        * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
          averaged.
        * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
          averaged.
          
        The following unusual options are also available for testing:
        
        * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
          4 frames per group and a gap of 8 frames between each group.
        * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
          averaged and default ngroups=1.
        * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
          averaged and default ngroups=4.
        * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
          4 frames per group and a gap of 8 frames between each group.
          
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the detector properties default.
    subarray: string, optional
        Detector subarray mode for output. This can be one 'FULL',
        'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
        'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
        full-frame and the other modes read out portions of the detector
        as described in the MIRI Operational Concept Document.
        The simulator can also accept the other subarray modes defined
        in detector_properties.py for test purposes.
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file (unless the input data
        specifies a non-standard subarray). Failing that, it will be
        obtained from the detector properties default (FULL).
    frame_time: float, optional
        The detector frame time, in seconds.
        If specified, this parameter overrides the frame time obtained
        from the detector readout mode and subarray. This can be
        used, for example, to simulate the timing of full frame
        exposures with subarray-sized data (to save memory). If None,
        the frame time obtained from the current readout mode and
        subarray is used.
    inttime: float, optional
        The integration time in seconds to be simulated. This parameter
        will be ignored if ngroups is provided as well.
        *NOTE: Any requested integration time will be rounded up to
        the time resulting from the nearest whole number of groups.*
    ngroups: int, optional
        The number of readout groups. If None, this is derived from
        the integration time and readout mode. Otherwise the parameter
        overrides the integration time and sets the number of groups
        directly. There must be at least one group.
    nints: int, optional
        The number of integrations per exposure. It must be at least 1.
        If None, the default set by the readout mode is used.
        The total exposure time is nints x inttime.
    start_time: string or float, optional
        The required clock time of the start of the exposure (MJD days).
        Strings other than 'NOW' are converted to floating point.
        If set to 'NOW', the current date-time is obtained from
        the system clock. By default, the internally stored clock
        time is used.
    wait_time: float, optional
        The wait time in seconds that has elapsed since the previous
        exposure. The default is 0.0 seconds.
    temperature: float, optional, default=detector target temperature
        The temperature of the detector, which will determine the dark
        current and readout noise.
    cosmic_ray_mode: string, optional, default=cosmic ray properties
        The cosmic ray environment mode to be simulated. Available
        modes are:
        
        * 'NONE' - No cosmic rays.
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)

        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the cosmic ray properties default.
    fileformat: string, optional, default='STScI'
        The kind of file format to be written.
            
        * 'STScI' - use the STScI level 1 model for the JWST
          DMS pipeline.
        * 'FITSWriter' - emulate the format written by the FITSWriter
          during MIRI VM, FM and CV tests and read by the DHAS.
              
    datashape: string, optional, default='hypercube'
        The SCI data shape to be written.
            
        * 'hypercube' - write the SCI data to a 4 dimensional FITS image
          with separate columns x rows x groups x integrations
          dimensions.
        * 'cube' - append the groups and integrations to make a 3
          dimensional FITS image with columns x rows x (groups and
          integrations) dimensions.
          
    include_pixeldq: boolean, optional, default=True
        A flag that may be used to switch on the inclusion of the
        PIXELDQ data array in the output exposure data.
        By default, this array is included (and contains the bad pixel
        mask used to generate the simulated data).
    include_groupdq: boolean, optional, default=False
        A flag that may be used to switch on the inclusion of the
        GROUPDQ data array in the output exposure data.
        By default, this array is not included.
    overwrite: bool, optional, default=False
        Parameter passed to pyfits.HDUlist.writeto
    qe_adjust: boolean, optional, default=True
        A flag that may be used to switch off quantum efficiency
        adjustment (for example to observe what effects in a simulation
        are caused by QE). *NOTE: When QE adjustment is turned off,
            the input illumination is assumed to be in electrons/second.*
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for example
        to observe what effects in a simulation are caused by Poisson
        noise).
    simulate_read_noise: boolean, optional, default=True
        A flag that may be used to switch off read noise (for example
        to observe what effects in a simulation are caused by read
        noise).
    simulate_ref_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of reference
        pixels in the data. 
    simulate_bad_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of bad
        pixels in the data, even if a bad pixel map containing bad pixels
        is specified in the detector properties.
    simulate_dark_current: boolean, optional, default=True
        A flag that may be used to switch off the addition of dark
        current (for example to observe what effects in a simulation are
        caused by dark current).
    simulate_flat_field: boolean, optional, default=True
        A flag that may be used to switch off the simulation of the pixel
        flat-field (for example to observe the effect of quantum
        efficiency in isolation).
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
    cdp_ftp_path: str, optional, default=None
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
    seedvalue: int, optional, default=None
        The seed to be sent to the np.random number generator before
        generating the test data.
        If not specified, a value of None will be sent, which
        randomises the seed.
    verbose: int, optional, default=2
        Verbosity level. Activates print statements when non-zero.
        
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
        Raised if any of the simulation parameters are out of range.
        Also raised if the value of a parameter is invalid or
        inappropriate.
    TypeError
        Raised if any of the simulation parameters are of the
        wrong type, size or shape.
    KeyError
        Raised if any of the simulation parameters do not contain
        a recognised keyword.
    IndexError
        Raised if an attempt is made to access beyond the bounds of
        the data.
    IOError
        Raised if a file cannot be read, interpreted or written.
    ImportError
        A delayed ImportError is raised if there is an attempt to use
        an optional library (such as matplotlib) which is not
        available. The software fails immediately if a compulsory
        library is not available.
    AttributeError
        Raised if simulation methods are executed in the wrong order,
        meaning that a necessary attribute is missing. A programming
        error.
        
    """
    strg = "\nsimulate_sca function is now deprecated. "
    strg += "Please use the SensorChipAssembly.simulate_files() class method."
    logger.warn(strg)
    # Create an sca object and then run the simuation.
    sca = SensorChipAssembly1()
    sca.simulate_files(inputfile, outputfile, detectorid, scale=scale,
        fringemap=fringemap, readout_mode=readout_mode, subarray=subarray,
        burst_mode=burst_mode,
        frame_time=frame_time, inttime=inttime, ngroups=ngroups, nints=nints,
        start_time=start_time, wait_time=wait_time, temperature=temperature,
        cosmic_ray_mode=cosmic_ray_mode, fileformat=fileformat,
        datashape=datashape, include_pixeldq=include_pixeldq,
        include_groupdq=include_groupdq, overwrite=overwrite, qe_adjust=qe_adjust,
        simulate_poisson_noise=simulate_poisson_noise,
        simulate_read_noise=simulate_read_noise,
        simulate_ref_pixels=simulate_ref_pixels,
        simulate_bad_pixels=simulate_bad_pixels,
        simulate_dark_current=simulate_dark_current,
        simulate_flat_field=simulate_flat_field,
        simulate_gain=simulate_gain,
        simulate_nonlinearity=simulate_nonlinearity,
        simulate_drifts=simulate_drifts,
        simulate_latency=simulate_latency,
        cdp_ftp_path=cdp_ftp_path,
        readnoise_version=readnoise_version,
        bad_pixels_version=bad_pixels_version,
        flat_field_version=flat_field_version, 
        linearity_version=linearity_version, 
        gain_version=gain_version,
        makeplot=makeplot, seedvalue=seedvalue, verbose=verbose)

def simulate_sca_list(inputfile, outputfile, detectorid, scale=1.0,
                      fringemap=None,
                      readout_mode=None, subarray=None, burst_mode=True,
                      frame_time=None, inttime=None, ngroups=None, nints=None,
                      start_time=None, wait_time=0.0, temperature=None,
                      cosmic_ray_mode=None, fileformat='STScI', 
                      datashape='hypercube', overwrite=False,
                      include_pixeldq=True, include_groupdq=False,
                      qe_adjust=True,
                      simulate_poisson_noise=True, simulate_read_noise=True,
                      simulate_ref_pixels=True, simulate_bad_pixels=True,
                      simulate_dark_current=True, simulate_flat_field=True,
                      simulate_gain=True, simulate_nonlinearity=True,
                      simulate_drifts=True, simulate_latency=True,
                      cdp_ftp_path=SIM_CDP_FTP_PATH,
                      readnoise_version='', bad_pixels_version='',
                      flat_field_version='', linearity_version='',
                      gain_version='',
                      makeplot=False, seedvalue=None, verbose=2, logger=LOGGER):
    """
    
    Runs the MIRI SCA simulator on the given list of input files,
    generating the given list of output files. Persistence effects
    are simulated between one exposure and the next.
                 
    :Parameters:
    
    inputfile: string or list of strings
        The name of the FITS file containing detector illumination data
        (normally created by another MIRI simulator).
    outputfile: string or list of strings
        The name of the FITSWriter or level 1 FITS file to contain the
        simulated SCA data. If an entry in the list is empty, no file
        will be written.
    detectorid: string
        The detector ID, identifying a particular detector.
        The MIRI instrument has three detectors: 'MIRIMAGE',
        'MIRIFULONG' and 'MIRIFUSHORT'. (These correspond to the imager
        and the long and short wavelength arms of the MRS respectively.)
    scale: float, optional, default=1.0
        An optional scale factor to apply to the intensity data.
        This is for debugging and testing only - the input data
        should already be scaled to photons/second/pixel.
    fringemap: string, optional
        The name of a file containing a fringe map to be used to
        modulate the input illumination. If not specified, no fringe
        map will be used. 
    readout_mode: string, optional
        Readout mode, which can be one of the following common options:
        
        * 'SLOW' - 10 samples per readout and defaults of ngroups=10
          and nints=1.
        * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
          nints=10.
        * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
          averaged.
        * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
          averaged.
          
        The following unusual options are also available for testing:
        
        * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
          4 frames per group and a gap of 8 frames between each group.
        * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
          averaged and default ngroups=1.
        * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
          averaged and default ngroups=4.
        * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
          4 frames per group and a gap of 8 frames between each group.
          
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the detector properties default.
    subarray: string, optional
        Detector subarray mode for output. This can be one 'FULL',
        'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
        'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
        full-frame and the other modes read out portions of the detector
        as described in the MIRI Operational Concept Document.
        The simulator can also accept the other subarray modes defined
        in detector_properties.py for test purposes.
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file (unless the input data
        specifies a non-standard subarray). Failing that, it will be
        obtained from the detector properties default (FULL).
    frame_time: float, optional
        The detector frame time, in seconds.
        If specified, this parameter overrides the frame time obtained
        from the detector readout mode and subarray. This can be
        used, for example, to simulate the timing of full frame
        exposures with subarray-sized data (to save memory). If None,
        the frame time obtained from the current readout mode and
        subarray is used.
    inttime: float, optional
        The integration time in seconds to be simulated. This parameter
        will be ignored if ngroups is provided as well.
        *NOTE: Any requested integration time will be rounded up to
        the time resulting from the nearest whole number of groups.*
    ngroups: int, optional
        The number of readout groups. If None, this is derived from
        the integration time and readout mode. Otherwise the parameter
        overrides the integration time and sets the number of groups
        directly. There must be at least one group.
    nints: int, optional
        The number of integrations per exposure. It must be at least 1.
        If None, the default set by the readout mode is used.
        The total exposure time is nints x inttime.
    start_time: string or float, optional
        The required clock time of the start of the first exposure (MJD days).
        Strings other than 'NOW' are converted to floating point.
        If set to 'NOW', the current date-time is obtained from
        the system clock. By default, the internally stored clock
        time is used.
    wait_time: float, optional
        The wait time in seconds elapsed in between exposures.
        The default is 0.0 seconds.
    temperature: float, optional, default=detector target temperature
        The temperature of the detector, which will determine the dark
        current and readout noise.
    cosmic_ray_mode: string, optional, default=cosmic ray properties
        The cosmic ray environment mode to be simulated. Available
        modes are:
        
        * 'NONE' - No cosmic rays.
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)

        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the cosmic ray properties default.
    fileformat: string, optional, default='STScI'
        The kind of file format to be written.
            
        * 'STScI' - use the STScI level 1 model for the JWST
          DMS pipeline.
        * 'FITSWriter' - emulate the format written by the FITSWriter
          during MIRI VM, FM and CV tests and read by the DHAS.
              
    datashape: string, optional, default='hypercube'
        The SCI data shape to be written.
            
        * 'hypercube' - write the SCI data to a 4 dimensional FITS image
          with separate columns x rows x groups x integrations
          dimensions.
        * 'cube' - append the groups and integrations to make a 3
          dimensional FITS image with columns x rows x (groups and
          integrations) dimensions.
          
    include_pixeldq: boolean, optional, default=True
        A flag that may be used to switch on the inclusion of the
        PIXELDQ data array in the output exposure data.
        By default, this array is included (and contains the bad pixel
        mask used to generate the simulated data).
    include_groupdq: boolean, optional, default=False
        A flag that may be used to switch on the inclusion of the
        GROUPDQ data array in the output exposure data.
        By default, this array is not included.
    overwrite: bool, optional, default=False
        Parameter passed to pyfits.HDUlist.writeto
    qe_adjust: boolean, optional, default=True
        A flag that may be used to switch off quantum efficiency
        adjustment (for example to observe what effects in a simulation
        are caused by QE). *NOTE: When QE adjustment is turned off,
            the input illumination is assumed to be in electrons/second.*
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for example
        to observe what effects in a simulation are caused by Poisson
        noise).
    simulate_read_noise: boolean, optional, default=True
        A flag that may be used to switch off read noise (for example
        to observe what effects in a simulation are caused by read
        noise).
    simulate_ref_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of reference
        pixels in the data. 
    simulate_bad_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of bad
        pixels in the data, even if a bad pixel map containing bad pixels
        is specified in the detector properties.
    simulate_dark_current: boolean, optional, default=True
        A flag that may be used to switch off the addition of dark
        current (for example to observe what effects in a simulation are
        caused by dark current).
    simulate_flat_field: boolean, optional, default=True
        A flag that may be used to switch off the simulation of the pixel
        flat-field (for example to observe the effect of quantum
        efficiency in isolation).
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
    cdp_ftp_path: str, optional, default=None
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
    seedvalue: int, optional, default=None
        The seed to be sent to the np.random number generator before
        generating the test data.
        If not specified, a value of None will be sent, which
        randomises the seed.
    verbose: int, optional, default=2
        Verbosity level. Activates print statements when non-zero.
        
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
        Raised if any of the simulation parameters are out of range.
        Also raised if the value of a parameter is invalid or
        inappropriate.
    TypeError
        Raised if any of the simulation parameters are of the
        wrong type, size or shape.
    KeyError
        Raised if any of the simulation parameters do not contain
        a recognised keyword.
    IndexError
        Raised if an attempt is made to access beyond the bounds of
        the data.
    IOError
        Raised if a file cannot be read, interpreted or written.
    ImportError
        A delayed ImportError is raised if there is an attempt to use
        an optional library (such as matplotlib) which is not
        available. The software fails immediately if a compulsory
        library is not available.
    AttributeError
        Raised if simulation methods are executed in the wrong order,
        meaning that a necessary attribute is missing. A programming
        error.
        
    """
    strg = "\nsimulate_sca_list function is now deprecated. "
    strg += "Please use the SensorChipAssembly.simulate_files() class method."
    logger.warn(strg)
    # Create an sca object and then run the simuation.
    sca = SensorChipAssembly1()
    sca.simulate_files(inputfile, outputfile, detectorid, scale=scale,
        fringemap=fringemap, readout_mode=readout_mode, subarray=subarray,
        burst_mode=burst_mode,
        frame_time=frame_time, inttime=inttime, ngroups=ngroups, nints=nints,
        start_time=start_time, wait_time=wait_time, temperature=temperature,
        cosmic_ray_mode=cosmic_ray_mode, fileformat=fileformat,
        datashape=datashape, overwrite=overwrite, include_pixeldq=include_pixeldq,
        include_groupdq=include_groupdq, qe_adjust=qe_adjust,
        simulate_poisson_noise=simulate_poisson_noise,
        simulate_read_noise=simulate_read_noise,
        simulate_ref_pixels=simulate_ref_pixels,
        simulate_bad_pixels=simulate_bad_pixels,
        simulate_dark_current=simulate_dark_current,
        simulate_flat_field=simulate_flat_field,
        simulate_gain=simulate_gain,
        simulate_nonlinearity=simulate_nonlinearity,
        simulate_drifts=simulate_drifts,
        simulate_latency=simulate_latency,
        cdp_ftp_path=cdp_ftp_path,
        readnoise_version=readnoise_version,
        bad_pixels_version=bad_pixels_version,
        flat_field_version=flat_field_version, 
        linearity_version=linearity_version, 
        gain_version=gain_version,
        makeplot=makeplot, seedvalue=seedvalue, verbose=verbose)
   
def simulate_sca_pipeline(illumination_map, scale=1.0,
                          fringemap=None, readout_mode=None, subarray=None,
                          burst_mode=True, frame_time=None,
                          inttime=None, ngroups=None, nints=None,
                          start_time=None, wait_time=0.0, temperature=None,
                          cosmic_ray_mode=None, overwrite=False,
                          include_pixeldq=True, include_groupdq=False,
                          qe_adjust=True,
                          simulate_poisson_noise=True, simulate_read_noise=True,
                          simulate_ref_pixels=True, simulate_bad_pixels=True,
                          simulate_dark_current=True, simulate_flat_field=True,
                          simulate_gain=True, simulate_nonlinearity=True,
                          simulate_drifts=True, simulate_latency=True,
                          cdp_ftp_path=SIM_CDP_FTP_PATH,
                          readnoise_version='', bad_pixels_version='',
                          flat_field_version='', linearity_version='',
                          gain_version='',
                          makeplot=False, seedvalue=None, verbose=2,
                          logger=LOGGER):
    """
    
    Runs the MIRI SCA simulator on the given MiriIluminationModel
    and returns a MiriExposureModel
    This is an alternative to the file-based functions simulate_sca()
    and simulate_sca_list(), which can be used by pipeline code which
    needs to convert one data model into another.

        
    NOTE: The clear_exposure_data() method can be used to delete
    the internal reference to the exposure data and save memory in
    between simulations.
                 
    :Parameters:
    
    illumination_map: MiriIlluminationModel object
        A MIRI illumination map object which describes the input
        illumination.
    scale: float, optional, default=1.0
        An optional scale factor to apply to the intensity data.
        This is for debugging and testing only - the input data
        should already be scaled to photons/second/pixel.
    fringemap: string, optional
        The name of a file containing a fringe map to be used to
        modulate the input illumination. If not specified, no fringe
        map will be used. 
    readout_mode: string, optional
        Readout mode, which can be one of the following common options:
        
        * 'SLOW' - 10 samples per readout and defaults of ngroups=10
          and nints=1.
        * 'FAST' - 1 sample per readout and defaults of ngroups=1 and
          nints=10.
        * 'FASTINTAVG' - same as FAST but with groups of 4 integrations
          averaged.
        * 'FASTGRPAVG' - same as FAST but with groups of 4 groups
          averaged.
          
        The following unusual options are also available for testing:
        
        * 'FASTGRPGAP' - a non-MIRI mode similar to FAST mode but with
          4 frames per group and a gap of 8 frames between each group.
        * 'SLOWINTAVG' - same as SLOW but with groups of 4 integrations
          averaged and default ngroups=1.
        * 'SLOWGRPAVG' - same as SLOW but with groups of 4 groups
          averaged and default ngroups=4.
        * 'SLOWGRPGAP' - a non-MIRI mode similar to SLOW mode but with
          4 frames per group and a gap of 8 frames between each group.
        
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the detector properties default.
    subarray: string, optional
        Detector subarray mode for output. This can be one 'FULL',
        'MASK1550', 'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',
        'SUB256', 'SUB128', 'SUB64' or 'SLITLESSPRISM', etc. 'FULL' is
        full-frame and the other modes read out portions of the detector
        as described in the MIRI Operational Concept Document.
        The simulator can also accept the other subarray modes defined
        in detector_properties.py for test purposes.
        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file (unless the input data
        specifies a non-standard subarray). Failing that, it will be
        obtained from the detector properties default (FULL).
    frame_time: float, optional
        The detector frame time, in seconds.
        If specified, this parameter overrides the frame time obtained
        from the detector readout mode and subarray. This can be
        used, for example, to simulate the timing of full frame
        exposures with subarray-sized data (to save memory). If None,
        the frame time obtained from the current readout mode and
        subarray is used.
    inttime: float, optional
        The integration time in seconds to be simulated. This parameter
        will be ignored if ngroups is provided as well.
        *NOTE: Any requested integration time will be rounded up to
        the time resulting from the nearest whole number of groups.*
    ngroups: int, optional
        The number of readout groups. If None, this is derived from
        the integration time and readout mode. Otherwise the parameter
        overrides the integration time and sets the number of groups
        directly. There must be at least one group.
    nints: int, optional
        The number of integrations per exposure. It must be at least 1.
        If None, the default set by the readout mode is used.
        The total exposure time is nints x inttime.
    start_time: string or float, optional
        The required clock time of the start of the exposure (MJD days).
        Strings other than 'NOW' are converted to floating point.
        If set to 'NOW', the current date-time is obtained from
        the system clock. By default, the internally stored clock
        time is used.
    wait_time: float, optional
        The wait time in seconds that has elapsed since the previous
        exposure. The default is 0.0 seconds.
    temperature: float, optional, default=detector target temperature
        The temperature of the detector, which will determine the dark
        current and readout noise.
    cosmic_ray_mode: string, optional
        The cosmic ray environment mode to be simulated. Available
        modes are:
        
        * 'NONE' - No cosmic rays.
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)

        If this parameter is not explicitly given it will be obtained
        from the FITS header of the input file. Failing that, it will
        be obtained from the cosmic ray properties default.
    include_pixeldq: boolean, optional, default=True
        A flag that may be used to switch on the inclusion of the
        PIXELDQ data array in the output exposure data.
        By default, this array is included (and contains the bad pixel
        mask used to generate the simulated data).
    include_groupdq: boolean, optional, default=False
        A flag that may be used to switch on the inclusion of the
        GROUPDQ data array in the output exposure data.
        By default, this array is not included.
    qe_adjust: boolean, optional, default=True
        A flag that may be used to switch off quantum efficiency
        adjustment (for example to observe what effects in a simulation
        are caused by QE). *NOTE: When QE adjustment is turned off,
        the input illumination is assumed to be in electrons/second.*
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for example
        to observe what effects in a simulation are caused by Poisson
        noise).
    simulate_read_noise: boolean, optional, default=True
        A flag that may be used to switch off read noise (for example
        to observe what effects in a simulation are caused by read
        noise).
    simulate_ref_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of reference
        pixels in the data. 
    simulate_bad_pixels: boolean, optional, default=True
        A flag that may be used to switch off the inclusion of bad
        pixels in the data, even if a bad pixel map containing bad pixels
        is specified in the detector properties.
    simulate_dark_current: boolean, optional, default=True
        A flag that may be used to switch off the addition of dark
        current (for example to observe what effects in a simulation are
        caused by dark current).
    simulate_flat_field: boolean, optional, default=True
        A flag that may be used to switch off the simulation of the pixel
        flat-field (for example to observe the effect of quantum
        efficiency in isolation).
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
    cdp_ftp_path: str, optional, default=None
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
    seedvalue: int, optional, default=None
        The seed to be sent to the np.random number generator before
        generating the test data.
        If not specified, a value of None will be sent, which
        randomises the seed.
    verbose: int, optional, default=2
        Verbosity level. Activates print statements when non-zero.
        
        * 0 - no output at all except for error messages
        * 1 - warnings and errors only
        * 2 - normal output
        * 3, 4, 5 - additional commentary
        * 6, 7, 8 - extra debugging information

    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
     
    :Returns:
    
    exposure_map: MiriExposureModel
        A MIRI exposure data model containing the simulated data.
    
    :Raises:
    
    ValueError
        Raised if any of the simulation parameters are out of range.
        Also raised if the value of a parameter is invalid or
        inappropriate.
    TypeError
        Raised if any of the simulation parameters are of the
        wrong type, size or shape.
    KeyError
        Raised if any of the simulation parameters do not contain
        a recognised keyword.
    IndexError
        Raised if an attempt is made to access beyond the bounds of
        the data.
    IOError
        Raised if a file cannot be read, interpreted or written.
    ImportError
        A delayed ImportError is raised if there is an attempt to use
        an optional library (such as matplotlib) which is not
        available. The software fails immediately if a compulsory
        library is not available.
    AttributeError
        Raised if simulation methods are executed in the wrong order,
        meaning that a necessary attribute is missing. A programming
        error.
        
    """
    strg = "\nsimulate_sca_pipeline function is now deprecated. "
    strg += "Please use the SensorChipAssembly.simulate_pipe() class method."
    logger.warn(strg)
    # Create an sca object and then run the simuation.
    sca = SensorChipAssembly1()
    exposure_data = sca.simulate_pipe(illumination_map, scale=scale,
        fringemap=fringemap, readout_mode=readout_mode, subarray=subarray,
        burst_mode=burst_mode,
        frame_time=frame_time, inttime=inttime, ngroups=ngroups, nints=nints,
        start_time=start_time, wait_time=wait_time, temperature=temperature,
        cosmic_ray_mode=cosmic_ray_mode, fileformat=fileformat,
        datashape=datashape, overwrite=overwrite, include_pixeldq=include_pixeldq,
        include_groupdq=include_groupdq, qe_adjust=qe_adjust,
        simulate_poisson_noise=simulate_poisson_noise,
        simulate_read_noise=simulate_read_noise,
        simulate_ref_pixels=simulate_ref_pixels,
        simulate_bad_pixels=simulate_bad_pixels,
        simulate_dark_current=simulate_dark_current,
        simulate_flat_field=simulate_flat_field,
        simulate_gain=simulate_gain,
        simulate_nonlinearity=simulate_nonlinearity,
        simulate_drifts=simulate_drifts,
        simulate_latency=simulate_latency,
        cdp_ftp_path=cdp_ftp_path,
        readnoise_version=readnoise_version,
        bad_pixels_version=bad_pixels_version,
        flat_field_version=flat_field_version, 
        linearity_version=linearity_version, 
        gain_version=gain_version,
        makeplot=makeplot, seedvalue=seedvalue, verbose=verbose)
    return exposure_data


#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_sensor_chip_assembly module.
# However, the tests here are made in verbose mode and include more file
# I/O and plotting than the unit tests, so they show what is happening
# in more detail.
#
if __name__ == '__main__':
    print( "Testing the Sensor Chip Assembly class" )
    
    import random
    random.seed()
    
    test_input_file_name  = './data/SCATestHorseHead1024.fits'
#     test_input_file_name  = './data/SCATestInput80x64.fits'
#    test_input_file_name  = './data/DoesNotExist.fits'
    test_output_stub = './data/SCATestOutput'

    # Which tests to run?
    TEST_MODELS = True
    TEST_BANDS = True
    TEST_SUBARRAYS = True
    TEST_READOUTS = True
    
    # Which simulations to include?
    QE_ADJUST = True,
    SIMULATE_POISSON_NOISE = True
    SIMULATE_READ_NOISE = True
    SIMULATE_REF_PIXELS = True
    SIMULATE_BAD_PIXELS = True
    SIMULATE_DARK_CURRENT = True
    SIMULATE_FLAT_FIELD = True
    SIMULATE_GAIN = True
    SIMULATE_NONLINEARITY = True
    SIMULATE_DRIFTS = True
    SIMULATE_LATENCY = True
    
    # MODIFY THESE TWO VARIABLES TO CONTROL THE DEGREE OF INTERACTION
    VERBOSE = 1
    PLOTTING = False
    INCLUDE_PIXELDQ = True

    if TEST_MODELS:
        # ------------------------------------------------------------------
        # Test a simulation with the intensity array defined directly.
        print("\nTesting simulate_pipe from data ...\n" + (50 * "*"))
        list1 = (100000,0,100000,0,100000,0,100000,0,100000,0,100000,0,100000,0,100000,0)
        list2 = (list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1,list1)
        a = np.array(list2)
         
        # The correct way to invoke SCASim from a data array is now to convert
        # that data array into a data model and use simulate_pipe
        illumination_map = MiriIlluminationModel(intensity=a)
        illumination_map.set_instrument_metadata('MIRIFULONG')
        illumination_map.set_wcs_metadata(wcsaxes=2, cunit=['pixels','pixels'])
        if VERBOSE > 1:
            print( illumination_map )
        sca = SensorChipAssembly3()
        exposure_map = sca.simulate_pipe(illumination_map, scale=10.0,
                            readout_mode='SLOW', nints=1, ngroups=10,
                            start_time=42.0, makeplot=PLOTTING,
                            cosmic_ray_mode='NONE',
                            include_pixeldq=INCLUDE_PIXELDQ, qe_adjust=QE_ADJUST,
                            simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                            simulate_read_noise=SIMULATE_READ_NOISE,
                            simulate_ref_pixels=SIMULATE_REF_PIXELS,
                            simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                            simulate_dark_current=SIMULATE_DARK_CURRENT,
                            simulate_flat_field=SIMULATE_FLAT_FIELD,
                            simulate_gain=SIMULATE_GAIN,
                            simulate_nonlinearity=SIMULATE_NONLINEARITY,
                            simulate_drifts=SIMULATE_DRIFTS,
                            simulate_latency=SIMULATE_LATENCY,      
                            verbose=VERBOSE)
        if VERBOSE > 1:
            print( exposure_map )
            (exposure_time, duration, start_time, mid_time, end_time) = \
                exposure_map.get_exposure_times()
            print("Exposure time %.2f, duration %.2f, start %.2f, end %.2f" % \
                (exposure_time, duration, start_time, end_time))
        exposure_map.save('./data/SCATestOutputFromData.fits', overwrite=True)
      
        # ------------------------------------------------------------------
        # Test a simulation with an intensity model to be
        # converted to an exposure model.
        print("\nTesting simulate_pipe from file...\n" + (50 * "*"))
        illumination_map = MiriIlluminationModel(test_input_file_name)
        illumination_map.set_instrument_metadata('MIRIFUSHORT')
        illumination_map.set_wcs_metadata(wcsaxes=2, cunit=['pixels','pixels'])
        if VERBOSE > 1:
            print( illumination_map )
        exposure_map = sca.simulate_pipe(illumination_map, scale=10.0,
                            readout_mode='SLOW', nints=1, ngroups=10,
                            makeplot=PLOTTING, cosmic_ray_mode='NONE',
                            include_pixeldq=INCLUDE_PIXELDQ,
                            qe_adjust=QE_ADJUST,
                            simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                            simulate_read_noise=SIMULATE_READ_NOISE,
                            simulate_ref_pixels=SIMULATE_REF_PIXELS,
                            simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                            simulate_dark_current=SIMULATE_DARK_CURRENT,
                            simulate_flat_field=SIMULATE_FLAT_FIELD,
                            simulate_gain=SIMULATE_GAIN,
                            simulate_nonlinearity=SIMULATE_NONLINEARITY,
                            simulate_drifts=SIMULATE_DRIFTS,
                            verbose=VERBOSE)
        if VERBOSE > 1:
            print( exposure_map )
            (exposure_time, duration, start_time, mid_time, end_time) = \
                exposure_map.get_exposure_times()
            print("Exposure time %.2f, duration %.2f, start %.2f, end %.2f" % \
                (exposure_time, duration, start_time, end_time))
        exposure_map.save('./data/SCATestOutputFromObject.fits', overwrite=True)
        # Remove the internal reference to the exposure data
        sca.clear_exposure_data()
        del illumination_map, exposure_map

    if TEST_BANDS:
        # ------------------------------------------------------------------
        print("\nTesting simulate_pipe with all MRS bands...\n" + (50 * "*"))
        mode = 'SLOW'
        ngroups = 4
        del sca
        sca = SensorChipAssembly2()
        ii = 1
        # Shuffle the list of subarrays so it is tested each time in a random order.
        # Although the final ordering is random, always start with a full frame test
        # for convenience of testing.
        for detector in ['MIRIFUSHORT', 'MIRIFULONG']:
            for band in ['SHORT', 'MEDIUM', 'LONG']:
                print( "MRS %s band %s." % (detector, band) )
                test_output_file_name = "%s_%s_%s_%d.fits" % \
                        (test_output_stub, detector, band, ii)
                illumination_map = MiriIlluminationModel(test_input_file_name)
                illumination_map.set_instrument_metadata(detector, band=band)
                illumination_map.set_wcs_metadata(wcsaxes=2, cunit=['pixels','pixels'])
                if VERBOSE > 1:
                    print( illumination_map )
                exposure_map = sca.simulate_pipe(illumination_map, scale=10.0,
                        readout_mode=mode, ngroups=ngroups, nints=1, start_time='NOW',
                        cosmic_ray_mode='SOLAR_MIN', makeplot=PLOTTING,
                        include_pixeldq=INCLUDE_PIXELDQ, qe_adjust=QE_ADJUST,
                        simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                        simulate_read_noise=SIMULATE_READ_NOISE,
                        simulate_ref_pixels=SIMULATE_REF_PIXELS,
                        simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                        simulate_dark_current=SIMULATE_DARK_CURRENT,
                        simulate_flat_field=SIMULATE_FLAT_FIELD,
                        simulate_gain=SIMULATE_GAIN,
                        simulate_nonlinearity=SIMULATE_NONLINEARITY,
                        simulate_drifts=SIMULATE_DRIFTS,
                        verbose=VERBOSE)
                if VERBOSE > 1:
                    print( exposure_map )
                    (exposure_time, duration, start_time, mid_time, end_time) = \
                        exposure_map.get_exposure_times()
                    print("Exposure time %.2f, duration %.2f, start %.2f, end %.2f" % \
                        (exposure_time, duration, start_time, end_time))
                exposure_map.save(test_output_file_name, overwrite=True)
                # Remove the internal reference to the exposure data
                sca.clear_exposure_data()
                del illumination_map, exposure_map
                ii += 1

    if TEST_SUBARRAYS:
        # ------------------------------------------------------------------
        print("\nTesting simulate_pipe with all subarray outputs...\n" + (50 * "*"))
        mode = 'FAST'
        ngroups = 16
        del sca
        sca = SensorChipAssembly3()
        ii = 1
        # Shuffle the list of subarrays so it is tested each time in a random order.
        # Although the final ordering is random, always start with a full frame test
        # for convenience of testing.
        random_list = list(detector_properties['STANDARD_SUBARRAYS'])
        random.shuffle(random_list)
        random_list = ['FULL'] + random_list
        for subarray in random_list:
            print( "Subarray mode FULL->%s." % subarray )
            test_output_file_name = "%s_%s_%d.fits" % \
                    (test_output_stub, subarray, ii)
            illumination_map = MiriIlluminationModel(test_input_file_name)
            illumination_map.set_instrument_metadata('MIRIMAGE')
            illumination_map.set_wcs_metadata(wcsaxes=2, cunit=['pixels','pixels'])
            if VERBOSE > 1:
                print( illumination_map )
            exposure_map = sca.simulate_pipe(illumination_map, scale=10.0,
                    readout_mode=mode, subarray=subarray,
                    ngroups=ngroups, nints=1, start_time='NOW',
                    cosmic_ray_mode='SOLAR_MIN', makeplot=PLOTTING,
                    include_pixeldq=INCLUDE_PIXELDQ,  qe_adjust=QE_ADJUST,
                    simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                    simulate_read_noise=SIMULATE_READ_NOISE,
                    simulate_ref_pixels=SIMULATE_REF_PIXELS,
                    simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                    simulate_dark_current=SIMULATE_DARK_CURRENT,
                    simulate_flat_field=SIMULATE_FLAT_FIELD,
                    simulate_gain=SIMULATE_GAIN,
                    simulate_nonlinearity=SIMULATE_NONLINEARITY,
                    simulate_drifts=SIMULATE_DRIFTS,
                    verbose=VERBOSE)
            if VERBOSE > 1:
                print( exposure_map )
                (exposure_time, duration, start_time, mid_time, end_time) = \
                    exposure_map.get_exposure_times()
                print("Exposure time %.2f, duration %.2f, start %.2f, end %.2f" % \
                    (exposure_time, duration, start_time, end_time))
            exposure_map.save(test_output_file_name, overwrite=True)
            # Remove the internal reference to the exposure data
            sca.clear_exposure_data()
            del illumination_map, exposure_map
            ii += 1
    
        # ------------------------------------------------------------------
        print("\nTesting simulate_pipe with all subarray inputs...\n" + (50 * "*"))
        mode = 'FAST'
        ngroups = 16
        # Shuffle the list of subarrays so it is tested each time in a random order.
        # Although the final ordering is random, always start with a full frame test
        # for convenience of testing.
        random_list = list(detector_properties['STANDARD_SUBARRAYS'])
        random.shuffle(random_list)
        random_list = ['FULL'] + random_list
        for subarray in random_list:
            print( "Subarray mode %s->FULL." % subarray )
            test_output_file_name = "%s_%s_FULL_%d.fits" % \
                    (test_output_stub, subarray, ii)
            subproperties = detector_properties.get('SUBARRAY', subarray)
            if subproperties is not None:
                # Simulate a flat illumination the same size and shape as
                # the subarray.
                subshape = [ subproperties[-2], subproperties[-1] ]
                intensity = 200.0 * np.ones( subshape )
                illumination_map = MiriIlluminationModel(intensity=intensity)
                illumination_map.set_instrument_metadata('MIRIMAGE')
                illumination_map.set_subarray_metadata( subarray )
                illumination_map.set_wcs_metadata(wcsaxes=2, cunit=['pixels','pixels'])
                if VERBOSE > 1:
                    print( illumination_map )
                # Place the subarray onto a FULL frame simulation.
                # Turn off flat-field simulation to prevent the subarray
                # being masked by the flat-field.
                exposure_map = sca.simulate_pipe(illumination_map, scale=10.0,
                        readout_mode=mode, subarray='FULL', ngroups=ngroups,
                        nints=1, cosmic_ray_mode='SOLAR_MIN', 
                        makeplot=PLOTTING, include_pixeldq=INCLUDE_PIXELDQ,
                        qe_adjust=False,
                        simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                        simulate_read_noise=SIMULATE_READ_NOISE,
                        simulate_ref_pixels=SIMULATE_REF_PIXELS,
                        simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                        simulate_dark_current=SIMULATE_DARK_CURRENT,
                        simulate_flat_field=False,
                        simulate_gain=SIMULATE_GAIN,
                        simulate_nonlinearity=SIMULATE_NONLINEARITY,
                        simulate_drifts=SIMULATE_DRIFTS,
                        verbose=VERBOSE)
                if VERBOSE > 1:
                    print( exposure_map )
                    (exposure_time, duration, start_time, mid_time, end_time) = \
                        exposure_map.get_exposure_times()
                    print("Exposure time %.2f, duration %.2f, start %.2f, end %.2f" % \
                        (exposure_time, duration, start_time, end_time))
                exposure_map.save(test_output_file_name, overwrite=True)
                # Remove the internal reference to the exposure data
                sca.clear_exposure_data()
                del illumination_map, exposure_map
                ii += 1

    if TEST_READOUTS:
        # ------------------------------------------------------------------
        # Try all possible readout modes
        mode = 'SLOW'
        inttime = 60.0
        print("\nTesting simulate_files with all readout modes...\n" + (50 * "*"))
        for mode in detector_properties['READOUT_MODE'].keys():    
            # Use a longer integration time for SLOW mode.
            if mode == 'SLOW':
                inttime = 60.0
            else:
                inttime = 10.0
                
            # Case 1: Specify the integration time and let the simulator
            # calculate the number of groups.
            print( "Readout mode %s with time --> ngroups." % mode )
            test_output_file_name = "%s_%s_TIME.fits" % (test_output_stub, mode)
            sca.simulate_files(test_input_file_name, test_output_file_name,
                    'MIRIFULONG', scale=10.0, readout_mode=mode, subarray='FULL',
                    inttime=inttime, nints=1, cosmic_ray_mode='SOLAR_MIN',
                    fileformat='STScI', datashape='hypercube',
                    include_pixeldq=INCLUDE_PIXELDQ, overwrite=True,
                    qe_adjust=QE_ADJUST,
                    simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                    simulate_read_noise=SIMULATE_READ_NOISE,
                    simulate_ref_pixels=SIMULATE_REF_PIXELS,
                    simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                    simulate_dark_current=SIMULATE_DARK_CURRENT,
                    simulate_flat_field=SIMULATE_FLAT_FIELD,
                    simulate_gain=SIMULATE_GAIN,
                    simulate_nonlinearity=SIMULATE_NONLINEARITY,
                    simulate_drifts=SIMULATE_DRIFTS,
#                     fileformat='FITSWriter', datashape='hypercube',
                    makeplot=PLOTTING, verbose=VERBOSE )
            
            # Case 2: Specify the number of groups and let the simulator
            # calculate the integration time.
            print( "Readout mode %s with ngroups --> time." % mode )
            test_output_file_name = "%s_%s_GROUPS.fits" % (test_output_stub, mode)
            sca.simulate_files(test_input_file_name, test_output_file_name,
                    'MIRIFUSHORT', scale=10.0, readout_mode=mode, ngroups=8, nints=2,
                    cosmic_ray_mode='SOLAR_MIN',
                    fileformat='STScI', datashape='hypercube',
                    include_pixeldq=INCLUDE_PIXELDQ, overwrite=True,
                    qe_adjust=QE_ADJUST,
                    simulate_poisson_noise=SIMULATE_POISSON_NOISE,
                    simulate_read_noise=SIMULATE_READ_NOISE,
                    simulate_ref_pixels=SIMULATE_REF_PIXELS,
                    simulate_bad_pixels=SIMULATE_BAD_PIXELS,
                    simulate_dark_current=SIMULATE_DARK_CURRENT,
                    simulate_flat_field=SIMULATE_FLAT_FIELD,
                    simulate_gain=SIMULATE_GAIN,
                    simulate_nonlinearity=SIMULATE_NONLINEARITY,
                    simulate_drifts=SIMULATE_DRIFTS,
#                     fileformat='FITSWriter', datashape='hypercube',
                    makeplot=PLOTTING, verbose=VERBOSE )
    del sca

    print( "Test finished." )
