#!/usr/bin/env python
#
# :History:
#
# 05 Jul 2010: Created
# 16 Jul 2010: Try importing the SensorChipAssembly module in two different
#              ways. Added temperature argument and --silent option.
# 29 Jul 2010: Added cosmic ray mode and help text.
# 30 Jul 2010: Optional arguments converted to parser options.
# 02 Aug 2010: Added subarray argument.
# 06 Aug 2010: Documentation formatting problems corrected.
# 17 Aug 2010: Focal plane module added to parameters list.
# 21 Sep 2010: Added flags to turn Poisson noise and read noise on and off.
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: plot option added and unwanted option strings removed.
#              Another flag to turn amplifier bias, gain and non-linearity
#              on and off. Option of creating a file in FITSWriter or
#              level 1 format.
# 18 Oct 2010: Added --noqe flag.
# 08 Nov 2010: Modified description of subarray mode.
# 18 Nov 2010: Added --nobadpixels option.
# 10 Jan 2011: Detector and amplifier properties classified and looked up
#              by SCA_ID rather than by FPM_ID. fpmodule parameter
#              renamed to scaid.
# 04 Mar 2011: Documentation tweaks to resolve bad formatting.
# 30 Mar 2011: Default subarray mode changed to None, so that a default
#              value may be picked up from the FITS header if not
#              explicitly defined here.
# 20 Apr 2011: Intensity scaling factor added - for debugging.
# 05 Oct 2011: Added --previousfile option so detector persistence
#              effects can be simulated.
# 11 Jan 2012: Renamed symbol to avoid potential name clash:
#              format-->fileformat.
# 13 Nov 2012: Major restructuring of the package folder.
#              Import statements updated.
# 18 Jun 2013: Modified to match the new parameter names used by
#              simulate_sca. 'STSCI' added to file format options.
# 14 Jul 2014: Detector properties are looked up by detector name
#              rather than sca_id. Make the default 'IM'.
# 18 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
#              scaid renamed to detector.
# 02 Mar 2015: Subarray burst mode parameter added.
# 06 Mar 2015: Modified the selection of output data format and shape.
# 14 Aug 2015: Tidied up line breaks.
# 08 Sep 2015: Made compatible with Python 3.
# 25 Sep 2015: Modified to use the SensorChipAssembly simulate() methods,
#              instead of the global functions. Corrected an import problem.
# 04 Dec 2015: Added the ability to redirect output to a Python logger.
# 25 Feb 2016: Cosmic ray simulation suffix parameter added, to choose between
#              different variations of cosmic ray simulation.
# 23 Mar 2016: Added --noflat option.
# 05 May 2016: simulate_amp_effects split into simulate_gain and
#              simulate_nonlinearity
# 06 Jun 2016: Added some missing simulate_ref_pixels and simulate_nonlinearity
#              flags. Added flags to allow the version numbers of imported CDPs
#              to be specified. Improved documentation. Subarray burst mode
#              made True by default.
# 14 Jul 2016: Added simulate_drifts and simulate_latency flags.
# 20 Jan 2017: Replaced "clobber" parameter with "overwrite".
# 07 Feb 2017: Added option to control whether output file contains
#              data quality information.
# 13 Jun 2017: Added cdp_ftp_path parameter.
# 26 Oct 2017: Set cdp_ftp_path to default if not specified.
# 
# @author: Steven Beard (UKATC)

"""

The 'scasim.py' script runs the MIRI SCA simulator using a detector
illumination information contained in a MiriIlluminationModel file
and writes an output file (by default) in level 1b ramp format,
compatible with the JWST DMS data pipeline. Other output data
formats are possible (see the format and datashape parameters).

A MiriIlluminationModel can be read from a FITS file containing
the following IMAGE extensions:

    INTENSITY: describes the intensity of light falling on each
        pixel of the detector (in photons/second).
    WAVELENGTH: describes the wavelength of light falling on
        each pixel of the detector (in microns).
    DIRECTION: describes the direction the light beam falling
        on each pixel of the detector (in degrees).
        
The WAVELENGTH and DIRECTION extensions are optional. Illumination
at multiple wavelengths can be described using 3-dimensional
INTENSITY or WAVELENGTH arrays.

JWST level 1b ramp data can be stored in a FITS file containing
the following extensions:

    SCI: describes the detector readout (in DN)
    PIXELDQ: describes the quality of the detector readout (optional)
    REFOUT: contains the detector reference output data (in DN)

The compulsory scasim.py command arguments are:

    inputfile
        The path+name of the detector illumination file.
        This must be a FITS file in a format compatible with
        the MiriIlluminationModel data model.
    outputfile
        The path+name of the output FITS file to be created.
        By default, this is a level 1b file compatible with
        the MiriRampModel data model and the JWST DMS pipeline.
        Legacy ground-based data formats may be created by
        specifying the --format and --datashape parameters.

The following optional parameters may be provided by keyword:

    --detector
        The identity of the detector module to be simulated.
        Options for MIRI are 'MIRIMAGE', 'MIRIFULONG' and 'MIRIFUSHORT'.
        (These correspond to the imager and the long and short wavelength
        arms of the MRS respectively.) If not specified, the SCA ID will
        be derived from the FITS header of the detector illumination file,
        or if that isn't possible a default of 'MIRIMAGE' (imager) will be
        assumed.
    --fringemap
        The path+name of a FITS file containing a fringe map to be used
        to modulate the detector illumination. If not specified, no
        fringe map will be applied.
    --rdmode
        The detector readout mode ('SLOW', 'SLOWGRPAVG', 'SLOWINTAVG',
        'SLOWGRPGAP', 'FAST', 'FASTGRPAVG', 'FASTINTAVG' or
        'FASTGRPGAP'). If not specified, the value will be taken from
        the FITS header of the input file, or failing that the default
        contained in the detector properties file ('FAST') will be used.
    --subarray
        The output subarray mode, if needed ('FULL', 'MASK1550',
        'MASK1140', 'MASK1065', 'MASKLYOT', 'BRIGHTSKY',  'SUB256',
        'SUB128', 'SUB64', 'SLITLESSPRISM', etc.). If not specified, it
        is assumed to be the same as the input data. If neither are
        specified (or the input data specifies a non-standard subarray),
        'FULL' will be assumed.
    --noburstmode
        Include this parameter if subarray data are not to be simulated
        in burst mode.
    --inttime
        The integration time in seconds. If not specified, the
        integration time will be derived from the read mode and
        --ngroups. Use shorter times when the test data is small and
        when the detector readout mode is FAST.
        *NOTE: Any requested integration time will be rounded up to
        the time resulting from the nearest whole number of groups.*
        *NOTE: It is not sensible to specify both --inttime and
        --ngroups. If both parameters are specified, --ngroups will
        set the integration time and the --inttime value will be
        ignored. If neither parameter is specified, the default values
        associated with the detector readout mode are used.*
    --ngroups
        The number of groups making up each integration. If not
        specified, ngroups will be derived from the detector readout
        mode and integration time. If both parameters are specified,
        groups takes priority over inttime.
    --nints
        The number of integrations. The total exposure time will be
        inttime x nints. If not specified, the default number of
        integrations associated with the detector readout mode will be
        used.
    --wait
        The time that elapses in between exposures. 
    --temperature
        The detector temperature in K. If not specified it will default
        to the expected target temperature contained in the detector
        properties file.
    --crmode
        The cosmic ray environment ('NONE', 'SOLAR_MIN', 'SOLAR_MAX'
        or 'SOLAR_FLARE'). If not specified, the value will be taken
        from the FITS header of the input file, or failing that the
        default contained in the cosmic ray properties file
        ('SOLAR_MIN') will be used.
        *NOTE: 'SOLAR_FLARE' simulations can take a long time because
        of the high number of cosmic ray events involved. Use very
        short exposures.*
        An optional suffix can be added to this parameter (separated
        by a '+' to select a particular variation of the cosmic ray
        simulation library. For example, 'SOLAR_MIN+IPC' selects the
        variant of the solar minimum library which includes 'IPC'
        charge crosstalk effects. 
    --format
        The file format required ('STSCI' or 'FITSWriter').
        In 'FITSWriter' format the data array is written to the primary
        HDU. In 'STSCI' format the data array is written to an 'SCI'
        extension HDU in a format compatible with the STScI level 1b
        ramp data model. If format is set to 'STSCI-DQ', the output
        file will not contain any data quality information.
        If not specified, the default format is 'STSCI'.
    --datashape
        The data shape required ('hypercube', 'cube' or 'slope'). Controls
        whether each integration is saved as a separate cube within a
        hypercube or whether integrations are concatenated to make one
        big cube. The 'slope' option applies a straight line fit to the
        data to generate one 2-D image per integration. The default for
        'STSCI' data format is 'hypercube', which matches level 1b ramp
        format. The default for other data formats is 'cube', which
        generates data in a similar format to ground-based test data. 
    --scale
        Scale factor to be applied to the intensity data (for debugging).
        If not specified, a scale of 1.0 will be used. This parameter
        can be used to make faulty input data usable.
    --cdp_ftp_path
        A list of folders (or folders) on the SFTP host to be searched
        for CDP files, consisting of a list of folder names separated by
        a ":" delimiter. Examples: 'CDP', 'CDPSIM', 'CDPSIM:CDP:CDPTMP'.
        Do not change this parameter unless you know which CDP files
        will be imported. Defaults to the standard simulator CDP
        repository at Leuven.
    --cdprelease
        The CDP release from which CDP files are to be imported.
        Do not change this parameter unless you know which CDP files
        will be imported. Defaults to the latest release.
    --previousfile
        The path+name of a detector illumination file containing a
        previous exposure. If defined, this file is processed first and
        will cause detector persistence effects in the main observation.
        By default this file name is undefined and no persistence effects
        are simulated. If specified, this must be a FITS file compatible
        with the MiriIlluminationModel data model.
    
The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output.
    --debug or -d:
        Generate maximum output.
    --plot:
        Generate plots.
    --overwrite or -o:
        Overwrite any existing FITS output file.
    --noqe:
        Simulation without quantum efficiency adjustment (input in electrons/s)
    --nopoisson:
        Simulation without Poisson noise.
    --noreadnoise:
        Simulation without read noise.
    --norefpixels:
        Simulation without reference pixels.
    --nobadpixels:
        Simulation without bad pixels.
    --nodark:
        Simulation without dark current or hot pixels.
    --noflat:
        Simulation without flat-field adjustment.
    --nogain:
        Simulation without amplifier bias and gain.
    --nolinearity:
        Simulation without non-linearity effects.
    --nodrifts:
        Simulation without detector drifts effects.
    --nolatency:
        Simulation without detector latency effects.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default logger

import optparse
import os, sys, time

from miri.simulators.scasim.exposure_data import get_file_header
from miri.simulators.scasim.sensor_chip_assembly import SensorChipAssembly1
from miri.simulators.scasim.detector import SIM_CDP_FTP_PATH

if __name__ == "__main__":
    # parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile outputfile"
    usage += "\n\t[--detector] [--fringemap] [--rdmode] [--subarray] "
    usage += "[--inttime] [--ngroups]"
    usage += "\n\t[--nints] [--temperature] [--crmode] [--format] [--datashape]"
    usage += "\n\t[--scale] [--cdp_ftp_path] [--cdprelease] [--previousfile]"
    usage += "[--nopoisson] [--noreadnoise] [--norefpixels] [--nobadpixels]"
    usage += "\n\t[--nodark] [--noflat] [--nogain] [--nolinearity] "
    usage += "[--nodrifts] [--nolatency]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--detector", dest="detector", type="string",
                     default='', help="Sensor Chip Assembly ID"
                     )
    parser.add_option("", "--fringemap", dest="fringemap", type="string",
                     default='', help="Path+name of fringe map FITS file"
                     )
    parser.add_option("", "--rdmode", dest="rdmode", type="string",
                     default=None, help="Detector readout mode"
                     )
    parser.add_option("", "--subarray", dest="subarray", type="string",
                     default=None, help="Detector subarray mode for output"
                     )
    parser.add_option("", "--noburstmode", dest="noburstmode", default=False,
                      action="store_true", help="Turn off subarray burst mode"
                     )
    parser.add_option("", "--inttime", dest="inttime", type="float",
                     default=None, help="Integration time in seconds"
                     )
    parser.add_option("", "--ngroups", dest="ngroups", type="int",
                     default=None, help="Number of groups"
                     )
    parser.add_option("", "--nints", dest="nints", type="int",
                     default=None, help="Number of integrations"
                     )
    parser.add_option("", "--wait", dest="wait", type="float",
                     default=None, help="Elapsed time in seconds"
                     )
    parser.add_option("", "--temperature", dest="temp", type="float",
                     default=None, help="Detector temperature in K"
                     )
    parser.add_option("", "--crmode", dest="crmode", type="string",
                     default=None, help="Cosmic ray environment mode"
                     )
    parser.add_option("", "--format", dest="fileformat", type="string",
                     default='STSCI', help="Output file format required"
                     )
    parser.add_option("", "--datashape", dest="datashape", type="string",
                     default=None, help="Data shape required"
                     )
    parser.add_option("", "--scale", dest="scale", type="float",
                     default=None, help="Intensity scale factor"
                     )
    parser.add_option("", "--cdp_ftp_path", dest="cdp_ftp_path", type="string",
                     default='', help="Search path for imported CDPs"
                     )
    parser.add_option("", "--cdprelease", dest="cdprelease", type="string",
                     default='', help="Release number for imported CDPs"
                     )
    parser.add_option("", "--previousfile", dest="previous", type="string",
                     default='',
                     help="Path+name of FITS file for previous exposure"
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
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot calibration and simulated data"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite output file if it exists"
                     )

    parser.add_option("-q", "--noqe", dest="noqe",
                      action="store_true", help="Turn off QE adjustment"
                     )
    parser.add_option("-a", "--nopoisson", dest="nopoisson",
                      action="store_true", help="Turn off Poisson noise"
                     )
    parser.add_option("-r", "--noreadnoise", dest="noreadnoise",
                      action="store_true", help="Turn off read noise"
                     )
    parser.add_option("-x", "--norefpixels", dest="norefpixels",
                      action="store_true", help="Turn off reference pixels"
                     )
    parser.add_option("-b", "--nobadpixels", dest="nobadpixels",
                      action="store_true", help="Turn off bad pixels"
                     )
    parser.add_option("-k", "--nodark", dest="nodark",
                      action="store_true", help="Turn off dark current"
                     )
    parser.add_option("-f", "--noflat", dest="noflat",
                      action="store_true", help="Turn off pixel flat-field"
                     )
    parser.add_option("-g", "--nogain", dest="nogain",
                      action="store_true",
                      help="Turn off bias and gain"
                     )
    parser.add_option("-l", "--nolinearity", dest="nolinearity",
                      action="store_true",
                      help="Turn off non-linearity effects"
                     )
    parser.add_option("-t", "--nodrifts", dest="nodrifts",
                      action="store_true",
                      help="Turn off detector drifts effects"
                     )
    parser.add_option("-y", "--nolatency", dest="nolatency",
                      action="store_true",
                      help="Turn off detector latency effects"
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
    detector = options.detector
    fringemap = options.fringemap
    rdmode = options.rdmode
    subarray = options.subarray
    noburstmode = options.noburstmode
    if noburstmode:
        burstmode = False
    else:
        burstmode = True
    time = options.inttime
    ngroups = options.ngroups
    nints = options.nints
    wait_time = options.wait
    temp = options.temp
    crmode = options.crmode
    fileformat = options.fileformat
    datashape = options.datashape
    # Determine whether data quality is to be included.
    include_dq = True
    if fileformat == 'STSCI-DQ' or fileformat == 'STScI-DQ':
        fileformat = 'STSCI'
        include_dq = False
    # The default shape depends on the file format.
    if datashape is None:
        if fileformat == 'STSCI' or fileformat == 'STScI':
            datashape = 'hypercube'
        else:
            datashape = 'cube'
    scale = options.scale
    if scale is None:
        scale = 1.0
    cdp_ftp_path = options.cdp_ftp_path
    if not cdp_ftp_path:
        cdp_ftp_path = SIM_CDP_FTP_PATH
    cdprelease = options.cdprelease
    previous = options.previous

    # Boolean flags.
    overwrite = options.overwrite
    verb = options.verb
    debug = options.debug
    silent = options.silent
    makeplot = options.makeplot
    noqe = options.noqe
    nopoisson = options.nopoisson
    noreadnoise = options.noreadnoise
    nobadpixels = options.nobadpixels
    norefpixels = options.norefpixels
    nodark = options.nodark
    noflat = options.noflat
    nogain = options.nogain
    nolinearity = options.nolinearity
    nodrifts = options.nodrifts
    nolatency = options.nolatency
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb: verbose = 2
    if debug: verbose = 4
    if silent: verbose = 0
    
    if noqe:
        qe_adjust = False
    else:
        qe_adjust = True
    if nopoisson:
        simulate_poisson_noise = False
    else:
        simulate_poisson_noise = True
    if noreadnoise:
        simulate_read_noise = False
    else:
        simulate_read_noise = True
    if norefpixels:
        simulate_ref_pixels = False
    else:
        simulate_ref_pixels = True
    if nobadpixels:
        simulate_bad_pixels = False
    else:
        simulate_bad_pixels = True
    if nodark:
        simulate_dark_current = False
    else:
        simulate_dark_current = True
    if noflat:
        simulate_flat_field = False
    else:
        simulate_flat_field = True
    if nogain:
        simulate_gain = False
    else:
        simulate_gain = True
    if nolinearity:
        simulate_nonlinearity = False
    else:
        simulate_nonlinearity = True
    if nodrifts:
        simulate_drifts = False
    else:
        simulate_drifts = True
    if nolatency:
        simulate_latency = False
    else:
        simulate_latency = True

    # If the SCA_ID has not been explicitly given, look it up in the
    # FITS header of the input file. Failing that, set it to a default.
    if not detector:
        input_header = get_file_header(inputfile)
        if input_header is not None and 'DETECTOR' in input_header:
            detector = input_header['DETECTOR']
    if not detector:
        detector = 'MIRIMAGE'
        if verbose > 0:
            LOGGER.info( "Assuming default detector ID: %s" % str(detector) )

    sca = SensorChipAssembly1(logger=LOGGER)
    if previous:
        # Previous exposure specified. Run a simulation on a list of files.
        infilelist = (previous, inputfile)
        outfilelist = ('', outputfile)
        sca.simulate_files(infilelist, outfilelist, detector, scale=scale,
                     fringemap=fringemap,
                     readout_mode=rdmode,
                     subarray=subarray, burst_mode=burstmode,
                     inttime=time, ngroups=ngroups, nints=nints,
                     wait_time=wait_time, temperature=temp,
                     cosmic_ray_mode=crmode,
                     fileformat=fileformat, datashape=datashape,
                     overwrite=overwrite, include_pixeldq=include_dq,
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
                     readnoise_version=cdprelease,
                     bad_pixels_version=cdprelease,
                     flat_field_version=cdprelease, 
                     gain_version=cdprelease,
                     makeplot=makeplot, verbose=verbose)
    else:
        # No previous exposure specified. Run a basic simulation.
        sca.simulate_files(inputfile, outputfile, detector, scale=scale,
                     fringemap=fringemap,
                     readout_mode=rdmode,
                     subarray=subarray, burst_mode=burstmode,
                     inttime=time, ngroups=ngroups, nints=nints,
                     wait_time=wait_time, temperature=temp,
                     cosmic_ray_mode=crmode,
                     fileformat=fileformat, datashape=datashape,
                     overwrite=overwrite, include_pixeldq=include_dq,
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
                     readnoise_version=cdprelease,
                     bad_pixels_version=cdprelease,
                     flat_field_version=cdprelease, 
                     gain_version=cdprelease,
                     makeplot=makeplot, verbose=verbose)
    del sca
    if verbose > 0:
        LOGGER.info( "Simulation finished." )
