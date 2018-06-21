#!/usr/bin/env python
#
# Script 'make_qe_fits' is based on the miri.tools 'make_filters_fits' script.
# It creates a FITS file describing a QuantumEfficiency filter from ASCII data
# and command-line arguments.
#
# :History:
# 
# 23 Aug 2010: Created
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: Unwanted option strings removed.
# 11 Jan 2011: Detector identified by SCA ID rather than by FPM ID.
# 15 Mar 2011: Moved from scasim to miri.tools.
# 29 Mar 2012: QuantumEfficiency class now imported from miri.tools.filters
#              module. Unnecessary imports removed.
# 13 Apr 2012: Major changes to implement the Filter class as a JwstDataProduct.
# 13 Nov 2012: Major restructuring of the package folder. Import statements
#              updated.
# 11 Jul 2014: Switch over to using the new MiriQuantumEfficiency class,
#              based on the jwst_lib data models.
# 08 Sep 2015: Made compatible with Python 3.
#
#@author: Steven Beard (UKATC)

"""

The 'make_qe_fits' script is based on the miri.tools 'make_filters_fits'
script. It creates a FITS file describing a MiriQuantumEfficiency
filter from ASCII data and command-line arguments. This utility is
useful for converting flight model test data for the MIRI focal plane
modules into FITS files. The compulsory command arguments are:

    filename
        The path+name of the file to be read. The FITS table file
        will have the same name but with a '.fits' extension.
    qeID
        An identifier for this quantum efficiency measurement
        (equivalent to the filterID used by miri.tools). For example
        'QE493' for the QE measurement of sensor chip assembly 493.
    instrument
        The name of the instrument (normally 'MIRI').
    detector
        The sensor chip assembly ID to which this quantum efficiency
        measurement refers (e.g. 493, 494 or 495).
    wavunit
        Wavelength unit (e.g. 'microns').
    temperature
        Temperature in K at which the QE measurement was made.

The following optional parameters may be provided by keyword:

    --version
        An optional string describing the version of the quantum
        efficiency measurement (recommended for change control of
        the measurements).
    --comment
        An optional comment describing the measurement.

The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output.
    --debug or -d:
        Generate maximum output.
    --plot or -p:
        Plot the quantum efficiency measurement.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""
# This module is now converted to Python 3.


import optparse
import os, sys, time

from miri.datamodels.miri_filters import MiriQuantumEfficiency, \
    ascii_to_filter

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] fname qeID instrument detector wavunit temperature"
    usage += "\n\t[--version] [--comment]\n"
    usage += "Creates a FITS file containing quantum efficiency info from "
    usage += "the provided ASCII file and arguments"
    parser = optparse.OptionParser(usage)

    # Optional arguments (long option strings only).
    parser.add_option("", "--version", dest="version", type="string",
                     default='', help="Version number of measurement"
                     )
    parser.add_option("", "--comment", dest="comment", type="string",
                     default='', help="Optional comment"
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
                      help="Plot filter transmission"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )

    (options, args) = parser.parse_args()

    # Compulsory arguments
    try:
        # TODO: Parse these command arguments by name rather than by position.
        fname = args[0].rstrip()
        author = 'MIRI Software Team'
        qeID = args[1]
        instrument = args[2]
        obsmode = 'any'     # FPA simulations don't depend on instrument mode.
        detector = int(args[3])
        wavunit = args[4]
        temperature = args[5]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    version = options.version
    comment = options.comment

    # Boolean flags.
    verb = options.verb
    debug = options.debug
    silent = options.silent
    makeplot = options.makeplot
    overwrite = options.overwrite

    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb: verbose = 2
    if debug: verbose = 4
    if silent: verbose = 0

# Read data from ASCII file. The name must end in ".txt"
    if fname[-3:] != "txt":
        fname = fname + ".txt"
    if verbose > 0:
        print( "Reading %s\n" % fname )
    qe_filt = ascii_to_filter(fname, filter_name='ANY', detector=detector,
                              temperature=temperature)
    if author or version:
        qe_filt.set_housekeeping_metadata('MIRI EC', author=author,
                                          version=version)
    qe_filt.set_instrument_metadata(detector, modelnam='FM',
                                filt='ANY', channel='', band='ANY',
                                ccc_pos='OPEN', deck_temperature=None,
                                detector_temperature=temperature)
    if verbose > 1:
        print( qe_filt, "\n" )

    # Plot filter transmission if required
    if makeplot:
        qe_filt.plot()

    # Save data in FITS file
    basename, extension = os.path.splitext(fname)
    newfname = basename + '.fits'
    qe_filt.save(newfname, overwrite=overwrite)
    if verbose > 0:
        print( "Data saved to %s\n" % newfname )
