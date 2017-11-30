#!/usr/bin/env python
#
# Script 'make_measurement_fits' reads one or more variable measurements
# from an ASCII format file and writes them to a FITS table file, together
# with supplementary information supplied as command line arguments..
#
# :History:
# 
# 23 Aug 2010: Created
# 27 Sep 2010: Python environment for windows verified.
# 12 Oct 2010: Unwanted option strings removed.
# 15 Nov 2010: Moved to miri.tools.
# 28 Mar 2012: BUG DISCOVERED! matplotlib import crashes Python 2.7 when
#              combined with the detector_properties import. Reason unknown.
#              Worked around by commenting out the offending code.
# 02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
#              under Windows until the problem has been solved.
# 04 Apr 2012: Improvements suggested by pylint.
# 14 May 2012: Major changes to implement the MeasuredVariable class as a
#              JwstDataProduct.
# 13 Nov 2012: Major restructuring of the package folder. Import statements
#              updated.
# 11 Jul 2014: Modified to use the new MiriMeasurement class instead of
#              the old MeasuredVariable class.
# 01 Jun 2015: Corrected to use ".save" instead of ".tofits".
# 08 Sep 2015: Made compatible with Python 3.
#
#@author: Steven Beard (UKATC)

"""

The 'make_measurements_file' script reads one or more variable
measurements from an ASCII format file and writes them to a FITS table
file, together with supplementary information supplied as command line
arguments. This utility is useful for converting flight model test data
for the MIRI focal plane modules into FITS files. The compulsory command
arguments are:

    inputfile
        The path+name of the file to be read. The FITS table file
        will have the same name but with a '.fits' extension.
    name
        The name of the quantity being measured (e.g. "Read noise").
    unit
        The units of the quantity being measured (e.g. "Electrons").
    vname
        The name of the variable against which the quantity is being
        measured (e.g. "Temperature").
    vunit
        The units of the variable against which the quantity is being
        measured (e.g. "K").

The following optional parameters may be provided by keyword:

    --suffix
        A string used, with column number appended, to distinguish
        between the data columns. The default value is '_ch' which will
        cause the strings '_ch1', '_ch2', '_ch3', etc... to be appended
        to the column names. If there is only one column, or this
        parameter is set to a null string, no suffix will be added.
    --interptype
        The type of interpolation to be used between the data. Possible
        values are:
        
        * 'LINEAR' - Apply linear interpolation
        * 'LOGLIN' - Interpolate the logarithm of the measurements
                     against the unchanged adjustable parameter.
        * 'LINLOG' - Interpolate the measurements linearly against the
                     logarithm of the adjustable parameter.
        * 'LOGLOG' - Interpolate the logarithm of the measurements
                     against the logarithm of the adjustable parameter.
        
        The default value is 'LINEAR'.
    --version
        An optional string describing the version of the variable
        measurement (recommended for change control of the measurements).
    --comment
        An optional comment describing the measurement.

The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.
    --plot:
        Generate plots.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time
import numpy as np

# Import the matlib if available and flag if this has been successful.
# It is only used for plotting and is therefore optional.
try:
    import matplotlib.pyplot as plt
    _PLOT_AVAILABLE = True
except ImportError:
    _PLOT_AVAILABLE = False

from miri.datamodels.miri_measurement import MiriMeasurement, \
    ascii_to_measurement

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] filename name unit pname punit "
    usage += "\n\t [--suffix] [--interptype] "
    usage += "[--version] [--comment]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments (long option strings only).
    parser.add_option("", "--suffix", dest="suffix", type="string",
                     default='_ch', help="Column name suffix"
                     )
    parser.add_option("", "--interptype", dest="interptype", type="string",
                     default='LINEAR', help="Interpolation type"
                     )
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
                      help="Plot the measurements"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite existing file"
                     )

    (options, args) = parser.parse_args()
    
    # Compulsory arguments
    try:
        filename = args[0]
        author = 'MIRI Software Team'
        name = args[1]
        unit = args[2]
        vname = args[3]
        vunit = args[4]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    suffix = options.suffix
    interptype = options.interptype
    version = options.version
    comment = options.comment

    # Boolean flags.
    overwrite = options.overwrite
    verb = options.verb
    debug = options.debug
    silent = options.silent
    if _PLOT_AVAILABLE:
        makeplot = options.makeplot
    elif options.makeplot:
        print( "Sorry, the matplotlib plotting library is not available." )
        makeplot = False
    else:
        makeplot = False
    
    # Set the verbosity level according to the --verbose and --silent
    # options. (Note that --debug wins over --verbose and --silent wins
    # over all the other options if they are provided together.)
    verbose = 1
    if verb:
        verbose = 2
    if debug: 
        verbose = 4
    if silent:
        verbose = 0

    # Read the given file and generate a list of MeasuredVariable objects.
    mv = ascii_to_measurement(filename, title=comment, name=name, unit=unit,
                              vname=vname, vunit=vunit, interptype=interptype)
    if author or version:
        mv.set_housekeeping_metadata("MIRI EC", author=author, version=version)
    ncols = len(mv.values)
    if verbose > 0:
        print( ncols, "columns read." )

        if verbose > 1:
            print( mv )
            
        if makeplot:
            mv.plot()

    # Construct the name of the output file by taking the input file
    # name and replacing the file extension with ".fits".
    basename, extension = os.path.splitext(filename)
    outputfile = basename + '.fits'
    mv.save(outputfile, overwrite=overwrite)
    if verbose > 0:
        print( "Measurement data saved to %s" % outputfile )
