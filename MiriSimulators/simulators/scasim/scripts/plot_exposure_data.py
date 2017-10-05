#!/usr/bin/env python
#
# Script `plot_exposure_data` reads a FITS file containing exposure data
# and displays its contents.
#
# :History:
# 
# 08 Sep 2010: Created
# 27 Sep 2010: Python environment for windows verified.
# 23 Mar 2011: ExposureData.plot() changed to ExposureData.plotfig()
# 13 Nov 2012: Major restructuring of the teams/miri folder.
#              Import statements updated.
# 03 Mar 2015: Simulator simplified by removing the ability to save data
#              to legacy file formats.
# 27 May 2015: Replaced pyfits with astropy.io.fits
# 08 Sep 2015: Made compatible with Python 3.
#
#@author: Steven Beard (UKATC)

"""

The plot_exposure_data script reads a FITS file containing exposure
data and displays it it various ways. The compulsory parameters are:

    inputfile

The following optional parameters may be provided by keyword:

    --plottype
        The type of plot to be generated:
        
        * IMAGE - Plot each integration and group as a set of images.
        * RAMP - Plot a ramp at a particular row and column.
        * AVERAGED - Plot an averaged ramp at a particular row or column.

    --row
        In RAMP mode, the row at which the ramp is to be displayed.
    --column
        In RAMP mode, the row at which the ramp is to be displayed.
    --description
        An optional string to be added as a plot title. If not
        specified, the input file name is used.

The command also takes the following options:

    --silent or -s:
        Generate no output.
    --verbose or -v:
        Generate more output and some data plots.
    --debug or -d:
        Generate maximum output.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time
import numpy as np
import astropy.io.fits as pyfits

# MIRI exposure data model.
#from miri.datamodels.sim import MiriExposureModel
from miri.datamodels.miri_measured_model import MiriRampModel

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] inputfile"
    usage += "\n\t[--plottype] [--row] [--column] [--description]"
    parser = optparse.OptionParser(usage)
    
    # Optional arguments
    
    parser.add_option("-p", "--plottype", dest="plottype", type="string",
                     default='IMAGE', help="Plot type"
                     )
    parser.add_option("-r", "--row", dest="row", type="int",
                     default=1, help="Row number at which to plot ramp"
                     )
    parser.add_option("-c", "--column", dest="column", type="int",
                     default=1, help="Column number at which to plot ramp"
                     )
    parser.add_option("-e", "--description", dest="description", type="string",
                     default='', help="Plot description"
                     )
     
    # Boolean options.
    parser.add_option("-d", "--debug", dest="debug", action="store_true",
                      help="Debugging mode"
                     )
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-s", "--silent", dest="silent", action="store_true",
                      help="Silent mode"
                     )

    (options, args) = parser.parse_args()
    
    # Compulsory arguments
    try:
        inputfile = args[0]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)
        
    # Optional arguments.
    plottype = options.plottype
    # The row and column are only relevant for a RAMP plot.
    if plottype == 'RAMP' or plottype == 'AVERAGED':
        row = int(options.row)
        column = int(options.column)
            
    # The input file name is used as a description if one is not provided.
    description = options.description
    if not description:
        description = inputfile

    # Boolean flags.
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

    # Create a new exposure data object from the file.
    exposure = MiriRampModel(inputfile)
    if verbose > 1:
        print( exposure )
        exposure.statistics()

    # Plot the exposure data.
    if plottype == 'RAMP':
        exposure.plot_ramp(row, column, description=description)
    elif plottype == 'AVERAGED':
        exposure.plot_ramp(row, column, averaged=True, description=description)
    else:
        exposure.plot(description=description)
        
    if verbose > 0:
        print( "Plot completed." )
