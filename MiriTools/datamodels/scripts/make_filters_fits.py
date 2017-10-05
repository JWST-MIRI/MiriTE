#!/usr/bin/env python
#
# Script 'make_filters_fits' creates a FITS file describing a filter from
# ASCII data and command-line arguments.
#
# :History:
# 
# 23 Aug 2010: Created
# 29 Mar 2012: History and docstring added and unwanted imports removed.
# 02 Apr 2012: Corrected documentation.
# 13 Apr 2012: Major changes to implement the Filter class as a JwstDataProduct.
# 13 Nov 2012: Major restructuring of the teams/miri folder. Import statements
#              updated.
# 01 Jun 2015: Corrected to use ".save" instead of ".tofits".
# 08 Sep 2015: Made compatible with Python 3. Corrected import errors.
#
#@author: Julien Morin (DIAS), Steven Beard (UKATC)

"""

Script `make_filters_fits` creates a FITStable file describing a filter
from ASCII data and command-line arguments. This utility is useful
for converting flight model test data for the MIRI focal plane modules
into FITS files. The following command arguments are defined by position,
and all must be supplied:

    filename[0]
        The path+name of the file to be read. The FITS table file
        will have the same name but with a '.fits' extension.
    filter[1]
        An identifier for the filter (e.g. 'F2100W'.
    filtertype[2]
        The filter type (e.g. 'BP' for a band pass filter).
    instrument[3]
        The name of the instrument (normally 'MIRI').
    obsmode[4]
        The observation mode (e.g. 'IMG').
    wavecentre[5]
        The central wavelength for a band pass filter.
    fwhm[6]
        The FWHM of the wavelength for a band pass filter.
    wavunit[7]
        Wavelength unit (e.g. 'microns').
    temperature[8]
        Temperature in K at which the QE measurement was made.
    version[9]
        The version code for this filter measurement (e.g. 'V1.0')

The command also takes the following options:

    --verbose or -v:
        Generate more output.
    --plot or -p:
        Plot the filter measurement.
    --overwrite or -o:
        Overwrite any existing FITS file.

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import optparse
import os, sys, time

from miri.datamodels.miri_filters import MiriFilter, MiriBandPassFilter

if __name__ == "__main__":
    # Parse arguments
    help_text = __doc__
    usage = "%prog [opt] fname filter instrument wavecentre fwhm wavunit\n"
    usage += "       temperature version\n"
    usage += "Creates a FITS file containing filter info from the provided"
    usage += "ASCII file and arguments"
    parser = optparse.OptionParser(usage)
    parser.add_option("-v", "--verbose", dest="verb", action="store_true",
                      help="Verbose mode"
                     )
    parser.add_option("-p", "--plot", dest="makeplot", action="store_true",
                      help="Plot filter transmission"
                     )
    parser.add_option("-o", "--overwrite", dest="overwrite", action="store_true",
                      help="Overwrite the FITS file if it already exists"
                     )

    (options, args) = parser.parse_args()

    try:
        # TODO: Parse these command arguments by name rather than by position.
        fname = args[0].rstrip()
        author = 'MIRI Software Team'
        filter = args[1]
        filtertype = args[2]
        instrument = args[3]
        obsmode = args[4]
        wavecentre = float(args[5])
        fwhm = float(args[6])
        wavunit = args[7]
        temperature = args[8]
        version = args[9]
    except IndexError:
        print( help_text )
        time.sleep(1) # Ensure help text appears before error messages.
        parser.error("Not enough arguments provided")
        sys.exit(1)

    verb = options.verb
    makeplot = options.makeplot
    overwrite = options.overwrite

# Read data from ASCII file. The name must end in ".txt"
    if fname[-3:] != "txt":
        fname = fname + ".txt"
    if verb:
        print( "Reading %s\n" % fname )
    if filtertype == 'BP':
        filt = MiriBandPassFilter(fname, filtername=filter, instrument=instrument,
                              obsmode=obsmode, wavecentre=wavecentre, fwhm=fwhm,
                              wavunit=wavunit,
                              temperature=temperature, version=version,
                              removeneg=True, removeout=5.0)
    else:
        filt = MiriFilter(fname, filtername=filter, filtertype=filtertype,
                      instrument=instrument, obsmode=obsmode,
                      wavunit=wavunit, temperature=temperature, version=version,
                      removeneg=True)
        
    if author or version:
        filt.set_housekeeping_metadata('MIRI EC', author=author,
                                       version=version)
    if verb:
        print( filt, "\n" )

# Plot filter transmission if required
    if makeplot:
        filt.plotfig()

# Save data in FITS file
    basename, extension = os.path.splitext(fname)
    newfname = basename + '.fits'
    filt.save(newfname, overwrite=overwrite)
    if verb:
        print( "Data saved to %s\n" % newfname )
