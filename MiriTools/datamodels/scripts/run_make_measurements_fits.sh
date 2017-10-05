#!/bin/sh
#
# An example script showing how to use make_measurements_fits.py to read an
# ASCII file of measurements and generate the corresponding fits file.
#
# This script assumes that make_measurements_fits.py has been installed
# into the PATH by running "setup.py install"
#

# This command converts the measurements described in the ASCII file
# example_measurement.txt into FITS format in example_measurement.fits.
# Note that additional information for the FITS header can be added to
# the command line.
#
make_measurements_fits.py --clobber ../data/example_measurement.txt "Read noise" "electrons" "Temperature" "K" --suffix "_ch" --version "0.2" --comment "Some example read noise measurements"
