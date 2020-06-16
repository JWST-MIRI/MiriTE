#!/bin/sh
#
# An example script showing how to use make_qe_fits.py to read an ASCII
# file containing a quantum efficiency measurement and generate the
# corresponding fits file.
#
# This script assumes that make_qe_fits.py has been installed into the
# PATH by running "setup.py install"
#

# This command converts the measurements described in the ASCII file
# example_qe.txt into FITS format in example_qe.fits.
# Note that additional information for the FITS header can be added to
# the command line.
make_qe_fits.py --clobber ../data/example_qe.txt "QE_Example" "MIRI" "495" "microns" 273.0 --version "0.2" --comment "An example QE measurement"
