#!/bin/sh
#
# Run make_qe_fits.py on the quantum efficiency measurements contained in
# ASCII files in the miri/simulators/data/detector folder to generate the
# corresponding fits files.
#
# This script assumes that make_qe_fits.py has been installed
# from miri.datamodels into the PATH by running "setup.py install"

# Convert the QE measurements for each of the MIRI focal plane modules.
#
make_qe_fits.py --clobber ../data/detector/qe_measurement493.txt "QE493" "MIRI" "493" "microns" 293.0 --version "Initial X1" --comment "QE prediction for FM FPM S/N 106"
make_qe_fits.py --clobber ../data/detector/qe_measurement494.txt "QE494" "MIRI" "494" "microns" 293.0 --version "Initial X1" --comment "QE prediction for FM FPM S/N 104"
make_qe_fits.py --clobber ../data/detector/qe_measurement495.txt "QE495" "MIRI" "495" "microns" 293.0 --version "Initial X1" --comment "QE prediction for FM FPM S/N 105"
