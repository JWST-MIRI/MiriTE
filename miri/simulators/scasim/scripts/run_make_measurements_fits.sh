#!/bin/sh
#
# Run miri/datamodels/make_measurements_fits.py on the ASCII files contained in
# miri/scasim/data/amplifiers and miri/scasim/data/detector to
# generated the corresponding fits files.
#
# This script assumes that make_measurements_fits.py has been installed
# into the PATH by running "setup.py install"
#

# First convert the ASCII read noise measurements to FITS. There is
# one ASCII file containing read noise measurements for the 5 amplifiers
# for each of the MIRI focal plane modules.
#
make_measurements_fits.py --overwrite ../../data/amplifiers/read_noiseIM.txt "Read noise" "electrons" "Temperature" "K" --suffix "_ch" --version "Initial X1" --comment "Read noise measurements from FPM S/N 106 end item data package"
make_measurements_fits.py --overwrite ../../data/amplifiers/read_noiseLW.txt "Read noise" "electrons" "Temperature" "K" --suffix "_ch" --version "Initial X1" --comment "Read noise measurements from FPM S/N 104 end item data package"
make_measurements_fits.py --overwrite ../../data/amplifiers/read_noiseSW.txt "Read noise" "electrons" "Temperature" "K" --suffix "_ch" --version "Initial X1" --comment "Read noise measurements from FPM S/N 105 end item data package"

# Now convert the ASCII dark current measurements to FITS. There is a
# separate file for each of the MIRI focal plane modules.
# Normalize the dark current to 1.0 at 6.7K.
#
make_measurements_fits.py --overwrite ../../data/detector/dark_currentIM.txt "Dark current" "electrons" "Temperature" "K" --interptype LOGLIN --normalizeat 6.7 --version "Initial X1" --comment "Dark current measurement from FPM S/N 106 end item data package"
make_measurements_fits.py --overwrite ../../data/detector/dark_currentLW.txt "Dark current" "electrons" "Temperature" "K" --interptype LOGLIN --normalizeat 6.7 --version "Initial X1" --comment "Dark current measurement from FPM S/N 104 end item data package"
make_measurements_fits.py --overwrite ../../data/detector/dark_currentSW.txt "Dark current" "electrons" "Temperature" "K" --interptype LOGLIN --normalizeat 6.7 --version "Initial X1" --comment "Dark current measurement from FPM S/N 105 end item data package"

