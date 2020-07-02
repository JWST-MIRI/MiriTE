#!/bin/sh
#
# Some examples of how to run the SCA simulator from the command line.        |
# This file also be used to test the installation of scasim.
#
# This script assumes that scasim.py has been installed
# into the PATH by running "setup.py install"
#
# An exposure on the 80x64 test data included with the release. The detector
# readout mode is SLOW and the exposure will be made of 12 groups and 2
# integrations. (It is better to give the number of groups explicitly
# because the 80x60 data has a very small frame time and can easily
# convert a small integration time into hundreds of groups). The detector
# temperature is 7.0K and the simulation will assume a cosmic ray
# environment of SOLAR_MAX. There will be verbose output (and data plots),
# and any existing file will be overwritten. Level 1 FITS data will be
# saved to SCATestOutput1.fits.
#
python scasim.py --verbose --overwrite ../data/SCATestInput80x64.fits SCATestOutput1.fits --rdmode SLOW --ngroups 12 --nints 2 --temperature 7.0 --crmode SOLAR_MAX

# As above but with FAST readout mode. Note that you don't need to specify "python".
#
scasim.py --verbose --overwrite ../data/SCATestInput80x64.fits SCATestOutput2.fits --rdmode FAST --ngroups 12 --nints 2 --temperature 7.0 --crmode SOLAR_MAX

# As above but with FASTGRPAVG readout mode. The resulting file will contain
# 3 sets of (4) averaged groups. By this point we will be fed up with
# verbose mode, so switch to normal output.
#
scasim.py --overwrite ../data/SCATestInput80x64.fits SCATestOutput3.fits --rdmode FASTGRPAVG --ngroups 12 --nints 2 --temperature 7.0 --crmode SOLAR_MAX

# As above but with FASTINTAVG readout mode. The 2 integrations specified
# will be rounded up to 4 and the file will contain 1 set of (4) averaged
# integrations. 
#
scasim.py --overwrite ../data/SCATestInput80x64.fits SCATestOutput4.fits --rdmode FASTINTAVG --ngroups 12 --nints 2 --temperature 7.0 --crmode SOLAR_MAX

# Back to the first combination but with the cosmic ray simulation using
# SOLAR_MIN and the detector temperature raised to 8.0K. This time
# there are 10 groups and 1 integration.
#
scasim.py --overwrite ../data/SCATestInput80x64.fits SCATestOutput5.fits --rdmode SLOW --ngroups 10 --nints 1 --temperature 8.0 --crmode SOLAR_MIN

# As above combination but with FAST readout and cosmic ray simulation
# using SOLAR_FLARE.
#
scasim.py --overwrite ../data/SCATestInput80x64.fits SCATestOutput6.fits --rdmode FAST --ngroups 10 --nints 1 --temperature 8.0 --crmode SOLAR_FLARE

# Turn off the cosmic ray simulation and the output altogether.
#
scasim.py --overwrite --silent ../data/SCATestInput80x64.fits SCATestOutput7.fits --rdmode FAST --ngroups 10 --nints 1 --temperature 8.0 --crmode NONE

# Back to the first command but specify a small integration time and let scasim calculate the
# number of groups (careful, this could get quite large for this small data set).
#
scasim.py --overwrite ../data/SCATestInput80x64.fits SCATestOutput8.fits --rdmode SLOW --inttime 0.1 --nints 1 --temperature 7.0 --crmode SOLAR_MAX

# Create some full frame test calibration data using make_sca_calibration.py
make_sca_calibration.py --overwrite 1024 1024 SCATestFullFrame1024.fits --pattern CONSTANT --constant 1.0

#
# A 120 second exposure on the full frame test data created above, with
# the detectors at 10K during solar minimum and the detectors reading out
# in SLOW mode with 1 integration. Any existing file will be overrwitten.
# Level 1 FITS data will at first be saved to SCATestOutput_FULL.fits
#
# With full frame data we can try out all the various subarray modes.
#
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_FULL.fits --rdmode SLOW --subarray FULL --inttime 120.0 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_MASK1550.fits --rdmode SLOW --subarray MASK1550 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_MASK1140.fits --rdmode SLOW --subarray MASK1140 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_MASK1065.fits --rdmode SLOW --subarray MASK1065 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_MASKLYOT.fits --rdmode SLOW --subarray MASKLYOT --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_BRIGHTSKY.fits --rdmode SLOW --subarray BRIGHTSKY --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_SUB256.fits --rdmode SLOW --subarray SUB256 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_SUB128.fits --rdmode SLOW --subarray SUB128 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_SUB64.fits --rdmode SLOW --subarray SUB64 --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
scasim.py --overwrite ../data/SCATestHorseHead1024.fits SCATestOutput_SLITLESSPRISM.fits --rdmode SLOW --subarray SLITLESSPRISM --ngroups 10 --nints 1 --temperature 7.0 --crmode NONE
