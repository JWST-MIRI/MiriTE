#!/bin/sh
#
# Run the scasim/make_sca_subarray.py utility to generate a
# test illumination for each unique subarray.
#
# This script assumes that make_sca_subarray.py has been installed
# into the PATH by running "setup.py install"
#

echo ""
echo "Creating test subarray data..."
#make_sca_subarray.py --overwrite 1024 1032 subarray_FULL.fits --subarray FULL
make_sca_subarray.py --overwrite 1024 1032 subarray_BRIGHTSKY.fits --subarray BRIGHTSKY
make_sca_subarray.py --overwrite 1024 1032 subarray_SUB256.fits --subarray SUB256
make_sca_subarray.py --overwrite 1024 1032 subarray_SUB128.fits --subarray SUB128
make_sca_subarray.py --overwrite 1024 1032 subarray_SUB64.fits --subarray SUB64
make_sca_subarray.py --overwrite 1024 1032 subarray_SLITLESSPRISM.fits --subarray SLITLESSPRISM
make_sca_subarray.py --overwrite 1024 1032 subarray_MASK1065.fits --subarray MASK1065
make_sca_subarray.py --overwrite 1024 1032 subarray_MASKLYOT.fits --subarray MASKLYOT
