#!/bin/sh
#
# Run SCASim to analyse the test subarray data stored in ../data.
#

echo ""
echo "Simulating test subarray data..."
scasim.py --detector MIRIMAGE --scale 200.0 --rdmode FAST --nints 1 --ngroups 30 --subarray BRIGHTSKY ../data/subarray_BRIGHTSKY.fits ../data/ramp_BRIGHTSKY.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 800.0 --rdmode FAST --nints 1 --ngroups 30 --subarray MASK1065 ../data/subarray_MASK1065.fits ../data/ramp_MASK1065.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 800.0 --rdmode FAST --nints 1 --ngroups 30 --subarray MASKLYOT ../data/subarray_MASKLYOT.fits ../data/ramp_MASKLYOT.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 800.0 --rdmode FAST --nints 1 --ngroups 30 --subarray SUB256 ../data/subarray_SUB256.fits ../data/ramp_SUB256.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 10000.0 --rdmode FAST --nints 1 --ngroups 40 --subarray SUB128 ../data/subarray_SUB128.fits ../data/ramp_SUB128.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 20000.0 --rdmode FAST --nints 1 --ngroups 40 --subarray SLITLESSPRISM ../data/subarray_SLITLESSPRISM.fits ../data/ramp_SLITLESSPRISM.fits --nodark --noflat --noreadnoise

echo "Simulating FULL to subarray data"
scasim.py --detector MIRIMAGE --scale 200.0 --rdmode FAST --nints 1 --ngroups 30 --subarray SUB256 ../data/subarray_FULL.fits ../data/ramp_FULL_256.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 200.0 --rdmode FAST --nints 1 --ngroups 30 --subarray SUB128 ../data/subarray_FULL.fits ../data/ramp_FULL_128.fits --nodark --noflat --noreadnoise

echo "Simulating subarray to FULL data"
scasim.py --detector MIRIMAGE --scale 200.0 --rdmode FAST --nints 1 --ngroups 30 --subarray FULL ../data/subarray_SUB256.fits ../data/ramp_256_FULL.fits --nodark --noflat --noreadnoise
scasim.py --detector MIRIMAGE --scale 200.0 --rdmode FAST --nints 1 --ngroups 30 --subarray FULL ../data/subarray_SUB128.fits ../data/ramp_128_FULL.fits --nodark --noflat --noreadnoise

