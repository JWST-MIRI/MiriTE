#!/bin/sh
#
# Run make_filters_fits.py on ASCII files contained in
# miri/datamodels/data/filters to generated the corresponding
# fits files.
#

echo 'Example filter: IM-01_F560EX'
make_filters_fits.py ../data/example_filter.txt -c F560EX BP MIRI IMG 5.6 1.2 um 7.0 1.3
