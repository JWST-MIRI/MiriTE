#!/bin/sh
#
# Run make_filters_fits.py on ASCII files contained in miri/simulators/data/filters to
# generated the corresponding fits files.
#

echo 'IM-01_F560W'
make_filters_fits.py ../data/filters/IM-01_F560W.txt -c F560W BandPass MIRI IMG 5.6 1.2 um 7.0 1.3
echo 'IM-02_F770W'
make_filters_fits.py ../data/filters/IM-02_F770W.txt -c F770W BandPass MIRI IMG 7.7 2.2 um 7.0 1.3
echo 'IM-03_F1000W'
make_filters_fits.py ../data/filters/IM-03_F1000W.txt -c F1000W BandPass MIRI IMG 10.0 2.0 um 7.0 1.3
echo 'IM-04_F1130W'
make_filters_fits.py ../data/filters/IM-04_F1130W.txt -c F1130W BandPass MIRI IMG 11.3 2.2 um 7.0 1.3
echo 'IM-05_F1280W'
make_filters_fits.py ../data/filters/IM-05_F1280W.txt -c F1280W BandPass MIRI IMG 12.8 2.4 um 7.0 1.3
echo 'IM-06_F1500W'
make_filters_fits.py ../data/filters/IM-06_F1500W.txt -c F1500W BandPass MIRI IMG 15.0 3.0 um 7.0 1.3
echo 'IM-07_F1800W'
make_filters_fits.py ../data/filters/IM-07_F1800W.txt -c F1800W BandPass MIRI IMG 18.0 3.0 um 7.0 1.3
echo 'IM-08_F2100W'
make_filters_fits.py ../data/filters/IM-08_F2100W.txt -c F2100W BandPass MIRI IMG 21.0 5.0 um 7.0 1.3
echo 'IM-09_F2550W'
make_filters_fits.py ../data/filters/IM-09_F2550W.txt -c F2550W BandPass MIRI IMG 25.5 5.0 um 7.0 1.3
echo 'IM-10_F2550WR'
make_filters_fits.py ../data/filters/IM-10_F2550W.txt -c F2550WR BandPass MIRI IMG 25.5 5.0 um 10.0 2.3
