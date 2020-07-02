#!/bin/sh
#
# Run the scasim/make_sca_calibration.py utility to generate
# some dummy calibration files. The command parameters can be varied
# to control the properties of the calibrations.
#
# This script assumes that make_sca_calibration.py has been installed
# into the PATH by running "setup.py install"
#

echo ""
echo "Creating test bad pixel masks..."
# Create some example bad pixel masks containing 100 randomly placed
# bad pixels.
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsEXAMPLEIM.fits --pattern BADPIXEL --nbad 100 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsEXAMPLELW.fits --pattern BADPIXEL --nbad 100 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsEXAMPLESW.fits --pattern BADPIXEL --nbad 100 --scaid MIRIFUSHORT

# Create some perfect bad pixel masks with no bad pixels.
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsPERFECTIM.fits --pattern BADPIXEL --nbad 0 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsPERFECTLW.fits --pattern BADPIXEL --nbad 0 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 bad_pixelsPERFECTSW.fits --pattern BADPIXEL --nbad 0 --scaid MIRIFUSHORT

echo ""
echo "Creating test dark maps..."
# Create some example dark maps with multiplier varying from x 0.5 to x 2.0
# and 50 randomly placed hot pixels
make_sca_calibration.py --overwrite 1024 1032 dark_mapEXAMPLEIM.fits --pattern DARKMAP --nbad 50 --min 0.5 --max 2.0 --badvalue 10000.0 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 dark_mapEXAMPLELW.fits --pattern DARKMAP --nbad 50 --min 0.5 --max 2.0 --badvalue 10000.0 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 dark_mapEXAMPLESW.fits --pattern DARKMAP --nbad 50 --min 0.5 --max 2.0 --badvalue 10000.0 --scaid MIRIFUSHORT

# Create some perfect dark maps with no variation in dark currect and no hot pixels.
make_sca_calibration.py --overwrite 1024 1032 dark_mapPERFECTIM.fits --pattern DARKMAP --nbad 0 --min 1.0 --max 1.0 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 dark_mapPERFECTLW.fits --pattern DARKMAP --nbad 0 --min 1.0 --max 1.0 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 dark_mapPERFECTSW.fits --pattern DARKMAP --nbad 0 --min 0.5 --max 1.0 --scaid MIRIFUSHORT

echo ""
echo "Creating test flat-fields..."
# Create some perfect flat-fields
make_sca_calibration.py --overwrite 1024 1032 flat_fieldPERFECTIM.fits --pattern CONSTANT --constant 1.0 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 flat_fieldPERFECTLW.fits --pattern CONSTANT --constant 1.0 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 flat_fieldPERFECTSW.fits --pattern CONSTANT --constant 1.0 --scaid MIRIFUSHORT

# Create some flat-fields containing a slope
make_sca_calibration.py --overwrite 1024 1032 flat_fieldSLOPEIM.fits --pattern SLOPE --constant 1.0 --rowslope 0.00005 --colslope -0.00005 --scaid MIRIMAGE
make_sca_calibration.py --overwrite 1024 1032 flat_fieldSLOPELW.fits --pattern SLOPE --constant 1.0 --rowslope -0.0001 --colslope 0.0 --scaid MIRIFULONG
make_sca_calibration.py --overwrite 1024 1032 flat_fieldSLOPESW.fits --pattern SLOPE --constant 1.0 --rowslope -0.0001 --colslope 0.00005 --scaid MIRIFUSHORT
