#!/usr/bin/env python

"""
Module test_miri_filters - unit tests for the MiriFilter classes.

:History:

11 Jul 2014: Created to replace the test_filters.py module, which was
             based on the old Filters class.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
08 Sep 2015: Made compatible with Python 3
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings
import numpy as np

from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal
from miri.datamodels.miri_filters import MiriFilter, \
    MiriBandPassFilter, MiriQuantumEfficiency

class TestMiriFilter(unittest.TestCase):

    def setUp(self):        
        # Create a MiriFilter object containing test data
        self.transmissions = [
                              (0.5, 0.5),
                              (1.0, 0.5),
                              (1.5, 0.5),
                              (2.0, 0.5),
                              (2.5, 0.5),
                              (3.0, 0.5),
                              (3.5, 0.5),
                              (4.0, 0.5),
                              (4.5, 0.5),
                              (5.0, 0.5),
                              (5.5, 0.5),
                              (6.0, 0.5),
                              (6.5, 0.5),
                              (7.0, 0.5),
                              (7.5, 0.5),
                              (8.0, 0.5),
                              (8.5, 0.5),
                              (9.0, 0.5),
                              (9.5, 0.5),
                              (10.0, 0.5)]
        self.filt = MiriFilter(filter_table=self.transmissions,
                               filter_name='ANY', filter_type='ANY')
        # Add some typical metadata
        self.filt.set_instrument_metadata(detector='MIRIMAGE', modelnam='FM',
                                filt='ANY', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=14.0,
                                detector_temperature=6.7)

        # Name of temporary file for testing FITS I/O.
        self.tempfile = 'test_miri_filter.fits'

    def tearDown(self):
        # Clean temporary files.
        if os.path.isfile(self.tempfile):
            try:
                os.remove(self.tempfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.tempfile + \
                        "\n   " + str(e)
                warnings.warn(strg)
        # Clean python variables
        del self.filt, self.transmissions

    def test_creation(self):
        # TBD
        pass
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr = str(self.filt)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.filt)
        self.assertIsNotNone(descr)
        del descr

        descr = str(self.filt.transmission)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.filt.wavelength)
        self.assertIsNotNone(descr)
        del descr


    def test_apply(self):
        # Test that a constant transmission is correctly applied
        flux = np.linspace(120.0, 150.0, len(self.transmissions))
        # The result should be the same as the efficiency array at those
        # corresponding wavelength values times the same constant used
        # at setUp.
        result = self.filt.apply_filter(flux)
        self.assertTrue(np.allclose(0.5*flux, result))
        
        # Same test but a wavelength array is provided and interpolation is
        # required
        wave = np.linspace(1.0, 9.0, 50)
        flux = np.linspace(120.0, 150.0, 50)
        # The result should be the same as the efficiency array at those
        # corresponding wavelength values times the same constant used
        # at setUp.
        result = self.filt.apply_filter(flux, wave)
        self.assertTrue(np.allclose(0.5*flux, result))

    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.filt.save(self.tempfile, overwrite=True)
            with MiriFilter(self.tempfile) as readback:
                assert_products_equal( self, self.filt, readback,
                                       arrays=[], tables='filter_table' )
                del readback


class TestMiriBandPassFilter(unittest.TestCase):

    def setUp(self):        
        # Create a MiriBandPassFilter object containing test data
        self.transmissions = [
                    (5.80, 0.8698470),
                    (5.81, 0.8759494),
                    (5.82, 0.8944225),
                    (5.83, 0.8899569),
                    (5.84, 0.8760563),
                    (5.85, 0.8726164),
                    (5.86, 0.8782486),
                    (5.87, 0.8753881),
                    (5.88, 0.8844002),
                    (5.89, 0.8682995),
                    (5.90, 0.8495247),
                    (5.91, 0.8289118),
                    (5.92, 0.8211463),
                    (5.93, 0.8199366),
                    (5.94, 0.8202344),
                    (5.95, 0.7952070),
                    (5.96, 0.7884885),
                    (5.97, 0.7938501),
                    (5.98, 0.7938051),
                    (5.99, 0.8033671),
                    (6.00, 0.7985086)
                    ]
        self.filt = MiriBandPassFilter(filter_table=self.transmissions,
                               filter_name='F560W', wavecent=5.6, fwhm=1.2)
        # Add some typical metadata
        self.filt.set_instrument_metadata(detector='MIRIMAGE', modelnam='FM',
                                filt='F560W', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=14.0,
                                detector_temperature=6.7)

        # Name of temporary file for testing FITS I/O.
        self.tempfile = 'test_miri_bandpass_filter.fits'

    def tearDown(self):
        # Clean temporary files.
        if os.path.isfile(self.tempfile):
            try:
                os.remove(self.tempfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.tempfile + \
                        "\n   " + str(e)
                warnings.warn(strg)
        # Clean python variables
        del self.filt, self.transmissions

    def test_creation(self):
        # TBD
        pass
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr = str(self.filt)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.filt)
        self.assertIsNotNone(descr)
        del descr

        descr = str(self.filt.transmission)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.filt.wavelength)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.filt.meta.instrument.filter_wavecent)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.filt.meta.instrument.filter_fwhm)
        self.assertIsNotNone(descr)
        del descr        

    def test_apply(self):
        # Test that a non-constant transmission is applied without problem.
        flux = np.linspace(120.0, 150.0, len(self.transmissions))
        # The result should at least be the same length.
        result = self.filt.apply_filter(flux)
        self.assertEqual(len(flux), len(result))
        
        # Same test but a wavelength array is provided and interpolation is
        # required
        wave = np.linspace(5.85, 5.95, 50)
        flux = np.linspace(120.0, 150.0, 50)
        # The result should at least be the same length.
        result = self.filt.apply_filter(flux, wave)
        self.assertEqual(len(flux), len(result))

    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.filt.save(self.tempfile, overwrite=True)
            with MiriBandPassFilter(self.tempfile) as readback:
                assert_products_equal( self, self.filt, readback,
                                       arrays=[], tables='filter_table' )
                del readback


class TestMiriQuantumEfficiency(unittest.TestCase):

    def setUp(self):        
        # Create a MiriQuantumEfficiency object containing test data
        self.efficiency = [
                    (5.80, 0.8698470),
                    (5.81, 0.8759494),
                    (5.82, 0.8944225),
                    (5.83, 0.8899569),
                    (5.84, 0.8760563),
                    (5.85, 0.8726164),
                    (5.86, 0.8782486),
                    (5.87, 0.8753881),
                    (5.88, 0.8844002),
                    (5.89, 0.8682995),
                    (5.90, 0.8495247),
                    (5.91, 0.8289118),
                    (5.92, 0.8211463),
                    (5.93, 0.8199366),
                    (5.94, 0.8202344),
                    (5.95, 0.7952070),
                    (5.96, 0.7884885),
                    (5.97, 0.7938501),
                    (5.98, 0.7938051),
                    (5.99, 0.8033671),
                    (6.00, 0.7985086)
                    ]
        self.filt = MiriQuantumEfficiency(qe_table=self.efficiency,
                               detector='MIRIMAGE', temperature=6.7)
        # Add some typical metadata
        self.filt.set_instrument_metadata(detector='MIRIMAGE', modelnam='FM',
                                filt='ANY', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=14.0,
                                detector_temperature=6.7)

        # Name of temporary file for testing FITS I/O.
        self.tempfile = 'test_miri_quantum_efficiency.fits'

    def tearDown(self):
        # Clean temporary files.
        if os.path.isfile(self.tempfile):
            try:
                os.remove(self.tempfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.tempfile + \
                        "\n   " + str(e)
                warnings.warn(strg)
        # Clean python variables
        del self.filt, self.efficiency

    def test_creation(self):
        # TBD
        pass
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr = str(self.filt)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.filt)
        self.assertIsNotNone(descr)
        del descr

        descr = str(self.filt.efficiency)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.filt.wavelength)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.filt.meta.instrument.detector_temperature)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.filt.meta.instrument.filter_fwhm)
        self.assertIsNotNone(descr)
        del descr        

    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.filt.save(self.tempfile, overwrite=True)
            with MiriBandPassFilter(self.tempfile) as readback:
                assert_products_equal( self, self.filt, readback,
                                       arrays=[], tables='filter_table' )
                del readback

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
