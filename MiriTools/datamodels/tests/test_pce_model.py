#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_pce_model - Contains the unit tests for the classes
in the datamodels.miri_pce_model module.

:History:

11 Jul 2014: Original version as test_miri_filters.py.
22 Jun 2016: Converted to test the new PCE model, MiriPceModel.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_pce_model import MiriPceModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriPceModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a MiriFilter object containing test data
        self.pce_table = [(0.5, 0.5, 1.0),
                          (1.0, 0.5, 1.0),
                          (1.5, 0.5, 1.0),
                          (2.0, 0.5, 1.0),
                          (2.5, 0.5, 1.0),
                          (3.0, 0.5, 1.0),
                          (3.5, 0.5, 1.0),
                          (4.0, 0.5, 1.0),
                          (4.5, 0.5, 1.0),
                          (5.0, 0.5, 1.0),
                          (5.5, 0.5, 1.0),
                          (6.0, 0.5, 1.0),
                          (6.5, 0.5, 1.0),
                          (7.0, 0.5, 1.0),
                          (7.5, 0.5, 1.0),
                          (8.0, 0.5, 1.0),
                          (8.5, 0.5, 1.0),
                          (9.0, 0.5, 1.0),
                          (9.5, 0.5, 1.0),
                          (10.0, 0.5, 1.0)]
        self.dataproduct = MiriPceModel(pce_table=self.pce_table,
                                           component='ANY', detector='ANY')
        # Add some typical metadata
        self.dataproduct.set_instrument_metadata(detector='MIRIMAGE',
                                modelnam='FM',
                                filt='ANY', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=14.0,
                                detector_temperature=6.7)

        # Name of temporary file for testing FITS I/O.
        self.testfile = "MiriPceModel_test.fits"
        
    def tearDown(self):
        # Clean temporary files.
        if os.path.isfile(self.testfile):
            try:
                os.remove(self.testfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.testfile + \
                        "\n   " + str(e)
                warnings.warn(strg)
        # Clean python variables
        del self.dataproduct, self.pce_table

    def test_referencefile(self):
        # Check that the data product contains the standard
        # reference file metadata.
        type1 = self.dataproduct.meta.model_type
        type2 = self.dataproduct.meta.reftype
        self.assertIsNotNone(type1)
        self.assertIsNotNone(type2)
        pedigree = self.dataproduct.meta.pedigree
        self.assertIsNotNone(pedigree)

    def test_creation(self):
        # Check that the field names in the class variable are the same
        # as the ones declared in the schema.
        class_names = list(MiriPceModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('pce_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")
 
        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriPceModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.pce_table = self.pce_table
        self.assertIsNotNone(nulldp.pce_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.pce_table)
        self.assertEqual( len(self.dataproduct.pce_table),
                          len(datacopy.pce_table) )
        table1 = np.asarray(self.dataproduct.pce_table)
        table2 = np.asarray(datacopy.pce_table)
        assert_recarray_equal(table1, table2)
        del datacopy
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr

        descr = str(self.dataproduct.wavelength)
        self.assertIsNotNone(descr)
        del descr        
        descr = str(self.dataproduct.efficiency)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.conversion)
        self.assertIsNotNone(descr)
        del descr

    def test_apply(self):
        # Test that a constant efficiency is correctly applied
        flux = np.linspace(120.0, 150.0, len(self.pce_table))
        # The result should be the same as the efficiency array at those
        # corresponding wavelength values times the same constant used
        # at setUp.
        result = self.dataproduct.apply_filter(flux)
        self.assertTrue(np.allclose(0.5*flux, result))
         
        # Same test but a wavelength array is provided and interpolation is
        # required
        wave = np.linspace(1.0, 9.0, 50)
        flux = np.linspace(120.0, 150.0, 50)
        # The result should be the same as the efficiency array at those
        # corresponding wavelength values times the same constant used
        # at setUp.
        result = self.dataproduct.apply_filter(flux, wave)
        self.assertTrue(np.allclose(0.5*flux, result))
 
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
         
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriPceModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=[], tables='pce_table' )
                del readback


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
