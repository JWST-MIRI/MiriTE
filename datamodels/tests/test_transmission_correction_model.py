#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_transmission_correction_model - Contains the unit tests for the classes
in the datamodels.miri_transmission_correction_model module.

:History:

10 May 2017: Original version.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_transmission_correction_model import MiriMrsTransmissionCorrectionModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriMrsTransmissionCorrectionModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a MiriFilter object containing test data
        self.tracorr_table = [(1, 4.87, 7.76, 96.7, 93.3, 91.0, 91.6),
                              (2, 7.45, 11.87, 96.2, 92.2, 91.9, 90.5)]
        self.dataproduct = MiriMrsTransmissionCorrectionModel( \
                                        tracorr_table=self.tracorr_table)
        # Add some typical metadata
        self.dataproduct.set_instrument_metadata(detector='MIRIMAGE',
                                modelnam='FM',
                                filt='ANY', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=14.0,
                                detector_temperature=6.7)

        # Name of temporary file for testing FITS I/O.
        self.testfile = "MiriMrsTransmissionCorrectionModel_test.fits"
        
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
        del self.dataproduct, self.tracorr_table

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
        class_names = list(MiriMrsTransmissionCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('tracorr_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")
 
        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsTransmissionCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.tracorr_table = self.tracorr_table
        self.assertIsNotNone(nulldp.tracorr_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.tracorr_table)
        self.assertEqual( len(self.dataproduct.tracorr_table),
                          len(datacopy.tracorr_table) )
        table1 = np.asarray(self.dataproduct.tracorr_table)
        table2 = np.asarray(datacopy.tracorr_table)
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

        descr = str(self.dataproduct.tracorr_table)
        self.assertIsNotNone(descr)
        del descr        
 
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
         
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsTransmissionCorrectionModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=[], tables='tracorr_table' )
                del readback


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
