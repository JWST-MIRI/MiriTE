#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_fringe_frequencies_model - Contains the unit tests for the classes
in the datamodels.miri_fringe_frequencies_model module.

:History:

10 May 2017: Original version.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Oct 2018: 'N/A' used as a metadata wildcard instead of 'ANY'.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_fringe_frequencies_model import \
    MiriMrsFringeFrequenciesModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriMrsFringeFrequenciesModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a MiriFilter object containing test data
        self.fringefreq_table = [('N/A', 11.0,  0.1, 42.0),
                      ('N/A', 12.0,  0.2, 22.0),
                      ('N/A', 13.0,  0.3, 32.0),
                      ('N/A', 14.0,  0.4, 12.0)]
        self.dataproduct = MiriMrsFringeFrequenciesModel( \
                                    fringefreq_table=self.fringefreq_table)
        # Add some typical metadata
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                         ccc_pos='OPEN', channel='N/A',
                                         band='N/A')
        self.dataproduct.set_subarray_metadata('FULL')
        self.dataproduct.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2015-11-20',
                                           description='Test data')

        # Name of temporary file for testing FITS I/O.
        self.testfile = "MiriMrsFringeFrequenciesModel_test.fits"
        
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
        del self.dataproduct, self.fringefreq_table

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
        class_names = list(MiriMrsFringeFrequenciesModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('fringefreq_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")
 
        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsFringeFrequenciesModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.fringefreq_table = self.fringefreq_table
        self.assertIsNotNone(nulldp.fringefreq_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.fringefreq_table)
        self.assertEqual( len(self.dataproduct.fringefreq_table),
                          len(datacopy.fringefreq_table) )
        table1 = np.asarray(self.dataproduct.fringefreq_table)
        table2 = np.asarray(datacopy.fringefreq_table)
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
 
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
         
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsFringeFrequenciesModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=[], tables='fringefreq_table' )
                del readback


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
