#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_readnoise_model - Contains the unit tests for the classes
in the datamodels.miri_readnoise_model module.

:History:

10 Oct 2014: Created.
31 Oct 2014: Removed ERR and DQ arrays to bring the tests up to date with
             changes made to the data model on 17 Oct 2014.
11 Sep 2015: Minimal data object creation test added.
07 Oct 2015: Made exception catching Python 3 compatible.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_readnoise_model import MiriReadnoiseModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal


class TestMiriReadnoiseModel(unittest.TestCase):
    
    # Test the MiriReadnoiseModel class.
    
    def setUp(self):
        # Create a typical readnoise data product.
        data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
        self.dataproduct = MiriReadnoiseModel( data=data3x3 )
        self.testfile = "MiriReadnoiseModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        # Remove temporary file, if able to.
        if os.path.isfile(self.testfile):
            try:
                os.remove(self.testfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.testfile + \
                    "\n   " + str(e)
                warnings.warn(strg)

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
        # It must be possible to create an empty data product.    
        nullproduct = MiriReadnoiseModel( )
        del nullproduct
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriReadnoiseModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data'] )
                del readback
        
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


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
