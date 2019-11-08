#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_reset_model - Contains the unit tests for the classes
in the datamodels.miri_reset_model module.

:History:

10 Oct 2014: Created.
07 Oct 2015: Made exception catching Python 3 compatible.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
07 Oct 2019: FIXME: dq_def removed from unit tests until data corruption
             bug fixed (Bug 589).

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_reset_model import \
    reset_reference_flags, MiriResetModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal


class TestMiriResetModel(unittest.TestCase):
    
    # Test the MiriResetModel class.
    
    def setUp(self):
        # Create a typical reset data product.
        data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
        err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
        data3x3x2 = [data3x3, data3x3]
        err3x3x2 = [err3x3, err3x3]
        data3x3x2x2 = [data3x3x2, data3x3x2]
        err3x3x2x2 = [err3x3x2, err3x3x2]
        dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
        dqdef = reset_reference_flags
        self.dataproduct = MiriResetModel( data=data3x3x2x2, err=err3x3x2x2,
                                           dq=dq3x3, dq_def=dqdef )
        self.testfile = "MiriResetModel_test.fits"
        
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
        # It must be possible to create an empty reset model.
        dataproduct = MiriResetModel(  )
        descr = str(dataproduct)
        self.assertIsNotNone(descr)
        del descr, dataproduct
        
        # Attempting to create a reset model other than 4-D must fail.
        data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
        err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
        dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
        dqdef = reset_reference_flags
        self.assertRaises(ValueError, MiriResetModel, data=data3x3, err=err3x3,
                          dq=dq3x3, dq_def=dqdef )
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriResetModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'])
                # FIXME: removed dq_def until data corruption bug fixed. Bug 589
                #                       tables='dq_def' )
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
