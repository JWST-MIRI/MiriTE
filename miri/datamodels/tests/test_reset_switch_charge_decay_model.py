#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_reset_switch_charge_decay_model - Contains the unit tests for the classes
in the datamodels.miri_reset_switch_charge_decay_model module.

:History:

19 Nov 2015: Created
15 Jun 2017: Do not set observation metadata. It isn't
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
24 Jul 2020: Reduced the number of tests to remove the dependency on the actual
             structure of the data model, since it changes frequently.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_reset_switch_charge_decay_model import MiriResetSwitchChargeDecayModel
from miri.datamodels.tests.util import assert_recarray_equal


class TestMiriResetSwitchChargeDecayModel(unittest.TestCase):
    
    # Test the MiriResetSwitchChargeDecayModel class.
    
    def setUp(self):
        # Create a minimal data product.
        self.dataproduct = MiriResetSwitchChargeDecayModel( )
        self.testfile = "MiriResetSwitchChargeDecayModel_test.fits"
        
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
        # It must be possible to create an empty data product and fill
        # in its contents later.
        nulldp = MiriResetSwitchChargeDecayModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        del nulldp, descr1

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.rscd_table)
        # NOTE: Test that the contents are equal is removed.
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriResetSwitchChargeDecayModel(self.testfile) as readback:
                self.assertIsNotNone(readback.rscd_table)
#                self.assertEqual( len(self.dataproduct.rscd_table),
#                                  len(readback.rscd_table) )
#                original = np.asarray(self.dataproduct.rscd_table)
#                duplicate = np.asarray(readback.rscd_table)
#                assert_recarray_equal(original, duplicate)
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
