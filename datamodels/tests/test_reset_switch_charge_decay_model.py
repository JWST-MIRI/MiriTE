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

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_reset_switch_charge_decay_model import MiriResetSwitchChargeDecayModel
from miri.datamodels.tests.util import assert_recarray_equal


class TestMiriResetSwitchChargeDecayModel(unittest.TestCase):
    
    # Test the MiriResetSwitchChargeDecayModel class.
    
    def setUp(self):
        # Create a typical rscd product.
        self.rscddata = [('FULL',          'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('FULL',          'FAST', 'EVEN', 1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('FULL',          'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('FULL',          'SLOW', 'EVEN', 1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1065',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1065',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1140',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1140',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1550',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASK1550',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASKLYOT',      'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('MASKLYOT',      'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('BRIGHTSKY',     'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('BRIGHTSKY',     'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('SLITLESSPRISM', 'FAST', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ('SLITLESSPRISM', 'SLOW', 'ODD',  1.0e-1,  2.0e-2, 3.0e-3, 4.0e-4),
                ]
        self.dataproduct = MiriResetSwitchChargeDecayModel( rscd_table=self.rscddata,
                                                            detector='MIRIMAGE' )
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIMAGE',
                                                 ccc_pos='CLOSED',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.testfile = "MiriResetSwitchChargeDecayModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        # TBD
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
        # Check that the field names in the class variable are the same
        # as the ones declared in the schema.
        class_names = list(MiriResetSwitchChargeDecayModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('rscd_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later.
        nulldp = MiriResetSwitchChargeDecayModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.rscd_table = self.rscddata
        self.assertIsNotNone(nulldp.rscd_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.rscd_table)
        self.assertEqual( len(self.dataproduct.rscd_table),
                          len(datacopy.rscd_table) )
        table1 = np.asarray(self.dataproduct.rscd_table)
        table2 = np.asarray(datacopy.rscd_table)
        assert_recarray_equal(table1, table2)
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
                self.assertEqual( len(self.dataproduct.rscd_table),
                                  len(readback.rscd_table) )
                original = np.asarray(self.dataproduct.rscd_table)
                duplicate = np.asarray(readback.rscd_table)
                assert_recarray_equal(original, duplicate)
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
