#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_jump_model - Contains the unit tests for the classes
in the datamodels.miri_jump_model module.

:History:

28 Mar 2014: Created.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
07 Oct 2015: Made exception catching Python 3 compatible.
15 Jun 2017: Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_jump_model import MiriJumpModel
from miri.datamodels.tests.util import assert_recarray_equal


class TestMiriJumpModel(unittest.TestCase):
    
    # Test the MiriJumpModel class.
    
    def setUp(self):
        # Create a typical jump product.
        finejumpdata = [('2PNTDIFF', 3.5),
                        ('YINT', 3.5)]
        crthresh = 5.0

        self.dataproduct = MiriJumpModel( finejump_table=finejumpdata,
                                          crthresh=crthresh )
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='CLOSED',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.testfile = "MiriJumpModel_test.fits"
        
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
        # Check that the field names in the class variable are the same
        # as the ones declared in the schema.
        class_names = list(MiriJumpModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('finejump_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later.
        jumpdata = [('2PNTDIFF', 3.7), ('YINT', 3.7)]
        nulldp = MiriJumpModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.finejump_table = jumpdata
        self.assertIsNotNone(nulldp.finejump_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.finejump_table)
        self.assertEqual( len(self.dataproduct.finejump_table),
                          len(datacopy.finejump_table) )
        table1 = np.asarray(self.dataproduct.finejump_table)
        table2 = np.asarray(datacopy.finejump_table)
        assert_recarray_equal(table1, table2)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriJumpModel(self.testfile) as readback:
                self.assertIsNotNone(readback.finejump_table)
                self.assertEqual( len(self.dataproduct.finejump_table),
                                  len(readback.finejump_table) )
                original = np.asarray(self.dataproduct.finejump_table)
                duplicate = np.asarray(readback.finejump_table)
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
