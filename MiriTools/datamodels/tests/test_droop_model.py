#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_droop_model - Contains the unit tests for the classes
in the datamodels.miri_droop_model module.

:History:

22 Feb 2013: Created.
26 Apr 2013: Corrected name of temporary file.
17 May 2013: Do not allow a blank table to be created.
22 Aug 2013: Check that the field names declared in the class variable
             match the schema.
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
12 Sep 2013: Test that the data product can be copied successfully.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: TYPE and REFTYPE are no longer identical.
11 Mar 2015: group_integration_time changed to group_time.
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

from miri.datamodels.miri_droop_model import MiriDroopModel
from miri.datamodels.tests.util import assert_recarray_equal


class TestMiriDroopModel(unittest.TestCase):
    
    # Test the MiriDroopModel class.
    
    def setUp(self):
        # Create a typical droop product.
        droopdata = [(0.012,  0.012)]
        self.dataproduct = MiriDroopModel( droop_table=droopdata )
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='CLOSED',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriDroopModel_test.fits"
        
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
        class_names = list(MiriDroopModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('droop_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later.
        droopdata = [(0.012,  0.012)]
        nulldp = MiriDroopModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.droop_table = droopdata
        self.assertIsNotNone(nulldp.droop_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.droop_table)
        self.assertEqual( len(self.dataproduct.droop_table),
                          len(datacopy.droop_table) )
        table1 = np.asarray(self.dataproduct.droop_table)
        table2 = np.asarray(datacopy.droop_table)
        assert_recarray_equal(table1, table2)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriDroopModel(self.testfile) as readback:
                self.assertIsNotNone(readback.droop_table)
                self.assertEqual( len(self.dataproduct.droop_table),
                                  len(readback.droop_table) )
                original = np.asarray(self.dataproduct.droop_table)
                duplicate = np.asarray(readback.droop_table)
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
