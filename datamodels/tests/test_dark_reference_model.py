#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_dark_reference_model - Contains the unit tests for the classes
in the datamodels.miri_dark_reference_model module.

:History:

16 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
26 Apr 2013: File closing problem has returned!
29 Jul 2013: stats() method added.
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
12 Sep 2013: Test that the data product can be copied successfully.
09 Jul 2014: field_def changed to dq_def.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
31 Oct 2014: Removed FITERR array to bring the tests up to date with
             changes made to the data model on 28 Oct 2014.
04 Nov 2014: Removed the legacy call to dataproduct.storage.close(),
             which caused a problem with r3073 of jwst_lib.
16 Jan 2015: A wrongly shaped array returns a TypeError rather than ValueError.
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

#import numpy as np

from miri.datamodels.miri_dark_reference_model import \
    dark_reference_flags, MiriDarkReferenceModel
from miri.datamodels.tests.util import assert_products_equal


class TestMiriDarkReferenceModel(unittest.TestCase):
    
    # Test the MiriDarkReferenceModel class.
    # Most of the relevant tests are already done in TestMiriMeasuredRampModel.
    # Only the additional tests specific to MiriDarkReferenceModel are
    # included here.
    
    def setUp(self):
        # Create a typical dark reference product.
        self.a1 = [[0.01,0.02,0.03,0.04], [0.05,0.05,0.07,0.08], [0.09,0.1,0.11,0.12]]
        self.b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        self.c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.dqdef = dark_reference_flags
        self.acube = [self.a1,self.a1,self.a1]
        self.bcube = [self.b1,self.b1,self.b1]
#         self.ccube = [self.c1,self.c1,self.c1]
#         self.ahyper = [self.acube,self.acube]#
#         self.bhyper = [self.bcube,self.bcube]
#         self.chyper = [self.ccube,self.ccube]
        self.dataproduct = MiriDarkReferenceModel(data=self.acube,
                                                  err=self.bcube,
                                                  dq=self.c1,
                                                  dq_def=self.dqdef)
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
        self.testfile = "MiriDarkReferenceModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
#         del self.ahyper, self.bhyper, self.chyper
        del self.acube, self.bcube
        del self.a1, self.b1, self.c1
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
        # 3-D data, 3-D err and 2-D dq must be accepted
        testproduct = MiriDarkReferenceModel(data=self.acube,
                                            err=self.bcube,
                                            dq=self.c1)
        descr = str(testproduct)
        del testproduct, descr
        
        # 3-D data withour err must be acceptable
        testproduct = MiriDarkReferenceModel(data=self.acube, dq=self.c1)
        descr = str(testproduct)
        del testproduct, descr

        # 3-D data withour err or dq be acceptable
        testproduct = MiriDarkReferenceModel(data=self.acube)
        descr = str(testproduct)
        del testproduct, descr
        
        # The main data array must be 3-D or 4-D. Other shapes must be rejected.
        data1d = [1, 2, 3, 4, 5, 6]
        self.assertRaises(TypeError, MiriDarkReferenceModel, data=data1d,
                          err=self.b1, dq=self.c1)

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq'],
                               tables='dq_def' )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriDarkReferenceModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables='dq_def' )
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
        descr = self.dataproduct.stats()
        self.assertIsNotNone(descr)
        del descr
        
        # Attempt to access the SCI, DQ and FITERR arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
