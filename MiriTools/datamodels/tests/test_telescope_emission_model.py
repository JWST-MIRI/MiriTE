#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_telescope_emission_model - Contains the unit tests for the classes
in the datamodels.miri_telescope_emission_model module.

:History:

16 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
26 Apr 2013: File closing problem has returned!
29 Jul 2013: stats() method added.
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Test that the data product can be copied successfully.
09 Jul 2014: field_def changed to dq_def.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
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

from miri.datamodels.miri_telescope_emission_model import \
    telem_reference_flags, MiriTelescopeEmissionModel
from miri.datamodels.tests.util import assert_products_equal


class TestMiriTelescopeEmissionModel(unittest.TestCase):
    
    # Test the MiriTelescopeEmissionModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriTelescopeEmissionModel are
    # included here.
    
    def setUp(self):
        # Create a typical telescope emission product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        dqdef = telem_reference_flags
        self.dataproduct = MiriTelescopeEmissionModel(data=a1, err=b1, dq=c1,
                                                      dq_def=dqdef,
                                                      filt='F560W',
                                                      temperature=40.0)
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIMAGE',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriTelescopeEmissionModel_test.fits"
        
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
        # It should be possible to set up an empty data product with
        # a specified shape. All three arrays should be initialised to
        # the same shape.
        emptydp = MiriTelescopeEmissionModel( (4,4) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (4,4))
        self.assertIsNotNone(emptydp.err)
        self.assertEqual(emptydp.err.shape, (4,4))
        self.assertIsNotNone(emptydp.dq)
        self.assertEqual(emptydp.dq.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr

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
            with MiriTelescopeEmissionModel(self.testfile) as readback:
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
        
        # Attempt to access the SCI, ERR and DQ arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.err)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
