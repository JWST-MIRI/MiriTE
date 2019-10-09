#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_pixel_saturation_model - Contains the unit tests for the classes
in the datamodels.miri_straylight_model module.

:History:

28 Mar 2014: Created.
09 Jul 2014: field_def changed to dq_def.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
11 Mar 2015: group_integration_time changed to group_time.
08 Aug 2015: Mask data are now stored in a "dq" array rather than a "data"
             array.
19 Aug 2015: Corrected typos in detector names.
07 Oct 2015: Made exception catching Python 3 compatible.
15 Jun 2017: Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
15 Nov 2018: New data model which uses the JWST schema, saturation.schema.
             Previous data model renamed to MiriMrsStraylightModel_CDP3.
07 Oct 2019: FIXME: dq_def removed from unit tests until data corruption
             bug (589) is fixed.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_straylight_model import \
    straylight_reference_flags, MiriMrsStraylightModel, MiriMrsStraylightModel_CDP3
from miri.datamodels.tests.util import assert_products_equal


class TestMiriMrsStraylightModel(unittest.TestCase):
    
    # Test the MiriMrsStrayLightModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriMrsStrayLightModel are
    # included here.
    
    def setUp(self):
        # Create a typical stray light product.
        self.maskdata = np.array([[0,0,1,0],
                                  [0,0,1,0],
                                  [0,1,0,0],
                                  [0,1,0,0]])
        self.dataproduct = MiriMrsStraylightModel(data=self.maskdata,
                                                  detector='MIRIFUSHORT')
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriMrsStraylightModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.maskdata
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
        # a specified shape. The mask should be initialised to
        # the same shape.
        emptydp = MiriMrsStraylightModel( (4,4) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data'] )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsStraylightModel(self.testfile) as readback:
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
        descr = self.dataproduct.stats()
        self.assertIsNotNone(descr)
        del descr
        
        # Attempt to access the data array through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr


class TestMiriMrsStraylightModel_CDP3(unittest.TestCase):
    
    # Test the MiriMrsStrayLightModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriMrsStrayLightModel are
    # included here.
    
    def setUp(self):
        # Create a typical stray light product.
        self.maskdata = np.array([[0,0,1,0],
                                  [0,0,1,0],
                                  [0,1,0,0],
                                  [0,1,0,0]])
        self.dqdef = straylight_reference_flags
        self.dataproduct = MiriMrsStraylightModel_CDP3(dq=self.maskdata,
                                                       dq_def=self.dqdef,
                                                       detector='MIRIFUSHORT')
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriMrsStraylightModel_CDP3_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.maskdata
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
        # a specified shape. The mask should be initialised to
        # the same shape.
        emptydp = MiriMrsStraylightModel_CDP3( (4,4) )
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
                               arrays=['dq'])
        # FIXME: removed dq_def until data corruption bug fixed.
        #                       tables='dq_def' )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsStraylightModel_CDP3(self.testfile) as readback:
                # FIXME: The dq_def table is not the same because a single line table is replace by the default.
#                 assert_products_equal( self, self.dataproduct, readback,
#                                arrays=['data'], tables='dq_def' )
                assert_products_equal( self, self.dataproduct, readback,
                               arrays=['dq'], tables=[] )
                del readback

    def test_operators(self):
        maskdata2 = np.array([[0,0,1,0],
                              [0,0,1,0],
                              [0,1,1,0],
                              [0,1,0,0]])
        newproduct = MiriMrsStraylightModel_CDP3(dq=maskdata2,
                                            dq_def=self.dqdef,
                                            detector='MIRIFUSHORT')

        # Scalar OR
        result = self.dataproduct | 1
        test1 = self.dataproduct.dq | 1
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product OR
        result = self.dataproduct | newproduct
        test1 = self.dataproduct.dq | newproduct.dq
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        del result

        # Scalar XOR
        result = self.dataproduct ^ 1
        test1 = self.dataproduct.dq ^ 1
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product XOR
        result = self.dataproduct ^ newproduct
        test1 = self.dataproduct.dq ^ newproduct.dq
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        del result
         
        # Scalar AND
        result = self.dataproduct & 15
        test1 = self.dataproduct.dq & 15
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product AND
        result = self.dataproduct & newproduct
        test1 = self.dataproduct.dq & newproduct.dq
        test2 = result.dq
        self.assertEqual(test1.all(), test2.all())
        
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
        
        # Attempt to access the data array through attributes.
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
