#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_linearity_model - Contains the unit tests for the classes
in the datamodels.miri_linearity_model module.

:History:

22 Feb 2013: Created.
14 Jun 2013: Corrected typo in temporary file name.
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
04 Nov 2014: Removed FITERR, MIN and MAX, to bring the unit test up to date
             with changes made on 29 Oct 2014.
12 Jan 2015: When testing that the model can be created with an empty
             array, the array shape must still have the correct number
             of dimensions.
11 Mar 2015: group_integration_time changed to group_time.
20 Aug 2015: Name changed from MiriNonlinearityModel to MiriLinearityModel,
             for consistency with STScI naming. minimage, maximage and fiterr
             arrays removed.
07 Oct 2015: Made exception catching Python 3 compatible.
15 Jun 2017: Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_linearity_model import \
    linearity_reference_flags, MiriLinearityModel
from miri.datamodels.tests.util import assert_products_equal


class TestMiriLinearityModel(unittest.TestCase):
    
    # Test the MiriLinearityModel class.
    
    def setUp(self):
        # Create a typical linearity product.
        data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
        err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
        dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
        data3x3x2 = [data3x3,data3x3]
        self.dataproduct = MiriLinearityModel( coeffs=data3x3x2, err=err3x3,
                                dq=dq3x3, dq_def=linearity_reference_flags )
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
        self.testfile = "MiriLinearityModel_test.fits"
        
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
        data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
        data3x3x2 = [data3x3,data3x3]
        # The coeffs and data attributes should point to the same data array.
        testdp = MiriLinearityModel( coeffs=data3x3x2 )
        self.assertIsNotNone(testdp.coeffs)
        self.assertIsNotNone(testdp.data)
        self.assertTrue(testdp.coeffs.any() == testdp.data.any())
        del testdp
        
        # It should be possible to set up an empty data product with
        # a specified shape. All three arrays should be initialised to
        # the same shape.
        emptydp = MiriLinearityModel( (4,4,2) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (4,4,2))
        self.assertIsNotNone(emptydp.err)
        self.assertEqual(emptydp.err.shape, (4,4,2))
        self.assertIsNotNone(emptydp.dq)
        self.assertEqual(emptydp.dq.shape, (4,4,2))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr
        
        # Attempting to create an array with the wrong number of
        # dimensions should cause an exception.
        self.assertRaises(ValueError, MiriLinearityModel, coeffs=data3x3)
        del data3x3, data3x3x2

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq'],
#                                        'minimage', 'maximage', 'fiterr'],
                               tables='dq_def' )
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLinearityModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
#                                                'minimage', 'maximage', 'fiterr'],
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
        # Check that the coeffs alias works
        data = self.dataproduct.data
        coeffs = self.dataproduct.coeffs
        self.assertTrue(data.all() == coeffs.all())


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
