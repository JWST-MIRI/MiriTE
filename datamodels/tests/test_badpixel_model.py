#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_badpixel_model - Contains the unit tests for the classes
in the datamodels.miri_badpixel_model module.

:History:

16 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
23 Apr 2013: mask element renamed data to keep up with changes
             in jwst_lib data model (where inclusion of a data
             element is now compulsory).
29 Jul 2013: stats() method added.
22 Aug 2013: columnnames renamed to dq_def_names. Check that the field
             names declared in the class variable match the schema.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
             Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Test that the data product can be copied successfully.
04 Oct 2013: Standard bad pixel flags declared with MiriBadPixelMask.
26 Jun 2014: field_def changed to dq_def.
29 Aug 2014: Use the JWST standard mask reference flags.
             Added test_referencefile.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
18 Aug 2015: Mask data are now stored in a "dq" array rather than a "data"
             array.
02 Sep 2015: Added subarray extraction function.
07 Oct 2015: Made exception catching Python 3 compatible.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_badpixel_model import \
    mask_reference_flags, MiriBadPixelMaskModel
from miri.datamodels.tests.util import assert_products_equal


class TestMiriBadPixelMaskModel(unittest.TestCase):
    
    # Test the MiriBadPixelMaskModel class.
        
    def setUp(self):       
        # Create a typical bad pixel mask product.
        self.maskdata = np.array([[0,1,0],[1,0,1],[0,1,0]])
        self.dataproduct = MiriBadPixelMaskModel( dq=self.maskdata,
                                                  dq_def=mask_reference_flags )
        self.testfile = "MiriBadPixelMask_test.fits"
        
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
        # Create a bad pixel product from a list of tuples.
        testmask = MiriBadPixelMaskModel( dq=self.maskdata,
                                          dq_def=mask_reference_flags )
        # Check that the field names in the class variable are the same
        # as the ones declared in the schema.
        class_names = list(MiriBadPixelMaskModel.dq_def_names)
        schema_names = list(testmask.get_field_names('dq_def'))
        self.assertEqual(class_names, schema_names,
                         "'dq_def_names' class variable does not match schema")
        # The product should contain a valid mask.
        self.assertIsNotNone(testmask.dq)
        descr = str(testmask)
        self.assertIsNotNone(descr)
        del testmask, descr

        # A mask with no field values should also be valid.
        testmask = MiriBadPixelMaskModel( dq=self.maskdata, dq_def=None )
        self.assertIsNotNone(testmask.dq)
        descr = str(testmask)
        del testmask, descr
        # An empty mask with a specified size is also valid.
        testmask = MiriBadPixelMaskModel( (3,3) )
        descr = str(testmask)
        del testmask, descr

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays='dq', tables='dq_def' )
        del datacopy

    def test_subarray(self):
        # Test that subarrays can be extracted from the data product.
        subdata1 = self.dataproduct.get_subarray( (1,1,2,2) )
        del subdata1
        subdata2 = self.dataproduct.get_subarray( (2,2,2,2) )
        del subdata2
        subdata3 = self.dataproduct.get_subarray( (1,1,3,3) )
        del subdata3
        subdata4 = self.dataproduct.get_subarray( (3,3,1,1) )
        del subdata4
        
        # Attempting to extract a subarray larger than the mask
        # should raise an exception.
        self.assertRaises(IndexError, self.dataproduct.get_subarray, (0,0,4,4) )
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriBadPixelMaskModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays='dq', tables='dq_def' )
                del readback

    def test_operators(self):
        newmask = np.array([[1,1,0],[0,1,1],[0,0,1]])
        newproduct = MiriBadPixelMaskModel( dq=newmask,
                                            dq_def=mask_reference_flags )

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
        
        # Attempt to access the mask array and dq_def table through attributes.
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq_def)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
