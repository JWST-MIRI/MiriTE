#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_illumination_model - Contains the unit tests for the classes
in the datamodels.miri_illumination_model module.

:History:

22 Aug 2013: Created.
13 Aug 2015: Added TestMiriIlluminationFringingModel
07 Oct 2015: Made exception catching Python 3 compatible.
18 Nov 2015: Added truncate() method.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_illumination_model import \
    MiriIlluminationModel, MiriIlluminationFringingModel


class TestMiriIlluminationModel(unittest.TestCase):
    
    # Test the MiriIlluminationModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriIlluminationModel are
    # included here.
    
    def setUp(self):
        # Create a typical illumination product.
        ii1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        ii3 = [ii1,ii1,ii1]
        ww = [[1,2,3,4], [5,6,7,8], [9,10,11,12]]
        ww3 = [ww,ww,ww]
        self.dataproduct = MiriIlluminationModel(intensity=ii3, wavelength=ww3)
        self.testfile = "MiriIlluminationModel_test.fits"
        
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

    def test_creation(self):
        ii1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        ii3 = [ii1,ii1,ii1]
        ww = [[1,2,3,4], [5,6,7,8], [9,10,11,12]]
        ww3 = [ww,ww,ww]
        # 1) Intensity array only. Intensity array must exist and be non-empty.
        # Wavelength array should exist and be the same size and shape as the
        # intensity array.
        newdp1 = MiriIlluminationModel(intensity=ii3)
        self.assertIsNotNone(newdp1.intensity)
        self.assertGreater(len(newdp1.intensity), 0)
        self.assertIsNotNone(newdp1.wavelength)
        self.assertEqual(newdp1.wavelength.shape, newdp1.wavelength.shape)
        del newdp1

        # 2) 3-D intensity array and 2-D wavelength array should work
        newdp2 = MiriIlluminationModel(intensity=ii3, wavelength=ww)
        self.assertIsNotNone(newdp2.intensity)
        self.assertIsNotNone(newdp2.wavelength)
        del newdp2

        # 3) 3-D intensity array and 1-D wavelength array should also work
        ww1 = [5,10,15]
        newdp3 = MiriIlluminationModel(intensity=ii3, wavelength=ww1)
        self.assertIsNotNone(newdp3.intensity)
        self.assertIsNotNone(newdp3.wavelength)
        del newdp3

        # 4) The creation should fail if the wavelength array has
        # an incompatible shape.
        wwbad = [ww,ww,ww,ww]
        newdp4 = self.assertRaises(TypeError, MiriIlluminationModel,
                                   intensity=ii3, wavelength=wwbad)
        del newdp4
        wwbad = [5, 10, 15, 20]
        newdp5 = self.assertRaises(TypeError, MiriIlluminationModel,
                                   intensity=ii3, wavelength=wwbad)
        del newdp5

        # 5) It must be possible to create an empty object
        nulldp = MiriIlluminationModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        
    def test_reshaping(self):
        # Query the illumination shape and then truncate the illumination
        # model to half that size.
        illum_shape = self.dataproduct.get_illumination_shape()
        self.assertTrue(isinstance(illum_shape, (tuple,list)))
        self.assertEqual(len(illum_shape), 2)
        
        new_shape = []
        new_shape.append( illum_shape[0] // 2 )
        new_shape.append( illum_shape[1] // 2 )
        self.dataproduct.truncate( new_shape )
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriIlluminationModel(self.testfile) as readback:
                self.assertTrue( np.allclose(self.dataproduct.intensity,
                                             readback.intensity) )
                self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                             readback.wavelength) )
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
        
        # Attempt to access the intensity and wavelength arrays through attributes.
        descr = str(self.dataproduct.intensity)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.wavelength)
        self.assertIsNotNone(descr)
        del descr


class TestMiriIlluminationFringingModel(unittest.TestCase):
    
    # Test the MiriIlluminationModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriIlluminationModel are
    # included here.
    
    def setUp(self):
        # Create a typical illumination product.
        ii1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        ii3 = [ii1,ii1,ii1]
        ww = [[1,2,3,4], [5,6,7,8], [9,10,11,12]]
        ww3 = [ww,ww,ww]
        dd1 = [[45,45,45,45], [45,45,45,45], [45,45,45,45]]
        dd2 = [[135,135,135,135], [135,135,135,135], [135,135,135,135]]
        dd3 = [dd1, dd2]
        dd4 = [dd3,dd3,dd3]
        self.dataproduct = MiriIlluminationFringingModel(intensity=ii3,
                                                         wavelength=ww3,
                                                         direction=dd4)
        self.testfile = "MiriIlluminationFringingModel_test.fits"
        
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

    def test_creation(self):
        # TBD
        pass
    
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriIlluminationFringingModel(self.testfile) as readback:
                self.assertTrue( np.allclose(self.dataproduct.intensity,
                                             readback.intensity) )
                self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                             readback.wavelength) )
                self.assertTrue( np.allclose(self.dataproduct.direction,
                                             readback.direction) )
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
        
        # Attempt to access the intensity and wavelength arrays through attributes.
        descr = str(self.dataproduct.intensity)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.wavelength)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.direction)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
