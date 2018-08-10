#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_distortion_model - Contains the unit tests for the classes
in the datamodels.miri_distortion_model module.

:History:

22 Feb 2013: Created.
29 Jul 2013: stats() method added.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
12 Sep 2013: Added more unit tests for creation and copy.
13 Sep 2013: Changed the (rmatrix,cmatrix) entries to (amatrix,bmatrix,
             tmatrix,mmatrix).
16 Sep 2013: Removed the ORDER parameter, since it is derivable from
             the size of the matrices.
30 Oct 2013: All MIRI distortion models (imaging, LRS, MRS) combined
             into one module. New model for LRS distortion and wavelength
             calibration.
31 Oct 2013: BETA array removed from MRS D2C model.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: TYPE and REFTYPE are no longer identical.
11 Mar 2015: group_integration_time changed to group_time.
20 May 2015: Corrected a mistake where tests were being run on a
             TestMiriMrsD2CModel object with null content.
02 Jul 2015: Major change to MRS distortion model. MiriMrsD2CModel (data array
             containing looking tables) replaced by MiriMrsDistortionModel
             (tables containing polynomial coefficients).
07 Oct 2015: Made exception catching Python 3 compatible.
10 Dec 2015: Old TestMiriMrsD2CModel tests removed.
15 Jun 2017: Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
10 Aug 2018: Updated MRS distortion models to reflect CDP-7 format.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_distortion_models import \
    MiriImagingDistortionModel, MiriLrsD2WModel, MiriMrsDistortionModel12, \
    MiriMrsDistortionModel34
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal


class TestMiriImagingDistortionModel(unittest.TestCase):
    
    # Test the MiriImagingDistortionModel class.
    
    def setUp(self):
        # Create a typical distortion product.
        self.bmatrix = [[0.1,0.2,0.3,0.4],
               [0.5,0.6,0.7,0.8],
               [0.8,0.7,0.6,0.5],
               [0.4,0.3,0.2,0.1]
               ]
        self.amatrix = [[0.1,0.0,0.0,0.0],
               [0.0,0.1,0.0,0.0],
               [0.0,0.0,0.1,0.0],
               [0.0,0.0,0.0,0.1]
               ]
        self.tmatrix = [[0.5,0.0,0.0],
               [0.2,0.5,0.2],
               [0.0,0.0,0.5]
               ]
        self.mmatrix = [[0.0,0.0,0.1],
               [0.0,0.1,0.0],
               [0.1,0.0,0.0]
               ]
        self.dataproduct = MiriImagingDistortionModel( amatrix=self.bmatrix,
                                                bmatrix=self.amatrix,
                                                tmatrix=self.tmatrix,
                                                mmatrix=self.mmatrix )
        # Add some example metadata.
        self.dataproduct.set_target_metadata(0.0, 0.0)
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
        self.testfile = "MiriImagingDistortionModel_test.fits"
        
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
        # Check the successful creation of a data product with all 4 matrices.
        dataproduct1 = MiriImagingDistortionModel( bmatrix=self.bmatrix,
                                            amatrix=self.amatrix,
                                            tmatrix=self.tmatrix,
                                            mmatrix=self.mmatrix )
        self.assertTrue( np.allclose(dataproduct1.bmatrix, self.bmatrix) )
        self.assertTrue( np.allclose(dataproduct1.amatrix, self.amatrix) )
        self.assertTrue( np.allclose(dataproduct1.tmatrix, self.tmatrix) )
        self.assertTrue( np.allclose(dataproduct1.mmatrix, self.mmatrix) )
        descr = str(dataproduct1)
        self.assertIsNotNone(descr)
        del dataproduct1

        # It should be possible to create a model containing A and B matrices
        # only.
        dataproduct2 = MiriImagingDistortionModel( bmatrix=self.bmatrix,
                                            amatrix=self.amatrix )
        self.assertTrue( np.allclose(dataproduct2.bmatrix, self.bmatrix) )
        self.assertTrue( np.allclose(dataproduct2.amatrix, self.amatrix) )
        descr = str(dataproduct2)
        self.assertIsNotNone(descr)
        del dataproduct2
        
        # Attempting to create a data model with arrays of the wrong shape
        # should fail.
        amatrix1D = [0.1,0.2,0.3,0.4]
        bmatrix1D = [0.1,0.0,0.0,0.0]
        self.assertRaises(TypeError, MiriImagingDistortionModel, amatrix=amatrix1D,
                          bmatrix=bmatrix1D )

        # It should be possible to set up an empty data product with
        # a specified shape. Both arrays should be initialised to
        # the same shape.
        emptydp = MiriImagingDistortionModel( (4,4) )
        self.assertIsNotNone(emptydp.amatrix)
        self.assertEqual(emptydp.amatrix.shape, (4,4))
        self.assertIsNotNone(emptydp.bmatrix)
        self.assertEqual(emptydp.bmatrix.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr
        
        # A null data product can also be created and populated
        # with data later.
        nulldp = MiriImagingDistortionModel( )
        descr1 = str(nulldp)
        nulldp.amatrix = np.asarray(self.amatrix)
        self.assertIsNotNone(nulldp.amatrix)
        nulldp.bmatrix = np.asarray(self.bmatrix)
        self.assertIsNotNone(nulldp.bmatrix)
        descr2 = str(nulldp)
        del nulldp, descr1, descr2
            
    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        self.assertTrue( np.allclose(self.dataproduct.amatrix,
                                     datacopy.amatrix) )
        self.assertTrue( np.allclose(self.dataproduct.bmatrix,
                                     datacopy.bmatrix) )
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriImagingDistortionModel(self.testfile) as readback:
                self.assertTrue( np.allclose(self.dataproduct.amatrix,
                                             readback.amatrix) )
                self.assertTrue( np.allclose(self.dataproduct.bmatrix,
                                             readback.bmatrix) )
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


class TestMiriLrsD2WModel(unittest.TestCase):
    
    # Test the MiriLrsD2WModel class.
    
    def setUp(self):
        # Create a typical wavelength_calibration product.
        wavedata = [
          (61.13976, 80.25328, 5.12652, 78.78318, 80.00265, 43.53634, 81.47993, 43.49634, 80.50392, 78.74317, 79.02664),
          (61.09973, 79.27727, 5.15676, 78.74317, 79.02664, 43.49634, 80.50392, 43.45628, 79.52791, 78.70311, 78.05063),
          (61.05965, 78.30126, 5.18700, 78.70311, 78.05063, 43.45628, 79.52791, 43.41617, 78.55190, 78.66300, 77.07462),
          (61.01951, 77.32526, 5.21725, 78.66300, 77.07462, 43.41617, 78.55190, 43.37601, 77.57589, 78.62285, 76.09861),
          (60.97933, 76.34925, 5.24749, 78.62285, 76.09861, 43.37601, 77.57589, 43.33580, 76.59988, 78.58264, 75.12260),
          (60.93910, 75.37324, 5.27773, 78.58264, 75.12260, 43.33580, 76.59988, 43.29554, 75.62387, 78.54238, 74.14659),
          (60.89881, 74.39723, 5.30797, 78.54238, 74.14659, 43.29554, 75.62387, 43.25523, 74.64786, 78.50207, 73.17058),
          (60.85848, 73.42122, 5.33821, 78.50207, 73.17058, 43.25523, 74.64786, 43.21487, 73.67186, 78.46171, 72.19458)
          ]
        self.dataproduct = MiriLrsD2WModel( wavelength_table=wavedata )
        # Add some example metadata.
        self.dataproduct.set_target_metadata(0.0, 0.0)
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
        self.testfile = "MiriLrsD2WModel_test.fits"
        
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
        class_names = list(MiriLrsD2WModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('wavelength_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.wavelength_table)
        self.assertEqual( len(self.dataproduct.wavelength_table),
                          len(datacopy.wavelength_table) )
        table1 = np.asarray(self.dataproduct.wavelength_table)
        table2 = np.asarray(datacopy.wavelength_table)
        assert_recarray_equal(table1, table2)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsD2WModel(self.testfile) as readback:
                self.assertIsNotNone(readback.wavelength_table)
                self.assertEqual( len(self.dataproduct.wavelength_table),
                                  len(readback.wavelength_table) )
                original = np.asarray(self.dataproduct.wavelength_table)
                duplicate = np.asarray(readback.wavelength_table)
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


class TestMiriMrsDistortionModel12(unittest.TestCase):
    
    # Test the MiriMrsDistortionModel class.
    
    def setUp(self):
        # Create a typical mrs_d2c product.
        slicenumber = [[1,2,3,4],
                       [1,2,3,4],
                       [1,2,3,4],
                       [1,2,3,4]
                       ]
        slicenumber3 = [slicenumber, slicenumber]
        fovdata = [(-2.95, 3.09),
                   (-2.96, 3.00)]
        d2cdata = [(100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0), 
                   (101.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0)]
        c2ddata = [(99.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0), 
                   (98.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    1.0, 22.0, 323.0, 24.0, 25.0)]
        transform = [('T_CH3C,V2', 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91),
                     ('T_CH3C,V3', 0.12, 0.22, 0.32, 0.42, 0.52, 0.62, 0.72, 0.82, 0.92)]
     
        self.dataproduct = MiriMrsDistortionModel12( slicenumber=slicenumber3,
                                 fov_ch1=fovdata, fov_ch2=fovdata,
                                 alpha_ch1=d2cdata, lambda_ch1=d2cdata,
                                 alpha_ch2=d2cdata, lambda_ch2=d2cdata,
                                 x_ch1=c2ddata, y_ch1=c2ddata,
                                 x_ch2=c2ddata, y_ch2=c2ddata,
                                 albe_xanyan=transform, xanyan_albe=transform,
                                 bzero1=-1.772, bdel1=0.177,
                                 bzero2=-2.238, bdel2=0.280)
        # Add some example metadata.
        self.dataproduct.set_target_metadata(0.0, 0.0)
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
        self.testfile = "MiriDistortionModel_test.fits"
        
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
        # It should be possible to set up an empty data product with
        # a specified shape. Both arrays should be initialised to
        # the same shape.
        emptydp = MiriMrsDistortionModel12( (4,4) )
        self.assertIsNotNone(emptydp.slicenumber)
        self.assertEqual(emptydp.slicenumber.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        # TODO: *** ADD MORE TESTS
        del emptydp, descr    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        self.assertTrue( np.allclose(self.dataproduct.slicenumber,
                                     datacopy.slicenumber) )
        table1 = np.asarray(self.dataproduct.fov_ch1)
        table2 = np.asarray(datacopy.fov_ch1)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.fov_ch2)
        table2 = np.asarray(datacopy.fov_ch2)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.alpha_ch1)
        table2 = np.asarray(datacopy.alpha_ch1)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.lambda_ch1)
        table2 = np.asarray(datacopy.lambda_ch1)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.alpha_ch2)
        table2 = np.asarray(datacopy.alpha_ch2)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.lambda_ch2)
        table2 = np.asarray(datacopy.lambda_ch2)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.x_ch1)
        table2 = np.asarray(datacopy.x_ch1)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.y_ch1)
        table2 = np.asarray(datacopy.y_ch1)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.x_ch2)
        table2 = np.asarray(datacopy.x_ch2)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.y_ch2)
        table2 = np.asarray(datacopy.y_ch2)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.albe_to_xanyan)
        table2 = np.asarray(datacopy.albe_to_xanyan)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.xanyan_to_albe)
        table2 = np.asarray(datacopy.xanyan_to_albe)
        assert_recarray_equal(table1, table2)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
 
            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsDistortionModel12(self.testfile) as readback:
                self.assertTrue( np.allclose(self.dataproduct.slicenumber,
                                             readback.slicenumber) )
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['slicenumber'],
                                       tables=['fov_ch1', 'fov_ch2',
                                               'alpha_ch1', 'lambda_ch1',
                                               'alpha_ch2', 'lambda_ch2',
                                               'x_ch1', 'y_ch1',
                                               'x_ch2', 'y_ch2',
                                               'albe_to_xanyan',
                                               'xanyan_to_albe'] )
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

class TestMiriMrsDistortionModel34(unittest.TestCase):
    
    # Test the MiriMrsDistortionModel class.
    
    def setUp(self):
        # Create a typical mrs_d2c product.
        slicenumber = [[1,2,3,4],
                       [1,2,3,4],
                       [1,2,3,4],
                       [1,2,3,4]
                       ]
        slicenumber3 = [slicenumber, slicenumber]
        fovdata = [(-2.95, 3.09),
                   (-2.96, 3.00)]
        d2cdata = [(100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0), 
                   (101.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0)]
        c2ddata = [(99.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    21.0, 22.0, 323.0, 24.0, 25.0), 
                   (98.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                    1.0, 22.0, 323.0, 24.0, 25.0)]
        transform = [('T_CH3C,V2', 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91),
                     ('T_CH3C,V3', 0.12, 0.22, 0.32, 0.42, 0.52, 0.62, 0.72, 0.82, 0.92)]
     
        self.dataproduct = MiriMrsDistortionModel34( slicenumber=slicenumber3,
                                 fov_ch3=fovdata, fov_ch4=fovdata,
                                 alpha_ch3=d2cdata, lambda_ch3=d2cdata,
                                 alpha_ch4=d2cdata, lambda_ch4=d2cdata,
                                 x_ch3=c2ddata, y_ch3=c2ddata,
                                 x_ch4=c2ddata, y_ch4=c2ddata,
                                 albe_xanyan=transform, xanyan_albe=transform,
                                 bzero3=-1.772, bdel3=0.177,
                                 bzero4=-2.238, bdel4=0.280)
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
        self.testfile = "MiriDistortionModel_test.fits"
        
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
        # It should be possible to set up an empty data product with
        # a specified shape. Both arrays should be initialised to
        # the same shape.
        emptydp = MiriMrsDistortionModel34( (4,4) )
        self.assertIsNotNone(emptydp.slicenumber)
        self.assertEqual(emptydp.slicenumber.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        # TODO: *** ADD MORE TESTS
        del emptydp, descr    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        self.assertTrue( np.allclose(self.dataproduct.slicenumber,
                                     datacopy.slicenumber) )
        table1 = np.asarray(self.dataproduct.fov_ch3)
        table2 = np.asarray(datacopy.fov_ch3)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.fov_ch4)
        table2 = np.asarray(datacopy.fov_ch4)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.alpha_ch3)
        table2 = np.asarray(datacopy.alpha_ch3)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.lambda_ch3)
        table2 = np.asarray(datacopy.lambda_ch3)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.alpha_ch4)
        table2 = np.asarray(datacopy.alpha_ch4)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.lambda_ch4)
        table2 = np.asarray(datacopy.lambda_ch4)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.x_ch3)
        table2 = np.asarray(datacopy.x_ch3)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.y_ch3)
        table2 = np.asarray(datacopy.y_ch3)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.x_ch4)
        table2 = np.asarray(datacopy.x_ch4)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.y_ch4)
        table2 = np.asarray(datacopy.y_ch4)
        assert_recarray_equal(table1, table2)
        
        table1 = np.asarray(self.dataproduct.albe_to_xanyan)
        table2 = np.asarray(datacopy.albe_to_xanyan)
        assert_recarray_equal(table1, table2)
        table1 = np.asarray(self.dataproduct.xanyan_to_albe)
        table2 = np.asarray(datacopy.xanyan_to_albe)
        assert_recarray_equal(table1, table2)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
 
            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsDistortionModel34(self.testfile) as readback:
                self.assertTrue( np.allclose(self.dataproduct.slicenumber,
                                             readback.slicenumber) )
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['slicenumber'],
                                       tables=['fov_ch3', 'fov_ch4',
                                               'alpha_ch3', 'lambda_ch3',
                                               'alpha_ch4', 'lambda_ch4',
                                               'x_ch3', 'y_ch3',
                                               'x_ch4', 'y_ch4',
                                               'albe_to_xanyan',
                                               'xanyan_to_albe'] )
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


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
