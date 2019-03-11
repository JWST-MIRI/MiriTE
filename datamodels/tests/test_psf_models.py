#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_psf_model - Contains the unit tests for the classes
in the datamodels.miri_psf_model module.

:History:

16 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
13 Feb 2013: Added psf_lut table.
26 Apr 2013: File closing problem has returned!
29 Jul 2013: stats() method added.
22 Aug 2013: Check that the field names declared in the class variable
             match the schema.
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Test that the data product can be copied successfully.
10 Oct 2013: psf_lut table removed from MiriPointSpreadFunctionModel. It now
             only applies to MiriImagingPointSpreadFunctionModel.
30 Oct 2013: Separate PSF models for imaging, LRS and MRS.
08 Apr 2014: Brought up to date with changes made on 28 Mar 2014.
             Wavelength extension removed from MiriLrsPointSpreadFunctionModel.
             Changed type back to 'PSF'.
09 Jul 2014: field_def changed to dq_def.
21 Jul 2014: IM, LW and SW detectors changed to MIRIMAGE, MIRIFULONG and
             MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
31 Oct 2014: Updated to match the schema changes made on 30 Oct 2014.
16 Jan 2015: Included extra creation tests.
11 Mar 2015: group_integration_time changed to group_time.
07 Oct 2015: Made exception catching Python 3 compatible.
03 Dec 2015: Added 4 columns to the imager PSF_LUT table: COL_FIELD, ROW_FIELD,
             XAN_FIELD, YAN_FIELD.
11 Dec 2015: Test the PSF-OOF and PSF-MONOCHROM variants.
15 Jun 2017: Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
27 Apr 2018: Temporarily comment out the TypeError test while converting
             to Python 3.
30 Jan 2019: Test that the REFTYPE and DATAMODL metadata is not altered
             when the data model is saved to a file.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_psf_models import psf_reference_flags, \
    MiriPointSpreadFunctionModel, MiriImagingPointSpreadFunctionModel, \
    MiriLrsPointSpreadFunctionModel, MiriMrsPointSpreadFunctionModel
from miri.datamodels.tests.util import assert_products_equal


class TestMiriPointSpreadFunctionModel(unittest.TestCase):
    
    # Test the MiriPointSpreadFunctionModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriPointSpreadFunctionModel are
    # included here.
    
    def setUp(self):
        # Create a typical 2-D PSF product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        dqdef = psf_reference_flags
        self.dataproduct = MiriPointSpreadFunctionModel(data=a1, err=b1, dq=c1,
                                                        dq_def=dqdef,
                                                        wavelen=10.0,
                                                        xfield=12.0, yfield=8.0)
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriPointSpreadFunctionModel_test.fits"
        
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
        emptydp = MiriPointSpreadFunctionModel( (4,4,2) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (4,4,2))
        self.assertIsNotNone(emptydp.err)
        self.assertEqual(emptydp.err.shape, (4,4,2))
        self.assertIsNotNone(emptydp.dq)
        self.assertEqual(emptydp.dq.shape, (4,4,2))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr
        
        # A null data product can also be created and populated
        # with data later.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        nulldp = MiriPointSpreadFunctionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.data = np.asarray(a1)
        self.assertIsNotNone(nulldp.err)
        self.assertIsNotNone(nulldp.dq)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2

        # Creating a data object with a data array or the wrong shape
        # should raise an exception.
        # NOTE: A bug in the JWST data model might cause an AttributeError
        # to be raised instead of a TypeError. If this happens, try a newer
        # version of the JWST data model library.
        a1d = [10,20,30,40]
        self.assertRaises(TypeError, MiriPointSpreadFunctionModel,
                          data=a1d)

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
            with MiriPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
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


class TestMiriImagingPointSpreadFunctionModel(unittest.TestCase):
    
    # Test the MiriImagingPointSpreadFunctionModel class.
    
    def setUp(self):
        # Create an imaging mode 3-D PSF product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        aa1 = [a1,a1,a1,a1]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        lut_im = [(1.0,  0.0, 0, 0.7, 0.7, 0.4, 0.4),
                  (1.1,  0.0, 1, 0.5, 0.5, 0.2, 0.2),
                  (1.0,  0.0, 2, 0.1, 0.1, 0.8, 0.8),
                  (1.1,  0.0, 3, 0.6, 0.6, 0.1, 0.1)
                  ]
        self.dataproduct = MiriImagingPointSpreadFunctionModel(data=aa1,
                                                               err=b1, dq=c1,
                                                               psf_lut=lut_im)
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
        self.testfile = "MiriImagingPointSpreadFunctionModel_test.fits"
        
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
        class_names = list(MiriImagingPointSpreadFunctionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('psf_lut'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # Try to create a data object with an invalid psf_lut should raise
        # an exception.
        # NOTE: A bug in the JWST data model might cause an AttributeError
        # to be raised instead of a TypeError. If this happens, try a newer
        # version of the JWST data model library.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        aa1 = [a1,a1,a1,a1]
        # FIXME: Causes a problem in Python 3 - exception not raised.
#         self.assertRaises(TypeError, MiriImagingPointSpreadFunctionModel,
#                           data=aa1, psf_lut=42)

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq'],
                               tables=['dq_def', 'psf_lut'] )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriImagingPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables=['dq_def', 'psf_lut'] )
                del readback

            # Check that other variations of the data product can be
            # written to a FITS file and read back again without changing
            # the data.
            a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
            aa1 = [a1,a1,a1,a1]
            b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
            c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
            lut_im = [(1.0,  0.0, 0, 0.7, 0.7, 0.4, 0.4),
                  (1.1,  0.0, 1, 0.5, 0.5, 0.2, 0.2),
                  (1.0,  0.0, 2, 0.1, 0.1, 0.8, 0.8),
                  (1.1,  0.0, 3, 0.6, 0.6, 0.1, 0.1)
                  ]
            testproduct = MiriImagingPointSpreadFunctionModel(data=aa1,
                                err=b1, dq=c1, psf_lut=lut_im,
                                psftype='PSF-OOF')
            testproduct.save(self.testfile, overwrite=True)
            with MiriImagingPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(testproduct.meta.reftype,
                                 readback.meta.reftype)
                assert_products_equal( self, testproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables=['dq_def', 'psf_lut'] )
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
        
        # Attempt to access the SCI, ERR and DQ arrays and PSF_LUT table
        # through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.err)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.psf_lut)
        self.assertIsNotNone(descr)
        del descr


class TestMiriLrsPointSpreadFunctionModel(unittest.TestCase):
    
    # Test the MiriLrsPointSpreadFunctionModel class.
    
    def setUp(self):
        # Create an LRS 3-D PSF product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        aa1 = [a1,a1,a1,a1]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.dataproduct = MiriLrsPointSpreadFunctionModel(data=aa1,
                                                           err=b1, dq=c1)
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
        self.testfile = "MiriMiriLrsPointSpreadFunctionModel_test.fits"
        
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
        # Try to create a data object with invalid data should raise
        # an exception. TBD.
        pass

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq'],
                               tables=['dq_def'] )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables=['dq_def',] )
                del readback

            # Check that other variations of the data product can be
            # written to a FITS file and read back again without changing
            # the data.
            a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
            aa1 = [a1,a1,a1,a1]
            b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
            c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
            testproduct = MiriLrsPointSpreadFunctionModel(data=aa1,
                                err=b1, dq=c1, psftype='PSF-MONOCHROM')
            testproduct.save(self.testfile, overwrite=True)
            with MiriImagingPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(testproduct.meta.reftype,
                                 readback.meta.reftype)
                assert_products_equal( self, testproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables=['dq_def',] )
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
        
        # Attempt to access the SCI, ERR and DQ arrays and WAVELENGTH array
        # through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        descr = str(self.dataproduct.err)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


class TestMiriMrsPointSpreadFunctionModel(unittest.TestCase):
    
    # Test the MiriMrsPointSpreadFunctionModel class.
    
    def setUp(self):
        # Create an MRS 3-D PSF product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        aa1 = [a1,a1,a1,a1]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.dataproduct = MiriMrsPointSpreadFunctionModel(data=aa1,
                                                           err=b1, dq=c1)
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFULONG',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriMiriMrsPointSpreadFunctionModel_test.fits"
        
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
        # Try to create a data object with invalid data should raise
        # an exception. TBD.
        pass

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq'],
                               tables=['dq_def'] )
        del datacopy
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsPointSpreadFunctionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables=['dq_def',] )
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
        
        # Attempt to access the SCI, ERR and DQ arrays
        # through attributes.
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
