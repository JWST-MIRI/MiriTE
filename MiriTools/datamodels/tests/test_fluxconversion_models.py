#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_fluxconversion_model - Contains the unit tests for the classes
in the datamodels.miri_fluxconversion_model module.

:History:

16 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
17 May 2013: Do not allow a blank table to be created.
22 Aug 2013: columnnames renamed to fieldnames. Check that the field
             names declared in the class variable match the schema.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst_lib
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
             Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Test that the data product can be copied successfully.
30 Oct 2013: LRS and MRS now use different flux conversion model classes.
09 Jul 2014: field_def changed to dq_def.
29 Aug 2014: Added test_referencefile.
07 Oct 2015: Made exception catching Python 3 compatible.
03 Dec 2015: Added MiriPowerlawColourCorrectionModel.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_fluxconversion_models import \
    MiriImagingFluxconversionModel,  MiriImagingColourCorrectionModel, \
    MiriPowerlawColourCorrectionModel, MiriLrsFluxconversionModel, \
    MiriMrsFluxconversionModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal


class TestMiriImagingFluxconversionModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a typical flux conversion product.
        
        self.flux = [('F560W',  1.0,  0.0),
                     ('F770W',  1.1,  0.0),
                     ('F1000W', 1.2,  0.01),
                     ('F1130W', 1.3,  0.0),
                     ('F1280W', 1.4,  0.0),
                     ('F1500W', 1.5,  0.02),
                     ('F1800W', 1.6,  0.0),
                     ('F2100W', 1.7,  0.03),
                     ('F2550W', 1.8,  0.0),
                     ]

        self.dataproduct = MiriImagingFluxconversionModel( \
                                flux_table=self.flux )
        self.testfile = "MiriImagingFluxconversion_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.flux
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
        class_names = list(MiriImagingFluxconversionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('flux_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriImagingFluxconversionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.flux_table = self.flux
        self.assertIsNotNone(nulldp.flux_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.flux_table)
        self.assertEqual( len(self.dataproduct.flux_table),
                          len(datacopy.flux_table) )
        table1 = np.asarray(self.dataproduct.flux_table)
        table2 = np.asarray(datacopy.flux_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriImagingFluxconversionModel(self.testfile) as readback:
                self.assertIsNotNone(readback.flux_table)
                self.assertEqual( len(self.dataproduct.flux_table),
                                  len(readback.flux_table) )
                original = np.asarray(self.dataproduct.flux_table)
                duplicate = np.asarray(readback.flux_table)
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
        
        # Attempt to access the flux table through attributes.
        descr = str(self.dataproduct.flux_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriImagingColourCorrectionModel(unittest.TestCase):
    
    # Test the MiriImagingColourCorrectionModel class.
        
    def setUp(self):
        # Create a typical imaging flux conversion product.
        self.flux = [(10.0, 'F560W',  1.0,  0.0),
                     (10.0, 'F770W',  1.1,  0.0),
                     (10.0, 'F1000W', 1.2,  0.01),
                     (10.0, 'F1130W', 1.3,  0.0),
                     (10.0, 'F1280W', 1.4,  0.0),
                     (10.0, 'F1500W', 1.5,  0.02),
                     (10.0, 'F1800W', 1.6,  0.0),
                     (10.0, 'F2100W', 1.7,  0.03),
                     (10.0, 'F2550W', 1.8,  0.0),
                     ]
        self.dataproduct = MiriImagingColourCorrectionModel( \
                                flux_table=self.flux )
        self.testfile = "MiriImagingColourCorrectionModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.flux
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
        class_names = list(MiriImagingColourCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('flux_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriImagingColourCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.flux_table = self.flux
        self.assertIsNotNone(nulldp.flux_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.flux_table)
        self.assertEqual( len(self.dataproduct.flux_table),
                          len(datacopy.flux_table) )
        table1 = np.asarray(self.dataproduct.flux_table)
        table2 = np.asarray(datacopy.flux_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriImagingColourCorrectionModel(self.testfile) as readback:
                self.assertIsNotNone(readback.flux_table)
                self.assertEqual( len(self.dataproduct.flux_table),
                                  len(readback.flux_table) )
                original = np.asarray(self.dataproduct.flux_table)
                duplicate = np.asarray(readback.flux_table)
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
        
        # Attempt to access the flux table through attributes.
        descr = str(self.dataproduct.flux_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriPowerlawColourCorrectionModel(unittest.TestCase):
    
    # Test the MiriImagingColourCorrectionModel class.
        
    def setUp(self):
        # Create a typical imaging flux conversion product.
        self.flux = [(1.1, 'F560W',  1.0,  0.0),
                     (1.1, 'F770W',  1.1,  0.0),
                     (1.2, 'F1000W', 1.2,  0.01),
                     (1.3, 'F1130W', 1.3,  0.0),
                     (1.4, 'F1280W', 1.4,  0.0),
                     (1.5, 'F1500W', 1.5,  0.02),
                     (1.6, 'F1800W', 1.6,  0.0),
                     (1.7, 'F2100W', 1.7,  0.03),
                     (1.8, 'F2550W', 1.8,  0.0),
                     ]
        self.dataproduct = MiriPowerlawColourCorrectionModel( \
                                flux_table=self.flux )
        self.testfile = "MiriPowerlawColourCorrectionModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.flux
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
        class_names = list(MiriPowerlawColourCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('flux_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriPowerlawColourCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.flux_table = self.flux
        self.assertIsNotNone(nulldp.flux_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.flux_table)
        self.assertEqual( len(self.dataproduct.flux_table),
                          len(datacopy.flux_table) )
        table1 = np.asarray(self.dataproduct.flux_table)
        table2 = np.asarray(datacopy.flux_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriPowerlawColourCorrectionModel(self.testfile) as readback:
                self.assertIsNotNone(readback.flux_table)
                self.assertEqual( len(self.dataproduct.flux_table),
                                  len(readback.flux_table) )
                original = np.asarray(self.dataproduct.flux_table)
                duplicate = np.asarray(readback.flux_table)
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
        
        # Attempt to access the flux table through attributes.
        descr = str(self.dataproduct.flux_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriLrsFluxconversionModel(unittest.TestCase):
    
    # Test the MiriLrsFluxconversionModel class.
        
    def setUp(self):
        # Create a typical LRS flux conversion product.
        self.flux = [( 2.0, 1.0,  0.0),
                     ( 4.0, 1.1,  0.0),
                     ( 6.0, 1.2,  0.01),
                     ( 8.0, 1.3,  0.0),
                     (10.0, 1.4,  0.0),
                     (12.0, 1.5,  0.02),
                     (14.0, 1.6,  0.0),
                     (16.0, 1.7,  0.03),
                     (18.0, 1.8,  0.0),
                     ]
        self.dataproduct = MiriLrsFluxconversionModel( \
                              flux_table=self.flux )
        self.testfile = "MiriLrsFluxconversion_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.flux
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
        class_names = list(MiriLrsFluxconversionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('flux_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLrsFluxconversionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.flux_table = self.flux
        self.assertIsNotNone(nulldp.flux_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.flux_table)
        self.assertEqual( len(self.dataproduct.flux_table),
                          len(datacopy.flux_table) )
        table1 = np.asarray(self.dataproduct.flux_table)
        table2 = np.asarray(datacopy.flux_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsFluxconversionModel(self.testfile) as readback:
                self.assertIsNotNone(readback.flux_table)
                self.assertEqual( len(self.dataproduct.flux_table),
                                  len(readback.flux_table) )
                original = np.asarray(self.dataproduct.flux_table)
                duplicate = np.asarray(readback.flux_table)
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
        
        # Attempt to access the flux table through attributes.
        descr = str(self.dataproduct.flux_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriMrsFluxconversionModel(unittest.TestCase):
    
    # Test the MiriMrsFluxconversionModel class.
        
    def setUp(self):
        # Create a typical MRS flux conversion product.        
        flux_plane = [(1.0, 1.1, 1.2),
                      (1.3, 1.4, 1.5),
                      (1.6, 1.7, 1.8)]
        self.flux = [flux_plane, flux_plane, flux_plane]
        err_plane = [(0.0, 0.01, 0.0),
                     (0.02, 0.0, 0.03),
                     (0.01, 0.04, 0.0)]
        self.err = [err_plane, err_plane, err_plane]
        self.dq =  [(1,0,0), (0,1,0), (1,0,1)]

        self.dataproduct = MiriMrsFluxconversionModel( \
                              data=self.flux, err=self.err, dq=self.dq )
        self.testfile = "MiriMrsFluxconversion_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.flux
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
        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsFluxconversionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.data = self.flux
        self.assertIsNotNone(nulldp.data)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.data)
        self.assertEqual( len(self.dataproduct.data),
                          len(datacopy.data) )
        flux1 = np.asarray(self.dataproduct.data)
        flux2 = np.asarray(datacopy.data)
        self.assertTrue( np.allclose(np.nan_to_num(flux1), np.nan_to_num(flux2)))
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsFluxconversionModel(self.testfile) as readback:
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
        
        # Attempt to access the flux data through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
