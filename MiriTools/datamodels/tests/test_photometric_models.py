#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_fluxconversion_model - Contains the unit tests for the classes
in the datamodels.miri_fluxconversion_model module.

:History:

10 Sep 2015: Created from test_fluxconversion_models.
07 Oct 2015: Made exception catching Python 3 compatible.
17 Nov 2015: SUBARRAY is GENERIC when photometric factors do not depend
             on subarray.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_photometric_models import \
    MiriPhotometricModel,  MiriImagingPhotometricModel, MiriPixelAreaModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

MAX_NLEM = 500 # Maximum size of the wavelength and relresponse arrays


class TestMiriPhotometricModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a typical photometric data product.
        # This product is a combination of imager and LRS, which should
        # exercise the flexibility of the data model.
        pixar_a2 = 0.136
        wavelength = []
        relresponse = []
        resp = 0
        nelm = 0
        for ii in range(300):
            wav = float(ii)/12.0
            resp = (resp + 1) % 30
            r10 = 0.1 + resp/35.0
            wavelength.append(wav)
            relresponse.append(r10)
            nelm += 1
        for ii in range(300,MAX_NLEM):
            wavelength.append(0.0)
            relresponse.append(0.0)
        self.phot_table = \
              [('F560W',  'GENERIC', 2.41,  0.26,  0, wavelength, relresponse),
               ('F770W',  'GENERIC', 1.32,  0.013, 0, wavelength, relresponse),
               ('F1000W', 'GENERIC', 1.76,  0.12,  0, wavelength, relresponse),
               ('F1130W', 'GENERIC', 5.76,  0.43,  0, wavelength, relresponse),
               ('F1280W', 'GENERIC', 2.11,  0.16,  0, wavelength, relresponse),
               ('F1500W', 'GENERIC', 1.84,  0.01,  0, wavelength, relresponse),
               ('F1800W', 'GENERIC', 2.68,  0.23,  0, wavelength, relresponse),
               ('F2100W', 'GENERIC', 2.04,  0.15,  0, wavelength, relresponse),
               ('F2550W', 'GENERIC', 4.25,  0.4,   0, wavelength, relresponse),
               ('F2550WR','GENERIC', 4.60,  0.24,  0, wavelength, relresponse),
               ('F1065C', 'GENERIC', 1.37,  0.1,   0, wavelength, relresponse),
               ('F1140C', 'GENERIC', 1.43,  0.11,  0, wavelength, relresponse),
               ('F1550C', 'GENERIC', 1.81,  0.13,  0, wavelength, relresponse),
               ('F2300C', 'GENERIC', 3.65,  0.23,  0, wavelength, relresponse),
               ('P750L',  'GENERIC',       1.0,  0.0,  nelm, wavelength, relresponse),
               ('P750L',  'SLITLESSPRISM', 0.9,  0.0,  nelm, wavelength, relresponse)
               ]

        self.dataproduct = MiriPhotometricModel( phot_table=self.phot_table, 
                                                 pixar_a2=pixar_a2 )
        self.testfile = "MiriPhotometric_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.phot_table
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
        class_names = list(MiriPhotometricModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('phot_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriPhotometricModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.phot_table = self.phot_table
        self.assertIsNotNone(nulldp.phot_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.phot_table)
        self.assertEqual( len(self.dataproduct.phot_table),
                          len(datacopy.phot_table) )
        table1 = np.asarray(self.dataproduct.phot_table)
        table2 = np.asarray(datacopy.phot_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriPhotometricModel(self.testfile) as readback:
                self.assertIsNotNone(readback.phot_table)
                self.assertEqual( len(self.dataproduct.phot_table),
                                  len(readback.phot_table) )
                original = np.asarray(self.dataproduct.phot_table)
                duplicate = np.asarray(readback.phot_table)
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
        descr = str(self.dataproduct.phot_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriImagingPhotometricModel(unittest.TestCase):
    
    # Test the MiriImagingFluxconversionModel class.
        
    def setUp(self):
        # Create a typical imaging data product.
        # The wavelength and relresponse arrays will be automatically zeroed.
        pixar_a2 = 0.136
        self.phot_table = \
             [('F560W',   '', 2.41,  0.26),
               ('F770W',  '', 1.32,  0.013),
               ('F1000W', '', 1.76,  0.12),
               ('F1130W', '', 5.76,  0.43),
               ('F1280W', '', 2.11,  0.16),
               ('F1500W', '', 1.84,  0.01),
               ('F1800W', '', 2.68,  0.23),
               ('F2100W', '', 2.04,  0.15),
               ('F2550W', '', 4.25,  0.4),
               ('F2550WR','', 4.60,  0.24),
               ('F1065C', '', 1.37,  0.1),
               ('F1140C', '', 1.43,  0.11),
               ('F1550C', '', 1.81,  0.13),
               ('F2300C', '', 3.65,  0.23)
               ]

        self.dataproduct = MiriImagingPhotometricModel( phot_table=self.phot_table, 
                                                 pixar_a2=pixar_a2 )
        self.testfile = "MiriImagingPhotometric_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.phot_table
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
        class_names = list(MiriPhotometricModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('phot_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriPhotometricModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        del nulldp, descr1 

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.phot_table)
        self.assertEqual( len(self.dataproduct.phot_table),
                          len(datacopy.phot_table) )
        table1 = np.asarray(self.dataproduct.phot_table)
        table2 = np.asarray(datacopy.phot_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriPhotometricModel(self.testfile) as readback:
                self.assertIsNotNone(readback.phot_table)
                self.assertEqual( len(self.dataproduct.phot_table),
                                  len(readback.phot_table) )
                original = np.asarray(self.dataproduct.phot_table)
                duplicate = np.asarray(readback.phot_table)
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
        descr = str(self.dataproduct.phot_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriPixelAreaModel(unittest.TestCase):
    
    # Test the MiriPixelAreaModel class.
    
    def setUp(self):
        # Create a typical pixel area data product.
        data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
        self.dataproduct = MiriPixelAreaModel( data=data3x3, pixar_a2=0.13 )
        self.testfile = "MiriPixelAreaModel_test.fits"
        
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
        # It must be possible to create an empty data product.    
        nullproduct = MiriPixelAreaModel( )
        del nullproduct
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriPixelAreaModel(self.testfile) as readback:
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


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
