#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_aperture_correction_model - Contains the unit tests for the classes
in the datamodels.miri_aperture_correction_model module.

:History:

13 Jan 2020: Created 
22 Nov 2021: TestMiriLrsPathlossCorrectionModel added.


@author:  Juergen Schreiber (MPIA), Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_aperture_correction_model import \
    MiriMrsApertureCorrectionModel, MiriLrsApertureCorrectionModel,\
    MiriLrsThroughputCorrectionModel, MiriLrsPositionCorrectionModel, MiriLrsPathlossCorrectionModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriMrsApertureCorrectionModel(unittest.TestCase):
    
    # Test the MiriMrsApertureCorrectionModel class.
        
    def setUp(self):
        # Create a typical photometric data product.
        apercorrdata = [(5.0, 'nominal', 4.87,  7.76,  8.57, 1.0,  0.0, 42.0, 0.1),
                   (10.0, 'nominal', 5.87,  8.76,  9.57, 1.0,  0.0, 32.0, 0.1),
                   (15.0, 'nominal', 6.87,  9.76, 10.57, 1.0,  0.0, 32.0, 0.1),
                   (15.0, 'nominal', 7.87, 10.76, 11.57, 0.7, 45.0, 12.0, 0.1)]
            
        self.apercorr_table = apercorrdata

        self.dataproduct = MiriMrsApertureCorrectionModel( apercorr_table=self.apercorr_table)
        self.testfile = "MiriMrsApCorr_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.apercorr_table
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
        class_names = list(MiriMrsApertureCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('apercorr_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsApertureCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.apercorr_table = self.apercorr_table
        self.assertIsNotNone(nulldp.apercorr_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.apercorr_table)
        self.assertEqual( len(self.dataproduct.apercorr_table),
                          len(datacopy.apercorr_table) )
        table1 = np.asarray(self.dataproduct.apercorr_table)
        table2 = np.asarray(datacopy.apercorr_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsApertureCorrectionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                self.assertIsNotNone(readback.apercorr_table)
                self.assertEqual( len(self.dataproduct.apercorr_table),
                                  len(readback.apercorr_table) )
                original = np.asarray(self.dataproduct.apercorr_table)
                duplicate = np.asarray(readback.apercorr_table)
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
        descr = str(self.dataproduct.apercorr_table)
        self.assertIsNotNone(descr)
        del descr

class TestMiriLrsApertureCorrectionModel(unittest.TestCase):
    
    # Test the MiriLrsApertureCorrectionModel class.
        
    def setUp(self):
        # Create a typical photometric data product.
        apcorr = np.zeros([388,40])
        apcorr_err = np.zeros([388,40])
        apcorr[:,:] = np.linspace(1.5,1,num = 40) 
        width = np.arange(40) + 1
        wave = np.linspace(3.5, 14, num = 388)
        
        apcorr_table=[]
        apcorr_table.append(("FULL", wave.tolist(), 388, width.tolist(), 40, apcorr.tolist(), apcorr_err.tolist()))
        apcorr_table.append(("SLITLESSPRISM", wave.tolist(), 388, width.tolist(), 40, apcorr.tolist(), apcorr_err.tolist()))
            
        self.apcorr_table = apcorr_table

        self.dataproduct = MiriLrsApertureCorrectionModel( apcorr_table=self.apcorr_table)
        self.testfile = "MiriLrsApCorr_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.apcorr_table
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
        class_names = list(MiriLrsApertureCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('apcorr_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLrsApertureCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.apcorr_table = self.apcorr_table
        self.assertIsNotNone(nulldp.apcorr_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.apcorr_table)
        self.assertEqual( len(self.dataproduct.apcorr_table),
                          len(datacopy.apcorr_table) )
        table1 = np.asarray(self.dataproduct.apcorr_table)
        table2 = np.asarray(datacopy.apcorr_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsApertureCorrectionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                self.assertIsNotNone(readback.apcorr_table)
                self.assertEqual( len(self.dataproduct.apcorr_table),
                                  len(readback.apcorr_table) )
                original = np.asarray(self.dataproduct.apcorr_table)
                duplicate = np.asarray(readback.apcorr_table)
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
        descr = str(self.dataproduct.apcorr_table)
        self.assertIsNotNone(descr)
        del descr

class TestMiriLrsThroughputCorrectionModel(unittest.TestCase):
    
    # Test the MiriLrsThroughputCorrectionModel class.    
    def setUp(self):
        throughcorrdata_lrs = \
            [( 4.0, 0.9,  0.01),
            ( 6.0, 0.85,  0.01),
            ( 8.0, 0.8,  0.01),
            (10.0, 0.75,  0.01),
            (12.0, 0.7,  0.01),
            (14.0, 0.65,  0.01),
            (16.0, 0.6,  0.01)]
    
            
        self.throughcorr_table = throughcorrdata_lrs

        self.dataproduct = MiriLrsThroughputCorrectionModel( throughcorr_table=self.throughcorr_table)
        self.testfile = "MiriLrsThroughCorr_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.throughcorr_table
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
        class_names = list(MiriLrsThroughputCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('throughcorr_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLrsThroughputCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.throughcorr_table = self.throughcorr_table
        self.assertIsNotNone(nulldp.throughcorr_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.throughcorr_table)
        self.assertEqual( len(self.dataproduct.throughcorr_table),
                          len(datacopy.throughcorr_table) )
        table1 = np.asarray(self.dataproduct.throughcorr_table)
        table2 = np.asarray(datacopy.throughcorr_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsThroughputCorrectionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                self.assertIsNotNone(readback.throughcorr_table)
                self.assertEqual( len(self.dataproduct.throughcorr_table),
                                  len(readback.throughcorr_table) )
                original = np.asarray(self.dataproduct.throughcorr_table)
                duplicate = np.asarray(readback.throughcorr_table)
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
        descr = str(self.dataproduct.throughcorr_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriLrsPositionCorrectionModel(unittest.TestCase):
    
    # Test the MiriLrsThroughputCorrectionModel class.    
    def setUp(self):
        wave = [5.,8.,11.]
        poscorr = np.zeros([81])
        poscorr = np.append(np.linspace(0.5,1,num = 41), np.linspace(0.99,0.5, num = 40))
     
        poscorr_table=[]
        poscorr_table.append((wave[0], poscorr.tolist()))
        poscorr_table.append((wave[1], poscorr.tolist()))
        poscorr_table.append((wave[2], poscorr.tolist()))
    
            
        self.poscorr_table = poscorr_table

        self.dataproduct = MiriLrsPositionCorrectionModel( poscorr_table=self.poscorr_table)
        self.testfile = "MiriLrsPosCorr_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.poscorr_table
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
        class_names = list(MiriLrsPositionCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('poscorr_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLrsPositionCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.poscorr_table = self.poscorr_table
        self.assertIsNotNone(nulldp.poscorr_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.poscorr_table)
        self.assertEqual( len(self.dataproduct.poscorr_table),
                          len(datacopy.poscorr_table) )
        table1 = np.asarray(self.dataproduct.poscorr_table)
        table2 = np.asarray(datacopy.poscorr_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsPositionCorrectionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                self.assertIsNotNone(readback.poscorr_table)
                self.assertEqual( len(self.dataproduct.poscorr_table),
                                  len(readback.poscorr_table) )
                original = np.asarray(self.dataproduct.poscorr_table)
                duplicate = np.asarray(readback.poscorr_table)
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
        descr = str(self.dataproduct.poscorr_table)
        self.assertIsNotNone(descr)
        del descr


class TestMiriLrsPathlossCorrectionModel(unittest.TestCase):
    
    # Test the MiriLrsPathlossCorrectionModel class.    
    def setUp(self):
        wave = [5.,8.,11.]
        pathloss = np.zeros([3,40,4])
        pathloss_err = np.zeros([3,40,4])
        pathloss_table=[]
        pathloss_table.append((wave[0], pathloss[0].tolist(), pathloss_err[0].tolist()))
        pathloss_table.append((wave[1], pathloss[1].tolist(), pathloss_err[1].tolist()))
        pathloss_table.append((wave[2], pathloss[2].tolist(), pathloss_err[2].tolist()))  
            
        self.pathloss_table = pathloss_table

        self.dataproduct = MiriLrsPathlossCorrectionModel( pathloss_table=self.pathloss_table)
        self.testfile = "MiriLrsPathlossCorr_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.pathloss_table
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
        class_names = list(MiriLrsPathlossCorrectionModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('pathloss_table'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLrsPathlossCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.pathloss_table = self.pathloss_table
        self.assertIsNotNone(nulldp.pathloss_table)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2    

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.pathloss_table)
        self.assertEqual( len(self.dataproduct.pathloss_table),
                          len(datacopy.pathloss_table) )
        table1 = np.asarray(self.dataproduct.pathloss_table)
        table2 = np.asarray(datacopy.pathloss_table)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLrsPathlossCorrectionModel(self.testfile) as readback:
                self.assertEqual(self.dataproduct.meta.reftype,
                                 readback.meta.reftype)
                self.assertEqual(self.dataproduct.meta.model_type,
                                 readback.meta.model_type)
                self.assertIsNotNone(readback.pathloss_table)
                self.assertEqual( len(self.dataproduct.pathloss_table),
                                  len(readback.pathloss_table) )
                original = np.asarray(self.dataproduct.pathloss_table)
                duplicate = np.asarray(readback.pathloss_table)
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
        
        # Attempt to access the flux table pathloss attributes.
        descr = str(self.dataproduct.pathloss_table)
        self.assertIsNotNone(descr)
        del descr



# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
