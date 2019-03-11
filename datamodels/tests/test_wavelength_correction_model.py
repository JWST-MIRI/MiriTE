#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_wavelength_correction_model - Contains the unit tests for
the datamodels.miri_wavelength_correction_model module.

:History:

10 May 2017: Original version.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Oct 2018: 'N/A' used as a metadata wildcard instead of 'ANY'.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_wavelength_correction_model import MiriMrsWavelengthCorrectionModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriMrsWavelengthCorrectionModel(unittest.TestCase):
    
    # Test the MiriMrsWavelengthCorrectionModel class.
        
    def setUp(self):
        # Create a typical MiriMrsWavelengthCorrectionModel data product.
        # Sub-bands
        self.wavcorr_optical = [('1A', 0.176, 4.87, 5.82, 3320.0, 3710.0),
                                ('1B', 0.176, 5.62, 6.73, 3190.0, 3750.0)]
        self.wavcorr_xslice =  [(0.2122, 0.3718)]
        self.wavcorr_shift =   [(0.000, 0.0, 0.0),
                                (0.005, -0.0460, -0.0687),
                                (0.010, -0.0924, -0.0687)]

        self.dataproduct = MiriMrsWavelengthCorrectionModel( \
                                        wavcorr_optical=self.wavcorr_optical,
                                        wavcorr_xslice=self.wavcorr_xslice,
                                        wavcorr_shift=self.wavcorr_shift )
        self.dataproduct.set_instrument_metadata('MIRIFUSHORT', modelnam='FM',
                    detsetng='N/A', filt='N/A', channel='12', band='N/A')
        self.dataproduct.set_subarray_metadata('GENERIC')
        self.dataproduct.meta.exposure.readpatt = 'N/A'
        self.dataproduct.meta.exposure.type = 'MIR_MRS'
        self.testfile = "MiriWavelengthCorrection_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.wavcorr_optical
        del self.wavcorr_xslice
        del self.wavcorr_shift
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
        class_names = list(MiriMrsWavelengthCorrectionModel.fieldnames_optical)
        schema_names = list(self.dataproduct.get_field_names('wavcorr_optical'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_optical' class variable does not match schema")
        class_names = list(MiriMrsWavelengthCorrectionModel.fieldnames_xslice)
        schema_names = list(self.dataproduct.get_field_names('wavcorr_xslice'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_xslice' class variable does not match schema")
        class_names = list(MiriMrsWavelengthCorrectionModel.fieldnames_shift)
        schema_names = list(self.dataproduct.get_field_names('wavcorr_shift'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_shift' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsWavelengthCorrectionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.wavcorr_optical = self.wavcorr_optical
        self.assertIsNotNone(nulldp.wavcorr_optical)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        nulldp.wavcorr_xslice = self.wavcorr_xslice
        self.assertIsNotNone(nulldp.wavcorr_xslice)
        descr3 = str(nulldp)
        self.assertIsNotNone(descr3)
        nulldp.wavcorr_shift = self.wavcorr_shift
        self.assertIsNotNone(nulldp.wavcorr_shift)
        descr4 = str(nulldp)
        self.assertIsNotNone(descr4)
        del nulldp, descr1, descr2, descr3, descr4 

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.wavcorr_optical)
        self.assertEqual( len(self.dataproduct.wavcorr_optical),
                          len(datacopy.wavcorr_optical) )
        table1 = np.asarray(self.dataproduct.wavcorr_optical)
        table2 = np.asarray(datacopy.wavcorr_optical)
        assert_recarray_equal(table1, table2)
        
        self.assertIsNotNone(datacopy.wavcorr_xslice)
        self.assertEqual( len(self.dataproduct.wavcorr_xslice),
                          len(datacopy.wavcorr_xslice) )
        table1 = np.asarray(self.dataproduct.wavcorr_xslice)
        table2 = np.asarray(datacopy.wavcorr_xslice)
        assert_recarray_equal(table1, table2)
        
        self.assertIsNotNone(datacopy.wavcorr_shift)
        self.assertEqual( len(self.dataproduct.wavcorr_shift),
                          len(datacopy.wavcorr_shift) )
        table1 = np.asarray(self.dataproduct.wavcorr_shift)
        table2 = np.asarray(datacopy.wavcorr_shift)
        assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsWavelengthCorrectionModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=[],
                                       tables=['wavcorr_optical',
                                               'wavcorr_xslice',
                                               'wavcorr_shift'] )
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
        
        # Attempt to access the tables through attributes.
        descr = str(self.dataproduct.wavcorr_optical)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.wavcorr_xslice)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.wavcorr_shift)
        self.assertIsNotNone(descr)
        del descr

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
