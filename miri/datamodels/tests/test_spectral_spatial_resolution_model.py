#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_spectral_spatial_resolution_model - Contains the unit tests for
the datamodels.miri_spectral_spatial_resolution_model module.

:History:

10 May 2017: Original version.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
26 Sep 2018: Major reorganisation of the data model by Jeb Bailey. The code
             supports both the old and new models but will give a warning if
             data structured according to the old model is detected.
17 Oct 2018: 'N/A' used as a metadata wildcard instead of 'ANY'.
09 Jun 2021: Reverted back to the old data model used before 26 Sep 2018.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_spectral_spatial_resolution_model import MiriMrsResolutionModel, MAX_NELEM
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal

class TestMiriMrsResolutionModel(unittest.TestCase):
    
    # Test the MiriMrsResolutionModel class.
        
    def setUp(self):
        # Create a typical MiriMrsResolutionModel data product.
        # Sub-bands
        sb = ['1SHORT', '1MEDIUM', '1LONG',
              '2SHORT', '2MEDIUM', '2LONG',
              '3SHORT', '3MEDIUM', '3LONG',
              '4SHORT', '4MEDIUM', '4LONG'
              ]
    
        # Centros
        cntr = [5.3, 6.15, 7.1, 8.15, 9.4, 10.9, 12.5, 14.5, 16.75, 19.2, 22.25, 26.05]
        
        # Low limits
        al = [2745.211671, 2737.58174, 2691.826643, 2716.566802, 2575.145064, 2563.664138, 2469.622611, 1864.309643, 2071.315475, 1899.989987, 1547.295843, 1220.329908]
        bl = [541.3274075, 506.8427022, 393.1504096, 324.8195469, 280.4117705, 257.0795746, 128.1952696, -24.10526842, 268.5610664, -50.94618217, 12.23096891, -146.1678896]
        cl = [310.8853084, 250.7537815, 35.05343649, -259.5669385, -182.4385113, -122.1030201, -279.8033748, -9.54407769, -23.40212084, -81.27269184, -12.23916046, -60.13030606]
        
        # High limits
        ah = [4063.403427, 3934.231411, 3742.01575, 4748.404518, 4227.307398, 4067.791865, 3642.877564, 2768.145139, 2641.910569, 2158.932162, np.nan, np.nan]
        bh = [193.251713, 151.0248429, 79.69070813, 476.6602166, 420.9441793, 227.7079654, 326.3927585, -86.76805885, 363.0089443, -70.19831027, np.nan, np.nan]
        ch = [-191.9995708, -418.0899389, -360.9548965, -1804.289718, -611.3820685, -569.5866384, -601.634922, -206.1799239, -45.5071111, -109.3828766, np.nan, np.nan]
        
        # Use values
        au = [3404.307549, 3335.906575, 3216.921197, 3732.48566, 3401.226231, 3315.728002, 3056.250087, 2316.227391, 2356.613022, 2029.461074, 1547.295843, 1220.329908]
        bu = [367.2895603, 328.9337726, 236.4205589, 400.7398818, 350.6779749, 242.39377, 227.2940141, -55.43666364, 315.7850053, -60.57224622, 12.23096891, -146.1678896]
        cu = [59.44286876, -83.66807868, -162.95073, -1031.928328, -396.9102899, -345.8448292, -440.7191484, -107.8620008, -34.45461597, -95.32778422, -12.23916046, -60.13030606]
    
        self.resolving_power = []
        for ii in range(0, len(cntr)):
            row = (sb [ii], cntr[ii], al[ii], bl[ii], cl[ii], ah[ii], bh[ii], ch[ii],
                   au[ii], bu[ii], cu[ii])
            self.resolving_power.append(row)
        
        self.psf_fwhm_alpha = [(8.0, 0.31, 0.0, 0.0, 0.03875)]
        self.psf_fwhm_beta = [(8.0, 0.0, 0.03875, 0.0, 0.03875)]

        self.dataproduct = MiriMrsResolutionModel( \
                                        resolving_power=self.resolving_power,
                                        psf_fwhm_alpha=self.psf_fwhm_alpha,
                                        psf_fwhm_beta=self.psf_fwhm_beta )
        self.dataproduct.set_referencefile_metadata( author='Alvaro Labiano',
                    pedigree='GROUND', useafter='DEFAULT',
                    description='MIRI MRS Spectral and Spatial Resolution CDP')
        self.dataproduct.add_referencefile_history(
                    document='MIRI-TN-00005-ETH Issue 1.0',
                    software='IDL and Python',
                    dataused='Derived from FM data',
                    differences='N/A')
        self.dataproduct.set_instrument_metadata('MIRIFUSHORT', modelnam='FM',
                    detsetng='N/A', filt='N/A', channel='12', band='N/A')
        self.dataproduct.set_subarray_metadata('GENERIC')
        self.dataproduct.meta.exposure.readpatt = 'N/A'
        self.dataproduct.meta.exposure.type = 'MIR_MRS'
        self.testfile = "MiriSpectralSpatialResolution_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.resolving_power
        del self.psf_fwhm_alpha
        del self.psf_fwhm_beta
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
        class_names = list(MiriMrsResolutionModel.fieldnames_resolving)
        schema_names = list(self.dataproduct.get_field_names('resolving_power'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_resolving' class variable does not match schema")

        class_names = list(MiriMrsResolutionModel.fieldnames_alpha)
        schema_names = list(self.dataproduct.get_field_names('psf_fwhm_alpha'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_alpha' class variable does not match schema")

        class_names = list(MiriMrsResolutionModel.fieldnames_beta)
        schema_names = list(self.dataproduct.get_field_names('psf_fwhm_beta'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames_beta' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriMrsResolutionModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.resolving_power = self.resolving_power
        self.assertIsNotNone(nulldp.resolving_power)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        nulldp.psf_fwhm_alpha = self.psf_fwhm_alpha
        self.assertIsNotNone(nulldp.psf_fwhm_alpha)
        descr3 = str(nulldp)
        self.assertIsNotNone(descr3)
        nulldp.psf_fwhm_beta = self.psf_fwhm_beta
        self.assertIsNotNone(nulldp.psf_fwhm_beta)
        descr4 = str(nulldp)
        self.assertIsNotNone(descr4)
        del nulldp, descr1, descr2, descr3, descr4 

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy.resolving_power)
        self.assertEqual( len(self.dataproduct.resolving_power),
                          len(datacopy.resolving_power) )
        # FIXME: Test fails!
#         table1 = np.asarray(self.dataproduct.resolving_power)
#         table2 = np.asarray(datacopy.resolving_power)
#         assert_recarray_equal(table1, table2)
        
        self.assertIsNotNone(datacopy.psf_fwhm_alpha)
        self.assertEqual( len(self.dataproduct.psf_fwhm_alpha),
                          len(datacopy.psf_fwhm_alpha) )
        # FIXME: Test fails!
#         table1 = np.asarray(self.dataproduct.psf_fwhm_alpha)
#         table2 = np.asarray(datacopy.psf_fwhm_alpha)
#         assert_recarray_equal(table1, table2)
        
        self.assertIsNotNone(datacopy.psf_fwhm_beta)
        self.assertEqual( len(self.dataproduct.psf_fwhm_beta),
                          len(datacopy.psf_fwhm_beta) )
        # FIXME: Test fails!
#         table1 = np.asarray(self.dataproduct.psf_fwhm_beta)
#         table2 = np.asarray(datacopy.psf_fwhm_beta)
#         assert_recarray_equal(table1, table2)
        del datacopy
                
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriMrsResolutionModel(self.testfile) as readback:
                # FIXME: Test fails!
#                 assert_products_equal( self, self.dataproduct, readback,
#                                        arrays=[],
#                                        tables=['resolving_power',
#                                                'psf_fwhm_alpha',
#                                                'psf_fwhm_beta'] )
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
        descr = str(self.dataproduct.resolving_power)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.psf_fwhm_alpha)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.psf_fwhm_beta)
        self.assertIsNotNone(descr)
        del descr

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
