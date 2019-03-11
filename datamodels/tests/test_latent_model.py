#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_latent_model - Contains the unit tests for the classes
in the datamodels.miri_latent_model module.

:History:

22 Aug 2013: Created.
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
12 Sep 2013: Test that the data product can be copied successfully.
24 Sep 2013: Reformatted model to include 4 latent decay tables instead of 1.
25 Sep 2013: Modified the constructor to accept a list of table objects.
             Corrected some typos and misunderstandings.
27 Sep 2013: Added some more creation tests.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
29 Aug 2014: Added test_referencefile.
25 Sep 2014: TYPE and REFTYPE are no longer identical.
11 Mar 2015: group_integration_time changed to group_time.
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

from miri.datamodels.miri_latent_model import \
    MiriLatentTable, MiriLatentDecayModel
from miri.datamodels.tests.util import assert_recarray_equal


class TestMiriLatentDecayModel(unittest.TestCase):
    
    # Test the MiriDroopModel class.
    
    def setUp(self):
        # Create a typical latent decay product.
        regime1 = "Unsaturated"
        exposure_id1 = 12029
        background1 = 8.34446
        latent_params1 = [('Fast',         355, 25000, 0.12,   8.33),
                          ('Intermediate', 365,    35, 0.01,   100),
                          ('Slow',         380,     9, 0.0029, 345)
                         ]
        self.latent_table1 = MiriLatentTable( regime1, exposure_id1,
                                              background1, latent_params1 )

        
        regime2 = "Full well"
        exposure_id2 = 12030
        background2 = 8.34446
        latent_params2 = [('Fast',        370, 28000, 0.12,   8.33),
                          ('Intermediate', 380,    16, 0.01,   100),
                          ('Slow',         380,    12, 0.0029, 345)
                          ]
        self.latent_table2 = MiriLatentTable( regime2, exposure_id2,
                                              background2, latent_params2 )
        
        latent_table_list = [self.latent_table1, self.latent_table2]
        self.dataproduct = MiriLatentDecayModel( latent_tables=latent_table_list )
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
        self.testfile = "MiriLatentDecayModel_test.fits"
        
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
        class_names = list(MiriLatentDecayModel.fieldnames)
        schema_names = list(self.dataproduct.get_field_names('latent1'))
        self.assertEqual(class_names, schema_names,
                         "'fieldnames' class variable does not match schema")

        # It must be possible to create an empty data product and fill
        # in its contents later. This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nulldp = MiriLatentDecayModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.latent1 = self.latent_table1.latent
        self.assertIsNotNone(nulldp.latent1)
        nulldp.latent2 = self.latent_table2.latent
        self.assertIsNotNone(nulldp.latent2)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2
        
        # Attempting to create a product from a corrupted table list
        # should result in an error.
        # NOTE: A bug in the JWST data model might cause an AttributeError
        # to be raised instead of a TypeError. If this happens, try a newer
        # version of the JWST data model library.
        self.assertRaises(TypeError, MiriLatentDecayModel,
                          latent_tables=42 )
        self.assertRaises(TypeError, MiriLatentDecayModel,
                          latent_tables=[1,2,3] )

        faulty_table = MiriLatentTable( 'faulty', 10.0, 25.0, 42 )
        latent_table_list = [faulty_table]
        self.assertRaises(TypeError, MiriLatentDecayModel,
                          latent_tables=latent_table_list )
        

    def test_copy(self):
        # Test that a copy can be made of the data product.
        # This will generate a warning.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        self.assertIsNotNone(datacopy.latent1)
        self.assertEqual( len(self.dataproduct.latent1),
                              len(datacopy.latent1) )
        tableA = np.asarray(self.dataproduct.latent1)
        tableB = np.asarray(datacopy.latent1)
        assert_recarray_equal(tableA, tableB)
        self.assertIsNotNone(datacopy.latent2)
        self.assertEqual( len(self.dataproduct.latent2),
                              len(datacopy.latent2) )
        tableA = np.asarray(self.dataproduct.latent2)
        tableB = np.asarray(datacopy.latent2)
        assert_recarray_equal(tableA, tableB)
        del datacopy
       
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriLatentDecayModel(self.testfile) as readback:
                # assertions TBD
                self.assertIsNotNone(readback.latent1)
                self.assertEqual( len(self.dataproduct.latent1),
                                  len(readback.latent1) )
                original = np.asarray(self.dataproduct.latent1)
                duplicate = np.asarray(readback.latent1)
                assert_recarray_equal(original, duplicate)
                self.assertIsNotNone(readback.latent2)
                self.assertEqual( len(self.dataproduct.latent2),
                                  len(readback.latent2) )
                original = np.asarray(self.dataproduct.latent2)
                duplicate = np.asarray(readback.latent2)
                assert_recarray_equal(original, duplicate)
                del readback
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error
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
