#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_measured_model - Contains the unit tests for the classes
in the datamodels.miri_measured_model module.

:History:

15 Jan 2013: Created.
21 Jan 2013: Warning messages controlled with Python warnings module.
05 Feb 2013: File closing problem solved by using "with" context manager.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
23 Apr 2013: Modified to keep up with behaviour of jwst_lib model.
             Uninitialised arrays now have the same size and shape as the
             data array but are full of default values.
26 Apr 2013: File closing problem has returned!
13 May 2013: Added MiriSlopeModel to describe MIRI slope data
             (which is different from "ImageModel" data because it
             preserves integrations). N.B. FINAL MODEL IS TBD.
04 Jun 2013: Shortened the names of the ramp, slope and image models.
10 Jun 2013: Added more metadata tests.
02 Jul 2013: MiriCubeModel added.
29 Jul 2013: stats() method added.
14 Aug 2013: Updated ramp model test to include groupdq and pixeldq
02 Sep 2013: Compare numpy record arrays in a way that it independent
             of the byte ordering.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
12 Sep 2013: Test that the data product can be copied successfully.
04 Oct 2013: Changed default field_def table to use MIRI reserved flags.
07 Oct 2013: GROUP_DEF table added to MIRI ramp data. Test MiriRampModel
             for masking and arithmetic operations.
24 Feb 2014: Instrument name (INSTRUME) changed from meta.instrument.type to
             meta.instrument.name.
27 Feb 2014: Added extra data arrays to MiriSlopeModel test.
04 Mar 2014: Added set_housekeeping_metadata.
25 Jun 2014: field_def and group_def changed to dq_def and groupdq_def.
             field_def for ramp data changed to pixeldq_def.
21 Jul 2014: IM, and LW detectors changed to MIRIMAGE and MIRIFULONG.
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
07 Nov 2014: The data model now raises an IOError when an invalid file
             path is provided.
11 Mar 2015: group_integration_time changed to group_time.
11 Jun 2015: Added a history record test.
09 Jul 2015: Reference output array (refout) added to MiriRampModel schema.
19 Aug 2015: Removed MiriImageModel and MiriCubeModel.
07 Oct 2015: Made exception catching Python 3 compatible.
08 Apr 2016: Removed obsolete FIXME statements.
04 May 2016: ERR array removed from ramp data model.
31 Aug 2016: Change exception detected when creating a data model with an
             invalid initialiser.
15 Jun 2017: Observation and target metadata is appropriate for ramp and
             slope data only.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
13 Sep 2017: Updated "not a file name" test to match the new behaviour of
             JWST pipeline version 0.7.8rc2
27 Apr 2018: Corrected bug in get_history() length test.
27 Jun 2018: Removed unused arrays.
15 Feb 2018: Check that the DQ_DEF table has the correct fieldnames.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import unittest
import warnings

import numpy as np

# Import the JWST master data quality flag definitions
from miri.datamodels.dqflags import master_flags, pixeldq_flags, \
    groupdq_flags

from miri.datamodels.miri_measured_model import MiriMeasuredModel, \
    MiriRampModel, MiriSlopeModel
from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal


class TestMiriMeasuredModel(unittest.TestCase):
    
    # Test the MiriMeasuredModel class

    def setUp(self):
        # Create a 64x64 simple MiriMeasuredModel object, with no error
        # or quality arrays.
        self.data = np.linspace(0.0, 100000.0, 64*64)
        self.data.shape = [64,64]
        self.simpleproduct = MiriMeasuredModel(data=self.data)
        # Add some example metadata.
        self.simpleproduct.set_housekeeping_metadata('MIRI EC', 'Joe Bloggs',
                                                     'V1.0')
        self.simpleproduct.set_instrument_metadata(detector='MIRIMAGE',
                                                   filt='F560W',
                                                   ccc_pos='OPEN',
                                                   deck_temperature=10.0,
                                                   detector_temperature=7.0)
        self.simpleproduct.set_exposure_metadata(readpatt='SLOW',
                                                 nints=1, ngroups=10,
                                                 frame_time=30.0,
                                                 integration_time=30.0,
                                                 group_time=300.0,
                                                 reset_time=0, frame_resets=3)
        
        # Create a more complex MiriMeasuredModel object from primary,
        # error and quality arrays.
        self.primary = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        self.error =   [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        self.quality = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.dataproduct = MiriMeasuredModel(data=self.primary,
                                             err=self.error,
                                             dq=self.quality,
                                             dq_def=master_flags)
        # Add some example metadata.
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)

        self.testfile1 = "MiriMeasuredModel_test1.fits"
        self.testfile2 = "MiriMeasuredModel_test2.fits"
        self.tempfiles = [self.testfile1, self.testfile2]
      
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        del self.primary, self.error, self.quality
        del self.simpleproduct
        del self.data
        # Remove temporary files, if they exist and if able to.
        for tempfile in self.tempfiles: 
            if os.path.isfile(tempfile):
                try:
                    os.remove(tempfile)
                except Exception as e:
                    strg = "Could not remove temporary file, " + tempfile + \
                        "\n   " + str(e)
                    warnings.warn(strg)
        del self.tempfiles
        
    def test_creation(self):
        # Check that the DQ_DEF field names in the class variable are the same
        # as the ones declared in the schema.
        dq_def_names = list(MiriMeasuredModel.dq_def_names)
        schema_names = list(self.dataproduct.get_field_names('dq_def'))
        self.assertEqual(dq_def_names, schema_names,
                         "'dq_def_names' class variable does not match schema")

        # Test that the error and quality arrays are optional.
        a2 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        b2 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c2 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
                
        # 1) Data array only. Data array must exist and be non-empty.
        # Other arrays should exist and be the same size and shape as the
        # data array. They should be full of default values.
        newdp1 = MiriMeasuredModel(data=a2)
        self.assertIsNotNone(newdp1.data)
        self.assertGreater(len(newdp1.data), 0)
        self.assertIsNotNone(newdp1.err)
        self.assertEqual(newdp1.err.shape, newdp1.data.shape)
        # Assumes default is 0.0 - see schema
        self.assertAlmostEqual(np.mean(newdp1.err), 0.0)
        self.assertIsNotNone(newdp1.dq)
        self.assertEqual(newdp1.dq.shape, newdp1.dq.shape)
        # Assumes default is 0 - see schema
        self.assertEqual(np.mean(newdp1.dq), 0)
        descr1 = str(newdp1)
        self.assertIsNotNone(descr1)
        del newdp1, descr1
        
        # 2) Data and error arrays only. Data and error arrays must exist
        # and be non-empty. Quality array should exist but be the same
        # size and shape as the data array. It should be full of default
        # values.
        newdp2 = MiriMeasuredModel(data=a2, err=b2)
        self.assertIsNotNone(newdp2.data)
        self.assertGreater(len(newdp2.data), 0)
        self.assertIsNotNone(newdp2.err)
        self.assertEqual(newdp2.err.shape, newdp2.data.shape)
        # The error array must not be full of default values.
        self.assertNotAlmostEqual(np.mean(newdp2.err), 0.0)
        self.assertIsNotNone(newdp2.dq)
        self.assertEqual(newdp2.dq.shape, newdp2.dq.shape)
        # Assumes default is 0 - see schema
        self.assertEqual(np.mean(newdp2.dq), 0)
        descr2 = str(newdp2)
        self.assertIsNotNone(descr2)
        del newdp2, descr2

        # 3) Data, error and quality arrays. All arrays must exist,
        # be non-empty and be the same size and shape.
        newdp3 = MiriMeasuredModel(data=a2, err=b2, dq=c2)
        self.assertIsNotNone(newdp3.data)
        self.assertGreater(len(newdp3.data), 0)
        self.assertIsNotNone(newdp3.err)
        self.assertEqual(newdp3.err.shape, newdp3.data.shape)
        # The error array must not be full of default values.
        self.assertNotAlmostEqual(np.mean(newdp3.err), 0.0)
        self.assertIsNotNone(newdp3.dq)
        self.assertEqual(newdp3.dq.shape, newdp3.dq.shape)
        # The quality array must not be full of default values.
        self.assertNotEqual(np.mean(newdp3.dq), 0)
        descr3 = str(newdp3)
        self.assertIsNotNone(descr3)
        del newdp3, descr3
        
        # It should be possible to set up an empty data product with
        # a specified shape. All three arrays should be initialised to
        # the same shape.
        emptydp = MiriMeasuredModel( (4,4) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (4,4))
        self.assertIsNotNone(emptydp.err)
        self.assertEqual(emptydp.err.shape, (4,4))
        self.assertIsNotNone(emptydp.dq)
        self.assertEqual(emptydp.dq.shape, (4,4))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr
        
        # A null data product can also be created and populated
        # with data later.
        nulldp = MiriMeasuredModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.data = np.asarray(a2)
        self.assertIsNotNone(nulldp.err)
        self.assertIsNotNone(nulldp.dq)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2
        
        # A scalar data product is possible, even if of little use.
        scalardp = MiriMeasuredModel( data=42 )
        self.assertEqual(scalardp.data, 42)
        self.assertIsNotNone(scalardp.err)
        self.assertIsNotNone(scalardp.dq)
        descr = str(scalardp)
        self.assertIsNotNone(descr)
        del scalardp, descr
        
        # Attempts to create a data product from invalid data types
        # and stupid values must be detected.
        # NOTE: A bug in the JWST data model might cause an AttributeError
        # to be raised instead of a ValueError. If this happens, try a newer
        # version of the JWST data model library.
        self.assertRaises(ValueError, MiriMeasuredModel, init=[])
        self.assertRaises(ValueError, MiriMeasuredModel, init=42)
        self.assertRaises(ValueError, MiriMeasuredModel, init='not a file name')
        self.assertRaises(IOError, MiriMeasuredModel, init='nosuchfile.fits')
        #self.assertRaises(ValueError, MiriMeasuredModel, init='')
        self.assertRaises(ValueError, MiriMeasuredModel, data='badstring')

    def test_metadata(self):
        # Check the dataproducts contain metadata
        # First test the basic STScI FITS keyword lookup method.
        kwstrg = self.simpleproduct.find_fits_keyword('TELESCOP',
                                                      return_result=True)
        self.assertIsNotNone(kwstrg)
        # kwstrg is a list - assume the first entry is what we want.
        telname = self.simpleproduct[kwstrg[0]]
        self.assertEqual(telname, 'JWST')
        # Accessing the tree structure directly should also work.
        telname = self.simpleproduct.meta.telescope
        self.assertEqual(telname, 'JWST')
        # An alternative lookup provided by the MIRI data model.
        telname = self.simpleproduct.get_fits_keyword('TELESCOP')
        self.assertEqual(telname, 'JWST')
        
        kwstrg = self.simpleproduct.find_fits_keyword('INSTRUME',
                                                      return_result=True)
        self.assertIsNotNone(kwstrg)
        insname = self.simpleproduct[kwstrg[0]]
        self.assertEqual(insname, 'MIRI')
        insname = self.simpleproduct.meta.instrument.name
        self.assertEqual(insname, 'MIRI')
        insname = self.simpleproduct.get_fits_keyword('INSTRUME')
        self.assertEqual(insname, 'MIRI')
        
        # Add some history records and check they exist.
        self.simpleproduct.add_history('History 1')
        self.simpleproduct.add_history('History 2')
        self.simpleproduct.add_history('History 3')
        self.assertGreaterEqual(len(self.simpleproduct.get_history()), 3)
        strg = self.simpleproduct.get_history_str()
        self.assertIsNotNone(strg)
        self.assertGreater(len(strg), 0)
    
    def test_content(self):
        # The data, err and dq attributes are aliases for the primary,
        # error and quality arrays
        self.assertTrue( np.allclose(self.primary, self.dataproduct.data) )
        self.assertTrue( np.allclose(self.error, self.dataproduct.err) )
        self.assertTrue( np.allclose(self.quality, self.dataproduct.dq) )

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
        
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.simpleproduct.save(self.testfile1, overwrite=True)
            with MiriMeasuredModel(self.testfile1) as readback:
                self.assertTrue( np.allclose(self.simpleproduct.data,
                                             readback.data) )
                del readback

            self.dataproduct.save(self.testfile2, overwrite=True)
            with MiriMeasuredModel(self.testfile2) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq'],
                                       tables='dq_def' )
                del readback

    def test_asciiio(self):
        # Check that the data products can be written to an ASCII
        # file and read back again without changing the data.
        # TODO: At the moment jwst_lib only supports FITS I/O
        pass
#         # Suppress metadata warnings
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             self.simpleproduct.save(self.testfile_ascii, overwrite=True)
#             with MiriMeasuredModel(self.testfile_ascii) as readback:
#                 self.assertTrue( np.allclose(self.simpleproduct.data,
#                                              readback.data) )
#                 del readback

    def test_masking(self):
        # The DQ array must mask off bad values in the SCI and ERR arrays.
        a2 = [[10,999,10,999], [999,10,10,999], [10,10,999,10]]
        b2 = [[1,99,1,99],     [99,1,1,99],     [1,1,99,1]]
        c2 = [[0,1,0,1],       [1,0,0,1],       [0,0,1,0]]

        # Without a DQ array (assuming the default quality value is 0)
        # the SCI and ERR arrays are not masked, so their averages
        # include the 999s and are greater than they ought to be.
        newdp = MiriMeasuredModel(data=a2, err=b2)
        meandata = np.mean(newdp.data_masked)
        self.assertGreater(meandata, 10)
        meanerr = np.mean(newdp.err_masked)
        self.assertGreater(meanerr, 1)
        
        # The addition of the quality data should cause the SCI and ERR
        # arrays to be masked off and give the correct average.
        newdp2 = MiriMeasuredModel(data=a2, err=b2, dq=c2)
        meandata2 = np.mean(newdp2.data_masked)
        self.assertAlmostEqual(meandata2, 10)
        meanerr2 = np.mean(newdp2.err_masked)
        self.assertAlmostEqual(meanerr2, 1)
        
        del newdp, newdp2

    def test_arithmetic(self):
        a2 = [[90,80,70,60],[50,40,30,20],[10,0,-10,-20]]
        b2 = [[1,2,3,4],[5,6,7,8],[9,10,11,12]]
        c2 = [[0,1,1,0],[0,2,0,2],[1,0,1,0]]
        newdp = MiriMeasuredModel(data=a2, err=b2, dq=c2)
        
        # Self-subtraction of the simple product. The result
        # should be zero.
        newsimple = self.simpleproduct - self.simpleproduct
        self.assertAlmostEqual(newsimple.data.all(), 0.0)
        del newsimple
        
        # Scalar addition
        result = self.dataproduct + 42
        test1 = self.dataproduct.data + 42
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product addition
        result = self.dataproduct + newdp
        test1 = self.dataproduct.data + newdp.data
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        # Test that error arrays are combined properly - at least for
        # a couple of unmasked points.
        expectedsq = self.error[1][0]*self.error[1][0] + b2[1][0]*b2[1][0]
        actualsq = result.err[1,0]*result.err[1,0]
        self.assertAlmostEqual(expectedsq, actualsq)
        expectedsq = self.error[2][1]*self.error[2][1] + b2[2][1]*b2[2][1]
        actualsq = result.err[2,1]*result.err[2,1]
        self.assertAlmostEqual(expectedsq, actualsq)
        del result
        
        # Scalar subtraction
        result = self.dataproduct - 42
        test1 = self.dataproduct.data - 42
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product subtraction
        result = self.dataproduct - newdp
        test1 = self.dataproduct.data - newdp.data
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        # Test that error arrays are combined properly - at least for
        # a couple of unmasked points.
        expectedsq = self.error[1][0]*self.error[1][0] + b2[1][0]*b2[1][0]
        actualsq = result.err[1,0]*result.err[1,0]
        self.assertAlmostEqual(expectedsq, actualsq)
        expectedsq = self.error[2][1]*self.error[2][1] + b2[2][1]*b2[2][1]
        actualsq = result.err[2,1]*result.err[2,1]
        self.assertAlmostEqual(expectedsq, actualsq)
        del result
        
        # Addition and subtraction should cancel each other out
        result = self.dataproduct + newdp - newdp
        test1 = self.dataproduct.data
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        del result
        
        # Scalar multiplication
        result = self.dataproduct * 3
        test1 = self.dataproduct.data * 3
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        del result

        # Data product multiplication
        result = self.dataproduct * newdp
        test1 = self.dataproduct.data * newdp.data
        test2 = result.data
        self.assertEqual(test1.all(), test2.all())
        err1 = self.dataproduct.err
        da1 = self.dataproduct.data
        err2 = newdp.err
        da2 = newdp.data
        expectedErr = np.sqrt(err1 * err1 * da2 * da2 + err2 * err2 * da1 * da1)
        self.assertTrue(np.array_equal(expectedErr, result.err))
        
        del result, da1, da2, err1, err2, expectedErr

        # Scalar division
        result = self.dataproduct / 3.0
        test1 = self.dataproduct.data / 3.0
        test2 = result.data
        self.assertAlmostEqual(test1.all(), test2.all())
        del test1, test2, result
        
        # Division by zero
        self.assertRaises(ValueError, self.dataproduct.__truediv__, 0.0)

        # Data product division
        #print("NOTE: The following test is expected to generate run time warnings.")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = self.dataproduct / newdp
            test1 = self.dataproduct.data / newdp.data
            test2 = result.data
            self.assertEqual(test1.all(), test2.all())
            # Test Juergen Schreiber error propagation
            dat = self.dataproduct.data[1][1]
            newdat = newdp.data[1][1]
            resultErr = result.err[1][1]
            dpErr = self.dataproduct.err[1][1]
            newdpErr = newdp.err[1][1]
            expectErr = np.sqrt( dpErr * dpErr/(newdat * newdat) + \
                                 newdpErr * newdpErr * dat * dat / \
                                 (newdat * newdat * newdat * newdat))            
 
            self.assertEqual(expectErr, resultErr)
            del test1, test2, result
        
            # More complex arithmetic should be possible.
            newdp2 = newdp * 2
            newdp3 = newdp * 3
            newdp4 = newdp2 + newdp3
            result = ((self.dataproduct - newdp) * newdp2 / newdp3) + newdp4
            del newdp, newdp2, newdp3, newdp4
            del result

    def test_broadcasting(self):
        # Test that operations where the broadcasting of one array
        # onto a similar shaped array work.
        a4x3 = [[90,80,70,60],[50,40,30,20],[10,0,-10,-20]]
        b4x3 = [[1,2,3,4],[5,6,7,8],[9,10,11,12]]
        #c4x3 = [[0,1,0,0],[0,0,1,0],[1,0,0,1]]
        
        a4x1 = [4,3,2,1]
        b4x1 = [1,2,1,2]
        c4x1 = [0,1,0,0]
        
        #a5x1 = [5,4,3,2,1]
        #b5x1 = [1,2,3,2,1]
        c5x1 = [0,1,0,0,1]
        
        # Create an object with 4x3 primary and error arrays but a 4x1
        # quality array. This should succeed because the quality array
        # is broadcastable.
        newdp1 = MiriMeasuredModel(data=a4x3, err=b4x3, dq=c4x1)
        self.assertTrue( np.allclose(a4x3, newdp1.data) )
        self.assertTrue( np.allclose(b4x3, newdp1.err) )
        self.assertTrue( np.allclose(c4x1, newdp1.dq) )

        # 5x1 is not broadcastable onto 4x3 and this statement should fail.
        # NOTE: Unfortunately this test also issues a warning message,
        # "'MiriMeasuredModel' object has no attribute '_real_cls'".
        # Turning off warnings does not stop this message from appearing.
        self.assertRaises(TypeError, MiriMeasuredModel, data=a4x3,
                          error=b4x3, quality=c5x1)
        
        # Combine two broadcastable object mathematically.
        # The + and - operations should be commutative and the result
        # should be saveable to a FITS file.
        newdp2 = MiriMeasuredModel(data=a4x1, err=b4x1, dq=c4x1)
        
        result1 = newdp1 + newdp2
        result2 = newdp2 + newdp1
        self.assertEqual(result1.data.shape, result2.data.shape)
        self.assertTrue( np.allclose(result1.data, result2.data) )
        self.assertTrue( np.allclose(result1.err, result2.err) )
        self.assertTrue( np.allclose(result1.dq, result2.dq) )
        result1.save(self.testfile1, overwrite=True)
        result2.save(self.testfile2, overwrite=True)
        del result1, result2

        result1 = newdp1 * newdp2
        result2 = newdp2 * newdp1
        self.assertEqual(result1.data.shape, result2.data.shape)
        self.assertTrue( np.allclose(result1.data, result2.data) )
        self.assertTrue( np.allclose(result1.err, result2.err) )
        self.assertTrue( np.allclose(result1.dq, result2.dq) )
        result1.save(self.testfile1, overwrite=True)
        result2.save(self.testfile2, overwrite=True)
        del result1, result2

        # The - and / operations are not commutative, but the data shape
        # should be consistent and the quality arrays should be combined
        # in the same way.
        result1 = newdp1 - newdp2
        result2 = newdp2 - newdp1
        self.assertEqual(result1.data.shape, result2.data.shape)
        self.assertTrue( np.allclose(result1.err, result2.err) )
        self.assertTrue( np.allclose(result1.dq, result2.dq) )
        result1.save(self.testfile1, overwrite=True)
        result2.save(self.testfile2, overwrite=True)
        del result1, result2

        result1 = newdp1 / newdp2
        result2 = newdp2 / newdp1
        self.assertEqual(result1.data.shape, result2.data.shape)
        # The errors resulting from division depend on the order
        # of the operation.
        self.assertTrue( np.allclose(result1.dq, result2.dq) )
        result1.save(self.testfile1, overwrite=True)
        result2.save(self.testfile2, overwrite=True)
        del result1, result2

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr = str(self.simpleproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.simpleproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = self.simpleproduct.stats()
        self.assertIsNotNone(descr)
        del descr

        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = self.dataproduct.stats()
        self.assertIsNotNone(descr)
        del descr
        
        # Attempt to access the SCI, ERROR and DQ arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.err)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


class TestMiriRampModel(unittest.TestCase):
    
    # Most of the necessary tests are already carried out by
    # the TestMiriMeasuredModel class.

    def setUp(self):
        # Create a ramp data product.
        # NOTE: A ramp product does not contain an ERR array.
        self.a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        self.c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.c2 = [[0,1,1,0],     [1,0,0,1],     [1,0,1,0]]
        self.acube = [self.a1,self.a1,self.a1]
        self.ccube = [self.c1,self.c2,self.c1]
        self.ahyper = [self.acube,self.acube]
        self.chyper = [self.ccube,self.ccube]
        self.refout = np.ones_like(self.chyper)
        self.dataproduct = MiriRampModel(data=self.ahyper, refout=self.refout,
                                         pixeldq=self.c1,
                                         dq_def=pixeldq_flags,
                                         groupdq=self.chyper,
                                         groupdq_def=groupdq_flags)
        # Add some example metadata.
        self.dataproduct.set_housekeeping_metadata('MIRI EC', 'Joe Bloggs',
                                                   'V1.0')
        self.dataproduct.set_observation_metadata()
        self.dataproduct.set_target_metadata(0.0, 0.0)
        self.dataproduct.set_instrument_metadata(detector='MIRIFULONG',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriRampModel_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.a1, self.c1, self.c2
        del self.acube, self.ccube
        del self.ahyper, self.chyper
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
        # Test that any of the quality arrays are optional.
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        bcube = [b1,b1,b1]
        bhyper = [bcube,bcube]
                
        # 1) Data array only. Data array must exist and be non-empty.
        # The quality arrays must be 2-D and 4-D.
        # Unspecified arrays must be filled with default values. 
        newdp1 = MiriRampModel(data=self.ahyper)
        self.assertIsNotNone(newdp1.data)
        self.assertGreater(len(newdp1.data), 0)
        # Assumes default is 0.0 - see schema
        self.assertIsNotNone(newdp1.pixeldq)
        self.assertTrue(newdp1.pixeldq.ndim == 2)
        # Assumes default is 0 - see schema
        # FIXME: The pixeldq array ends up containing null values.
        #self.assertEqual(np.mean(newdp1.pixeldq), 0)
        self.assertIsNotNone(newdp1.groupdq)
        self.assertTrue(newdp1.groupdq.ndim == 4)
        # Assumes default is 0 - see schema
        self.assertEqual(np.mean(newdp1.groupdq), 0)
        descr1 = str(newdp1)
        del newdp1, descr1

        # 2) Data and both quality arrays. All arrays must exist,
        # be non-empty and be the shape specified.
        newdp3 = MiriRampModel(data=self.ahyper, pixeldq=self.c1,
                               groupdq=self.chyper)
        self.assertIsNotNone(newdp3.data)
        self.assertGreater(len(newdp3.data), 0)
        # The pixeldq array must not be full of default values.
        self.assertIsNotNone(newdp3.pixeldq)
        self.assertTrue(newdp3.pixeldq.ndim == 2)
        self.assertNotEqual(np.mean(newdp3.pixeldq), 0)
        self.assertIsNotNone(newdp3.groupdq)
        self.assertTrue(newdp3.groupdq.ndim == 4)
        # The groupdq array must not be full of default values.
        self.assertNotEqual(np.mean(newdp3.groupdq), 0)
        descr3 = str(newdp3)
        del newdp3, descr3

        # 3) Data and pixeldq array only. All arrays must exist,
        # be non-empty and be the shape specified.
        newdp4 = MiriRampModel(data=self.ahyper, pixeldq=self.c1)
        self.assertIsNotNone(newdp4.data)
        self.assertGreater(len(newdp4.data), 0)
        # The pixeldq array must not be full of default values.
        self.assertIsNotNone(newdp4.pixeldq)
        self.assertTrue(newdp4.pixeldq.ndim == 2)
        self.assertNotEqual(np.mean(newdp4.pixeldq), 0)
        self.assertIsNotNone(newdp4.groupdq)
        self.assertTrue(newdp4.groupdq.ndim == 4)
        descr4 = str(newdp4)
        del newdp4, descr4

        # 4) Data and groupdq array only. All arrays must exist,
        # be non-empty and be the shape specified.
        newdp5 = MiriRampModel(data=self.ahyper, groupdq=self.chyper)
        self.assertIsNotNone(newdp5.data)
        self.assertGreater(len(newdp5.data), 0)
        self.assertIsNotNone(newdp5.pixeldq)
        self.assertTrue(newdp5.pixeldq.ndim == 2)
        # The groupdq array must not be full of default values.
        self.assertIsNotNone(newdp5.groupdq)
        self.assertTrue(newdp5.groupdq.ndim == 4)
        # The groupdq array must not be full of default values.
        self.assertNotEqual(np.mean(newdp5.groupdq), 0)
        descr5 = str(newdp5)
        del newdp5, descr5

        # It should be possible to set up an empty data product with
        # a specified 4-D shape. Data array should be
        # initialised to the same shape.
        emptydp = MiriRampModel( (2,2,2,2) )
        self.assertIsNotNone(emptydp.data)
        self.assertEqual(emptydp.data.shape, (2,2,2,2))
        self.assertIsNotNone(emptydp.pixeldq)
        #self.assertEqual(emptydp.pixeldq.shape, (2,2))
        self.assertIsNotNone(emptydp.groupdq)
        self.assertEqual(emptydp.groupdq.shape, (2,2,2,2))
        descr = str(emptydp)
        self.assertIsNotNone(descr)
        del emptydp, descr
        
        # A null data product can also be created and populated
        # with data later.
        nulldp = MiriRampModel( )
        descr1 = str(nulldp)
        self.assertIsNotNone(descr1)
        nulldp.data = np.asarray(self.ahyper)
        self.assertIsNotNone(nulldp.pixeldq)
        self.assertIsNotNone(nulldp.groupdq)
        descr2 = str(nulldp)
        self.assertIsNotNone(descr2)
        del nulldp, descr1, descr2
        
        # Creating an object with other than 4 dimensions must fail.
        a1d = [10,20,30,40]
        c1d = [1,0,0,0]
        self.assertRaises(ValueError, MiriRampModel, data=a1d, pixeldq=c1d)

        a2d = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        c2d = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        self.assertRaises(ValueError, MiriRampModel, data=a2d, groupdq=c2d)

        a3d = [a2d, a2d, a2d]
        c3d = [c2d, c2d, c2d]
        self.assertRaises(ValueError, MiriRampModel, data=a3d, pixeldq=c3d)
        self.assertRaises(ValueError, MiriRampModel, data=a3d, groupdq=c3d)

        # The pixeldq array must be 2-D.
        self.assertRaises(ValueError, MiriRampModel, data=self.ahyper,
                          pixeldq=self.ccube)
        # The groupdq array must be 4-D.
        self.assertRaises(ValueError, MiriRampModel, data=self.ahyper,
                          groupdq=self.c1)

    def test_masking(self):
        # Ramp data must have a dq array which gives a view of one
        # or both of the pixeldq and groupdq masks
        self.assertIsNotNone(self.dataproduct.dq)

        # Create a data product masked by the pixeldq array.
        # The dq and pixeldq arrays must be the same
        mask1 = MiriRampModel(data=self.ahyper, pixeldq=self.c1,
                              groupdq=self.chyper, maskwith='pixeldq')
        self.assertIsNotNone(mask1.pixeldq)
        self.assertGreater(len(mask1.pixeldq), 0)
        self.assertIsNotNone(mask1.dq)
        self.assertGreater(len(mask1.dq), 0)
        self.assertEqual(mask1.dq.shape, mask1.pixeldq.shape)
        self.assertTrue(np.all( mask1.dq == mask1.pixeldq ))
        del mask1

        # Create a data product masked by the groupdq array.
        # The dq and groupdq arrays must be the same
        mask2 = MiriRampModel(data=self.ahyper, pixeldq=self.c1,
                              groupdq=self.chyper, maskwith='groupdq')
        self.assertIsNotNone(mask2.groupdq)
        self.assertGreater(len(mask2.groupdq), 0)
        self.assertIsNotNone(mask2.dq)
        self.assertGreater(len(mask2.dq), 0)
        self.assertEqual(mask2.dq.shape, mask2.groupdq.shape)
        self.assertTrue(np.all( mask2.dq == mask2.groupdq ))
        del mask2

        # Create a data product masked by both pixeldq and groupdq arrays.
        # The result must have the same shape as the groupdq array but be
        # a combination of both masks.
        mask3 = MiriRampModel(data=self.ahyper, pixeldq=self.c1,
                              groupdq=self.chyper, maskwith='both')
        self.assertIsNotNone(mask3.pixeldq)
        self.assertGreater(len(mask3.pixeldq), 0)
        self.assertIsNotNone(mask3.groupdq)
        self.assertGreater(len(mask3.groupdq), 0)
        self.assertIsNotNone(mask3.dq)
        self.assertGreater(len(mask3.dq), 0)
        self.assertEqual(mask3.dq.shape, mask3.groupdq.shape)
        expected = mask3.groupdq | mask3.pixeldq
        self.assertTrue(np.all( mask3.dq == expected ))
        del mask3

    def test_arithmetic(self):
        # The ramp data model supports all the arithmetic operations
        # supported by the MiriMeasuredModel. The following are exceptions
        # specific to the ramp model.
        
        # Create a data model in which the DATA and DQ arrays have different
        # shapes.
        testdp = MiriRampModel(data=self.ahyper, pixeldq=self.c1,
                               groupdq=self.chyper, maskwith='both')
        descr = str(testdp)
        self.assertIsNotNone(descr)
        del descr
        
        # Suppress warning about the DQ array being propagated only from GROUPDQ
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check the product can be combined with itself
            double = testdp * 2.0
            self.assertIsNotNone(double.data)
            self.assertGreater(len(double.data), 0)
            expected = double.data * 2.0
            self.assertTrue(np.all( (double.data - expected) < 0.001 ))
            descr = str(double)
            self.assertIsNotNone(descr)
            del descr
        
            # When this is combined with another data product, the DATA
            # array is masked with both the pixeldq and groupdq arrays.
            warnings.simplefilter("ignore")
            result = self.dataproduct + testdp
            self.assertIsNotNone(result.data)
            self.assertGreater(len(result.data), 0)
            self.assertIsNotNone(result.dq)
            self.assertGreater(len(result.dq), 0)
            descr = str(result)
            self.assertIsNotNone(descr)
            del descr
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriRampModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'refout', 'pixeldq','groupdq'],
                                       tables=['pixeldq_def', 'groupdq_def'] )
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
        
        # Attempt to access the SCI, REFOUR and DQ arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.refout)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


class TestMiriSlopeModel(unittest.TestCase):
    
    # Most of the necessary tests are already carried out by
    # the TestMiriMeasuredModel class.

    def setUp(self):
        # Create a slope data product.
        a1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
        b1 = [[1,2,3,4],     [5,6,7,8],     [9,10,11,12]]
        c1 = [[1,0,0,0],     [0,1,0,1],     [1,0,1,0]]
        acube = [a1,a1,a1]
        bcube = [b1,b1,b1]
        ccube = [c1,c1,c1]
        dcube = [a1,b1,a1]
        self.dataproduct = MiriSlopeModel(data=acube, err=bcube,
                                          dq=ccube, dq_def=master_flags,
                                          zeropt=dcube, fiterr=dcube)
        # Add some example metadata.
        self.dataproduct.set_housekeeping_metadata('MIRI EC', 'Joe Bloggs',
                                                   'V1.0')
        self.dataproduct.set_observation_metadata()
        self.dataproduct.set_target_metadata(0.0, 0.0)
        self.dataproduct.set_instrument_metadata(detector='MIRIMAGE',
                                                 filt='F2550W',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='SLOW',
                                               nints=3, ngroups=10,
                                               frame_time=1.0,
                                               integration_time=100.0,
                                               group_time=1000.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriSlopeModel_test.fits"
        
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
        # Creating an object with other than 3 dimensions must fail.
        a1d = [10,20,30,40]
        b1d = [1,2,3,4]
        c1d = [1,0,0,0]
        self.assertRaises(ValueError, MiriSlopeModel, data=a1d, err=b1d,
                          dq=c1d)

        a2d = [a1d, a1d, a1d]
        b2d = [b1d, b1d, b1d]
        c2d = [c1d, c1d, c1d]
        self.assertRaises(ValueError, MiriSlopeModel, data=a2d, err=b2d,
                          dq=c2d)

        a3d = [a2d, a2d]
        b3d = [b2d, b2d]
        c3d = [c2d, c2d]
        a4d = [a3d, a3d]
        b4d = [b3d, b3d]
        c4d = [c3d, c3d]
        self.assertRaises(ValueError, MiriSlopeModel, data=a4d, err=b4d,
                          dq=c4d)

    def test_copy(self):
        # Test that a copy can be made of the data product.
        datacopy = self.dataproduct.copy()
        self.assertIsNotNone(datacopy)
        assert_products_equal( self, self.dataproduct, datacopy,
                               arrays=['data', 'err', 'dq',
                                       'nreads', 'readsat', 'ngoodseg',
                                       'zeropt', 'fiterr'],
                               tables='dq_def' )
        del datacopy

    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with MiriSlopeModel(self.testfile) as readback:
                assert_products_equal( self, self.dataproduct, readback,
                                       arrays=['data', 'err', 'dq',
                                               'nreads', 'readsat', 'ngoodseg',
                                               'zeropt', 'fiterr'],
                                       tables='dq_def' )
                del readback
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For this test to pass these only need to run without error.
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = self.dataproduct.stats()
        self.assertIsNotNone(descr)
        del descr
        
        # Attempt to access the SCI and DQ arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.dq)
        self.assertIsNotNone(descr)
        del descr


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
