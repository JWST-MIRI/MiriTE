#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module test_exposure_model - Contains the unit tests for the classes
in the datamodels.miri_exposure_model module.

:History:

07 Mar 2013: Created.
29 Jul 2013: stats() method added.
14 Aug 2013: Updated test to include groupdq and pixeldq
20 Jun 2014: Added set_simulation_metadata.
09 Jul 2014: field_def changed to dq_def.
02 Mar 2015: Subarray burst mode added to simulator metadata.
06 Aug 2015: groupdq array removed from the constructor method.
07 Oct 2015: Made exception catching Python 3 compatible.
03 Nov 2015: Added exposure time tests.
04 May 2016: ERR array removed from ramp data model.
03 Aug 2016: Included a check on the integrity of the size of the reference
             output data.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings

import numpy as np

from miri.datamodels.miri_exposure_model import MiriExposureModel


class TestMiriExposureModel(unittest.TestCase):
    
    # Test the MiriExposureModel class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to MiriExposureModel are
    # included here.
    
    def setUp(self):
        # Create a typical exposure product.
        self.nrows = 16
        self.ncolumns = 16
        self.ngroups = 4
        self.nints = 2
        dqdata = np.zeros((self.nrows,self.ncolumns))
        self.dataproduct = MiriExposureModel(rows=self.nrows,
                                             columns=self.ncolumns,
                                             ngroups=self.ngroups,
                                             nints=self.nints, readpatt='FAST',
                                             refrows=4, refcolumns=4, grpavg=1,
                                             intavg=1, nframes=1, groupgap=0,
                                             pixeldq=dqdata, dq_def=None)
        self.dataproduct.set_simulation_metadata( 'SOLAR_MIN', burst_mode=True )
        # Give the product some dummy data.
        expdata = np.ones((self.nints,self.ngroups,self.nrows,self.ncolumns)) * 21.0
        self.dataproduct.set_exposure( expdata )
        self.testfile = "MiriExposureModel_test.fits"
        
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
        # Check that exposure data can be created with a pixel quality array.
        dqdata2d = np.zeros((self.nrows,self.ncolumns))
        dp2 = MiriExposureModel(rows=self.nrows, columns=self.ncolumns,
                                ngroups=self.ngroups, nints=self.nints,
                                readpatt='FAST', refrows=0, refcolumns=0,
                                grpavg=1, intavg=1, nframes=1,
                                groupgap=0, pixeldq=dqdata2d, dq_def=None)
        del dp2
        
        # Check that stupid settings are rejected
        # Zero rows and columns
        self.assertRaises(ValueError, MiriExposureModel, rows=0, columns=0,
                          ngroups=0, nints=0, readpatt='FAST',
                          refrows=0, refcolumns=0, grpavg=1, intavg=1,
                          nframes=1, groupgap=0, pixeldq=dqdata2d, dq_def=None)
        # Invalid data averaging.
        self.assertRaises(ValueError, MiriExposureModel, rows=self.nrows,
                          columns=self.ncolumns, ngroups=self.ngroups,
                          nints=self.nints, readpatt='FAST',
                          refrows=0, refcolumns=0, grpavg=13, intavg=13,
                          nframes=1, groupgap=0, pixeldq=dqdata2d, dq_def=None)
        # Very invalid data averaging.
        self.assertRaises(ValueError, MiriExposureModel, rows=self.nrows,
                          columns=self.ncolumns, ngroups=self.ngroups,
                          nints=self.nints, readpatt='FAST',
                          refrows=0, refcolumns=0, grpavg=0, intavg=0,
                          nframes=1, groupgap=0, pixeldq=dqdata2d, dq_def=None)
        # Invalid reference output dimensions
        self.assertRaises(ValueError, MiriExposureModel, rows=self.nrows,
                          columns=self.ncolumns, ngroups=self.ngroups,
                          nints=self.nints, readpatt='FAST',
                          refrows=5, refcolumns=3, grpavg=1, intavg=1,
                          nframes=1, groupgap=0, pixeldq=dqdata2d, dq_def=None)
        
        
    def test_set_exposure(self):
        # Test the setting of the exposure data all at once.
        # Try it with and without a GROUPDQ array.
        expdata = np.ones((self.nints,self.ngroups,self.nrows,self.ncolumns)) * 42.0
        self.dataproduct.set_exposure( expdata )
        dqdata4d = np.zeros((self.nints,self.ngroups,self.nrows,self.ncolumns))
        self.dataproduct.set_exposure( expdata, dq=dqdata4d )
        
        # Data arrays of the wrong size are rejected.
        baddata = np.ones((self.nrows,self.ncolumns)) * 42.0
        self.assertRaises(TypeError, self.dataproduct.set_exposure, baddata )
        
        # But 3-D data is converted to 4-D if it has the correct size.
        # Try it with and without a GROUPDQ array.
        # NOTE: The data3d array is defined twice because its shape is
        # mangled by the first test.
        data3d = np.ones((self.nints * self.ngroups,self.nrows,self.ncolumns)) * 16.0
        self.dataproduct.set_exposure( data3d )
        data3d = np.ones((self.nints * self.ngroups,self.nrows,self.ncolumns)) * 16.0
        dqdata3d = np.zeros((self.nints * self.ngroups,self.nrows,self.ncolumns))
        self.dataproduct.set_exposure( data3d, dq=dqdata3d )

    def test_set_integrations(self):
        # Test the setting of the data one integration at a time.
        # Try the test with and without a GROUPDQ array.
        basedata = np.ones((self.ngroups,self.nrows,self.ncolumns)) * 20.0
        dqdata3d = np.zeros((self.ngroups,self.nrows,self.ncolumns))
        for intg in range(0, self.nints):
            intdata = basedata + intg
            self.dataproduct.set_integration(intdata, intg)
        for intg in range(0, self.nints):
            intdata = basedata + intg
            self.dataproduct.set_integration(intdata, intg, dq=dqdata3d)

    def test_set_groups(self):
        # Test the setting of the data one group at a time.
        # Try the test with and without a GROUPDQ array.
        basedata = np.ones((self.nrows,self.ncolumns)) * 10.0
        dqdata2d = np.zeros((self.nrows,self.ncolumns))
        for intg in range(0, self.nints):
            for group in range(0, self.ngroups):
                grpdata = basedata + group + intg
                self.dataproduct.set_group(grpdata, group, intg)
        for intg in range(0, self.nints):
            for group in range(0, self.ngroups):
                grpdata = basedata + group + intg
                self.dataproduct.set_group(grpdata, group, intg, dq=dqdata2d)
                
    def test_averaging(self):
        # Test the generation of averaged data.
        # When grpavg=1 and intavg=1 the averaged data should be the
        # same as the original data.
        averaged = self.dataproduct.data_averaged
        self.assertIsNotNone(averaged)
        self.assertTrue( np.allclose(self.dataproduct.data, averaged) )

        # Now define a product where grpavg=4 and intavg=4.
        biggroups = self.ngroups * 4
        bigints = self.nints * 4
        avgproduct = MiriExposureModel(rows=self.nrows, columns=self.ncolumns,
                                       ngroups=biggroups, nints=bigints,
                                       readpatt='FAST', refrows=0, refcolumns=0,
                                       grpavg=4, intavg=4, nframes=1, groupgap=0,
                                       pixeldq=None, dq_def=None)
        # Give the product some dummy data.
        bigdata = np.ones((bigints,biggroups,self.nrows,self.ncolumns)) * 33.0
        avgproduct.set_exposure( bigdata )
        # In this case the averaged data will be smaller than the
        # original data
        averaged = self.dataproduct.data_averaged
        self.assertIsNotNone(averaged)
        self.assertTrue( averaged.size < avgproduct.data.size )
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS file.
            # (It isn't normally read back.)
            self.dataproduct.save(self.testfile, overwrite=True)
        
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
        
        # Attempt to access the SCI and DQ arrays through attributes.
        descr = str(self.dataproduct.data)
        self.assertIsNotNone(descr)
        del descr
        descr = str(self.dataproduct.pixeldq)
        self.assertIsNotNone(descr)
        del descr
        
    def test_exposure_times(self):
        # Test the getting and setting of exposure times
        exposure_time = 10.0
        self.dataproduct.set_exposure_times(exposure_time=exposure_time,
                                            duration=exposure_time+1.0,
                                            start_time='NOW')
        (exposure_time, duration, start_time, mid_time, end_time) = \
            self.dataproduct.get_exposure_times()
        self.assertGreater(exposure_time, 0.0, \
            "Exposure time is not greater than 0.0")
        self.assertGreaterEqual(duration, exposure_time, \
            "Duration must be at least as long as exposure time.")
        self.assertGreater(end_time, start_time, \
            "Exposure end time is not later than exposure start time")

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
