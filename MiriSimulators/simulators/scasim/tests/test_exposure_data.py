#!/usr/bin/env python

"""

Module test_exposure_data - Contains the unit tests for the
ExposureData class.
NOTE: Since this module depends on the quantum_efficiency module, that module
should be executed first.

:History:

08 Sep 2015: Created
03 Nov 2015: Added exposure time tests.
09 Mar 2016: Make the message about exceeding the data size limit a
             warning rather than an exception.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os, time
import unittest
import numpy as np

from miri.simulators.scasim.exposure_data import ExposureData

class TestExposureData(unittest.TestCase):
     
    def setUp(self):
        # Create an ExposureData object for a 60 row x 80 column detector
        self.rows = 60
        self.columns = 80
        self.nints = 2
        self.ngroups = 2
        readpatt = 'FAST'
        self.exposure_data = ExposureData(self.rows, self.columns, self.ngroups,
                                          self.nints, readpatt)
          
    def tearDown(self):
        # Tidy up
        del self.exposure_data
     
    def test_creation(self):
        # Attempt to create a data set with a large number of integrations
        # and check that a warning is issued.
#         self.rows = 1024
#         self.columns = 1024
#         self.nints = 100
#         self.ngroups = 10
#         readpatt = 'FAST'
#         exposure_data = ExposureData(self.rows, self.columns, self.ngroups,
#                                           self.nints, readpatt)        
#         # Generate and apply an exposure in one go.
#         data = np.ones([self.nints, self.ngroups, self.rows, self.columns])
#         exposure_data.set_exposure(data)
        pass
     
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        descr = self.exposure_data.__str__()
        self.assertIsNotNone(descr)

    def test_metadata(self):
        # Test the exposure data metadata functions
        self.exposure_data.set_exposure_times(10.0, 12.0, 'NOW')

    def test_groups(self):
        # Generate a series of integrations and groups with predictable data values
        # and build up the exposure data.    
        for intg in range(0,self.nints):
            for group in range(0,self.ngroups):
                data = 1 * np.indices([self.rows,self.columns])[0] * \
                    (group+1) * (intg+2)
                self.exposure_data.set_group(data, group, intg)
                del data
                
        slope_data = self.exposure_data.slope_data(grptime=2.785)
        self.assertIsNotNone(slope_data)
        del slope_data

    def test_exposure(self):
        # Generate an exposure in one go.    
        data = np.ones([self.nints, self.ngroups, self.rows, self.columns])
        self.exposure_data.set_exposure(data)
        del data
        
    def test_exposure_times(self):
        # Test the getting and setting of exposure times
        exposure_time = 10.0
        self.exposure_data.set_exposure_times(exposure_time=exposure_time,
                                            duration=exposure_time+1.0,
                                            start_time='NOW')
        (exposure_time, duration, start_time, mid_time, end_time) = \
            self.exposure_data.get_exposure_times()
        self.assertGreater(exposure_time, 0.0, \
            "Exposure time is not greater than 0.0")
        self.assertGreaterEqual(duration, exposure_time, \
            "Duration must be at least as long as exposure time.")
        self.assertGreater(end_time, start_time, \
            "Exposure end time is not later than exposure start time")
        

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
