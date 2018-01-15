#!/usr/bin/env python

"""

Module test_measured_variable - Contains the unit tests for the
MiriMeasurement class.

:History:

11 Jul 2014: Created to replace the test_measured_variable.py module,
             which was based on the old MeasuredVariable class.
08 Sep 2015: Made compatible with Python 3
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os
import unittest
import warnings
import numpy as np

from miri.datamodels.tests.util import assert_recarray_equal, \
    assert_products_equal
from miri.datamodels.miri_measurement import MiriMeasurement

_NUMBER_OF_MEASUREMENTS = 10

class TestMiriMeasurement(unittest.TestCase):
    
    def setUp(self):
        # Set up a linearly spaced arrays of _NUMBER_OF_MEASUREMENTS
        # floating point parameter values from 1.0 to 10.0.
        self.parameters = np.linspace(1.0, 10.0, _NUMBER_OF_MEASUREMENTS)
        
        # Set up a linearly spaced array of corresponding variable
        # measurements from 1.0 to 10.0. Also set up an exponential
        # set of measurements.
        self.vallin = np.linspace(1.0, 10.0, _NUMBER_OF_MEASUREMENTS)
        self.valexp = 10.0 ** self.vallin
        
        measurement_table_lin = []
        measurement_table_exp = []
        for ii in range(0, len(self.parameters)):
            recordlin = (self.parameters[ii], self.vallin[ii], 1.0, 1.0, 1.0, 1.0)
            recordexp = (self.parameters[ii], self.valexp[ii], 1.0, 1.0, 1.0, 1.0)
            measurement_table_lin.append(recordlin)
            measurement_table_exp.append(recordexp)

        # Create two kinds of MiriMeasurement object: A linearly
        # interpolated object containing the linear values and a
        # logarithmically interpolated object containing the
        # exponential values.
        self.mv_lin = MiriMeasurement(measurement_table=measurement_table_lin,
                                      title='LINEAR test object',
                                      name='Widgets', vname='Thingies',
                                      interptype='LINEAR')
        self.mv_log = MiriMeasurement(measurement_table=measurement_table_exp,
                                      title='LOGLIN test object',
                                      name='Widgets', vname='Thingies',
                                       interptype='LOGLIN')

        # Parameters for test_tofits and test_fromfits
        self.tempfile1 = 'tmp_lin.fits'
        self.tempfile2 = 'tmp_exp.fits'
   
    def tearDown(self):
        # Clean temporary files.
        if os.path.isfile(self.tempfile1):
            try:
                os.remove(self.tempfile1)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.tempfile1 + \
                        "\n   " + str(e)
                warnings.warn(strg)
        if os.path.isfile(self.tempfile2):
            try:
                os.remove(self.tempfile2)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.tempfile2 + \
                        "\n   " + str(e)
                warnings.warn(strg)

        # Tidy up
        del self.mv_lin
        del self.mv_log
        del self.valexp
        del self.vallin
        del self.parameters
        
    def test_creation(self):
        pass # TBD
#         # Check for exceptions when creating bad MiriMeasurement objects.
#         # Both arrays must be the same size.
#         params = [1,2,3,4,5]
#         vals = [1,2,3]
#         self.assertRaises(TypeError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '')
# 
#         # Specifying a 2-D values array should be fine, as long as the
#         # size of each array listed in the values array matches the
#         # parameters array, and each array listed has a matching column
#         # name and column units.
#         vals = [[1,2,3,4,5],[6,7,8,9,10]]
#         temp = MiriMeasurement(params, 'param', '', vals, ['val1','val2'],
#                                 ['unit1','unit2'])
#         del temp
#         # It is possible to give all the columns the same unit by
#         # specifying a single string.
#         temp = MiriMeasurement(params, 'param', '', vals, ['val1','val2'],
#                                 'unit')
#         del temp
#         # An exception should be raised if the values array has the wrong
#         # shape compared to params.
#         vals = [[1,2,3,4],[6,7,8,9]]
#         self.assertRaises(TypeError, MiriMeasurement,
#                           params, 'param', '', vals, ['val1','val2'],
#                           ['',''])
#         # A non-rectangular values array should also raise an exception.
#         vals = [[1,2,3,4],[6,7,8,9,10]]
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, ['val1','val2'],
#                           ['',''])
# 
#         # Providing a string instead of a list of strings for the column
#         # names should be ok. The object will make its own unique column
#         # names.
#         vals = [[1,2,3,4,5],[6,7,8,9,10]]
#         temp = MiriMeasurement(params, 'param', '', vals, 'val', '')
#         del temp
#         # The column names and column units lists must be the same length.
#         self.assertRaises(TypeError, MiriMeasurement,
#                           params, 'param', '', vals, ['val1', 'val2'],
#                           ['unit2'])
# 
#         # The parameters array must be monotonically increasing
#         # and no value should appear more than once.
#         params = [1,2,4,3,5]
#         vals = [1,2,3,4,5]
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '')        
#         params = [1,2,4,4,5]
#         vals = [1,2,3,4,5]
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '')        
# 
#         # Single values or garbage for the parameters or values arrays
#         # should raise an exception.
#         self.assertRaises(TypeError, MiriMeasurement,
#                           42, 'param', '', 1, 'val', 'unit')
#         params = [1,2,3,4,5]
#         self.assertRaises(TypeError, MiriMeasurement,
#                           params, 'param', '', 1, 'val', 'unit')
#         self.assertRaises(TypeError, MiriMeasurement,
#                           [], 'param', '', [], 'val', '')
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', 'forty two', 'val', '')
# 
#         # Logarithmic interpolation cannot be used if the corresponding
#         # array contains negative values.
#         params = [1,2,3,4,5]
#         vals = [-1,0,1,2,3]
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '',
#                           interptype='LOGLIN')
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '',
#                           interptype='LOGLOG')
# 
#         params = [-1,0,1,2,3]
#         vals = [1,2,3,4,5]
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '',
#                           interptype='LINLOG')
#         self.assertRaises(ValueError, MiriMeasurement,
#                           params, 'param', '', vals, 'val', '',
#                           interptype='LOGLOG')
# 
#         # Linear interpolation with negative values should be ok
#         params = [-1,0,1,2,3]
#         vals = [-2,-1,0,1,2]
#         temp = MiriMeasurement(params, 'param', '', vals, 'val', '',
#                           interptype='LINEAR')
#         result = temp.lookup(-0.5)
#         expected = -1.5
#         self.assertAlmostEqual(expected, result)
#         del temp


    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these need to run without error
        # and generate non-null strings.
        descr1 = self.mv_lin.__str__()
        descr2 = self.mv_log.__str__()
        self.assertIsNotNone(descr2)
        self.assertIsNotNone(descr2)

    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
            # Check that the data products can be written to a FITS
            # file and read back again without changing the data.
            self.mv_lin.save(self.tempfile1, overwrite=True)
            with MiriMeasurement(self.tempfile1) as readback1:
                assert_products_equal( self, self.mv_lin, readback1,
                                       arrays=[], tables='measurement_table' )
                del readback1
            self.mv_lin.save(self.tempfile2, overwrite=True)
            with MiriMeasurement(self.tempfile2) as readback2:
                assert_products_equal( self, self.mv_lin, readback2,
                                       arrays=[], tables='measurement_table' )
                del readback2

    def test_lookup(self):
        # Take each element from the parameter array and look it up in
        # both of the MiriMeasurement objects. The value returned must
        # be the same as the corresponding element in the values array.
        for ii in range(0,len(self.parameters)):
            # Linear
            pvalue = self.parameters[ii]
            expected = self.vallin[ii]
            result = self.mv_lin.lookup(pvalue)
            self.assertEqual(expected, result)
            # Logarithmic
            expected = self.valexp[ii]
            result = self.mv_log.lookup(pvalue)
            # Use AlmostEqual because logarithmic interpolation involves
            # floating point calculations, even at the defined points.
            self.assertAlmostEqual(expected, result)
            
        # An attempt to look up a column outside the normal range
        # should raise an exception.
        ncolumns = self.mv_lin.meta.measurement_table.mcolumns
        pvalue = self.parameters[0]
        result = self.mv_lin.lookup(pvalue, column=0)
        self.assertRaises(ValueError, self.mv_lin.lookup, pvalue,
                          column=-1)
        self.assertRaises(ValueError, self.mv_lin.lookup, pvalue,
                          column=ncolumns)
            
    def test_interpolation(self):
        # Test that the interpolation works as expected when looking
        # up values that don't correspond to the measured points.
        # This time look up values exactly half way between measurements.
        logvals = np.log10(self.valexp)
        for ii in range(0,len(self.parameters)-1):
            # With the linear interpolation each result should also be half
            # way between the corresponding values.
            pvalue = (self.parameters[ii] + self.parameters[ii+1])/2.0
            expected = (self.vallin[ii] + self.vallin[ii+1])/2.0
            result = self.mv_lin.lookup(pvalue)
            self.assertAlmostEqual(expected, result)
            # With logarithmic interpolation the result should be half
            # way between the logarithms of the corresponding values.
            logexpected = (logvals[ii] + logvals[ii+1])/2.0
            expected = 10.0 ** logexpected
            result = self.mv_log.lookup(pvalue)
            self.assertAlmostEqual(expected, result)
            
    def test_extremes(self):
        # Looking up values beyond the range of measurements simply
        # returns the first or last measurement.
        toosmall = self.parameters.min() - 10.0
        toolarge = self.parameters.max() + 10.0
        # Linear
        expected = self.vallin[0]
        result = self.mv_lin.lookup(toosmall)
        self.assertEqual(expected, result)
        expected = self.vallin[-1]
        result = self.mv_lin.lookup(toolarge)
        self.assertEqual(expected, result)
        # Logarithmic
        expected = self.valexp[0]
        result = self.mv_log.lookup(toosmall)
        self.assertAlmostEqual(expected, result)
        expected = self.valexp[-1]
        result = self.mv_log.lookup(toolarge)
        self.assertAlmostEqual(expected, result)


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
