#!/usr/bin/env python

"""

Module test_spec_tools - Contains the unit tests for the module called spec_tools.

:History:

04 Feb 2014: Created
09 Dec 2014: print statements commented out. Units tests are
             supposed to generate minimal output.
08 Nov 2017: Moved from LRS pipeline package to general purpose spec_tools
             package.
04 Oct 2019: Test disabled until the following problem can be fixed.

FIXME: Test fails with the following error.

**** BAD FIT ****
Parameters were:  (0, 40.0, array([3]), 2)
Chi-Squared/dof for these parameter values = 10975.81063, CDF =    0.00000%
Uncertainties not calculated.

Try a different initial guess for the fit parameters.
Or if these parameters appear close to a good fit, try giving
    the fitting program more time by increasing the value of maxfev.
scipy.optimize.curve_fit failed

10 Oct 2019: Test enabled again.

@author:  Juergen Schreiber

"""
# This module is now converted to Python 3.


import unittest
import numpy as np

from miri.tools.spec_tools import condense, convGauss
from miri.tools.spec_tools import lrs2D_spextract, lrs_extract_spec, \
   lrs_extract_spec_with_fit, get_psf_fit, optimalSpecExtraction, \
   subtractBackground, interpolWaveOnRows, interpolLin, interpolSpline

from miri.datamodels.miri_measured_model import MiriMeasuredModel


class TestCondense(unittest.TestCase):
    def setUp(self):
        # Create an instance of Class1 and prepare it to be tested.
        # The setUp function is automatically executed before each test.
        self.spec = np.arange(50, 100, 0.2)
        self.x = np.linspace(0, 49)


    def tearDown(self):
        # Clean the files and objects created.
        # The tearDown function is automatically executed after each test.
        del self.spec
        del self.x
        del self.result

    
    def test_condense(self):
        # Test the method Class1.method1
        # See http://docs.python.org/library/unittest.html#test-cases for further
        # information
#         print "run condense"
        self.result = condense(self.spec, 5)
        self.assertTrue(len(self.spec)== len(self.result))
        self.result = condense(self.spec, 10, method = "Mean",running = False)
        self.assertTrue (len(self.spec)/10 == len(self.result))
        
        
    def test_convGauss(self):
#         print "run convGauss"
        self.result = convGauss(self.x, 5, 1)
 
        self.assertTrue ( np.mean(self.x) == np.mean(self.result))
        self.assertTrue( len(self.result) == len(self.x) - 4)
class TestLrsExtract(unittest.TestCase):
    def setUp(self):
        # Create an instance of Class1 and prepare it to be tested.
        # The setUp function is automatically executed before each test.
        a = [[10, 10, 20, 40,  20, 10, 10],
         [10, 10, 20, 50,  20, 10, 10],
         [10, 10, 20, 150,  20, 10, 10],
         [10, 10, 20, 70,  50, 25, 10],
         [10, 10, 20, 50,  20, 10, 10],
         [10, 10, 40, 70,  30, 10, 10],
         [10, 10, 20, 170,  20, 10, 10],
         [10, 10, 20, 60,  20, 10, 10],
         [10, 10, 20, 40,  20, 10, 10]]
        err = [[0.5, 0.5, 0.2, 0.05, 0.2, 0.1, 0.5],
         [0.6, 0.5, 0.2, 0.05, 0.2, 0.5, 0.6],
         [0.6, 0.5, 0.2, 0.05, 0.2, 0.5, 0.6],
         [0.5, 0.4, 0.2, 0.05, 0.2, 0.4, 0.5],
         [0.4, 0.5, 0.2, 0.02, 0.2, 0.5, 0.5],
         [0.1, 0.4, 0.2, 0.02, 0.2, 0.6, 0.4],
         [0.6, 0.5, 0.2, 0.03, 0.2, 0.5, 0.6],
         [0.6, 0.5, 0.2, 0.02, 0.2, 0.5, 0.5],
         [0.5, 0.5, 0.2, 0.05, 0.2, 0.5, 0.4]]
        self.dataproduct = MiriMeasuredModel(data = a, err = err)
        self.assertTrue(hasattr(self.dataproduct,'data'))

    def tearDown(self):
        # Clean the files and objects created.
        # The tearDown function is automatically executed after each test.
        del self.dataproduct
        
    def test_lrs_extract(self):
        # Test the method Class1.method1
        # See http://docs.python.org/library/unittest.html#test-cases for further
        # information
#         print "run lrs2D_spextract"
        prod = self.dataproduct.copy()
        result1 = lrs2D_spextract(prod, xmin = 0, xmax = 2, ymin = 0, ymax = 2)

        # All the assert* methods available in the module unittest.TestCase for
        # python >=2.5 are listed below with a usage example
        self.assertTrue(np.shape(result1.data)[0] ==2)
        self.assertTrue(np.shape(result1.data)[1] ==2)
#         print "run subtractBackground"
        sub = subtractBackground(result1, result1)
        self.assertTrue(sub.data[0,0] == 0.)
#         print "run lrs_extract_spec"
        spec = lrs_extract_spec(sub)
        self.assertTrue(len(spec.data)== 2)
        self.assertTrue(len(spec.err)== 2)
#         print "run lrs_extract_spec_with_fit"
        spec1 = lrs_extract_spec_with_fit(prod.copy())
        self.assertTrue( len(spec1) == 3)
#         print "run optimalSpecExtract"
        spec1 = optimalSpecExtraction(prod.copy())
        self.assertTrue( len(spec1.data) == 9)        
        self.assertTrue( len(spec1.err) == 9)  
#         print "run get_psf_fit"
        spec1 = get_psf_fit(prod.copy())
        self.assertTrue( len(spec1) == 4)
        
#         print "run interpolWaveOnRows"
        wave = np.arange(5, 10, 0.5)
        pos = np.arange(0, 5, 0.5)
        rows, new_wave = interpolWaveOnRows(pos, wave)
        self.assertTrue(len(rows)==5)
        self.assertTrue(len(new_wave)==5)
        
        self.assertTrue(np.alltrue(np.equal(np.mod(new_wave,1),0)))
        
#         print "run interpolLin"
        spec = np.arange(100, 105, 0.5)
        new_spec = interpolLin(wave, spec, new_wave)
        self.assertTrue(len(new_spec)==5)
        self.assertTrue(np.alltrue(np.equal(np.mod(new_spec,1),0)))
        
#         print "run interpolSpline"
        spec = np.arange(100, 105, 0.5)
        new_spec = interpolSpline(wave, spec, new_wave)
        self.assertTrue(len(new_spec)==5)             
        
        del result1, spec, spec1, sub, prod, new_spec

if __name__ == '__main__':
    unittest.main()

