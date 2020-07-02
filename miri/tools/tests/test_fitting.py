#!/usr/bin/env python

"""

Module test_fitting - Contains the unit tests for the module called fitting.

:History:

09 Jan 2013: Created
02 Sep 2013: Debugged. Suppress unwanted output.
08 Nov 2017: Moved from LRS pipeline package to general purpose fitting
             package.

@author:  Juergen Schreiber

"""
# This module is now converted to Python 3.


import unittest
import numpy

from miri.tools.fitting import gaussian, gaussPlusPoly, nonLinFit, LinFit


class TestFitting(unittest.TestCase):
    def setUp(self):
        # Create an instance of Class1 and prepare it to be tested.
        # The setUp function is automatically executed before each test.
        self.x_data = numpy.arange(0,100)
        self.p = (1., 100., 50., 10.)
        self.y_data = gaussian(self.x_data, *self.p)
        self.p_guess = (2, 110, 40, 7)
        self.p_poly = (1., 100., 50., 10., 1., 1.)
        self.p_poly_guess = (2, 149, 40, 7, 1.5, 1.5)
        self.y_data_poly = gaussPlusPoly(self.x_data, *self.p_poly)
        self.fitter = LinFit(self.x_data, self.x_data)

    def tearDown(self):
        # Clean the files and objects created.
        # The tearDown function is automatically executed after each test.
        del self.x_data
        del self.y_data
        del self.p
        del self.p_guess
        del self.p_poly
        del self.p_poly_guess
        del self.y_data_poly
        del self.fitter
 
    
    def test_fit(self):
        # Test the method Class1.method1
        # See http://docs.python.org/library/unittest.html#test-cases for further
        # information
        #print "run gaussian fit"
        #print gaussian
        fit, popt, perr = nonLinFit(gaussian, self.x_data, self.y_data, p_guess=self.p_guess, verbose=False)
        diff = max(fit) - (self.p[1]+self.p[0])
        self.assertTrue(abs(diff) < 0.0001)
        for i in range(len(self.p)):
            self.assertAlmostEqual(self.p[i],popt[i])
        
        
        #print "run gauss + polynomial fit"
        fitPoly, poptPoly, perr = nonLinFit(gaussPlusPoly, self.x_data, self.y_data_poly, p_guess=self.p_poly_guess, verbose=False) 

        for i in range(len(self.p)):
            self.assertAlmostEqual(self.p_poly[i],poptPoly[i])

        
        coeff = self.fitter.fit()
        self.assertTrue(coeff[0] == 1.)
        
        error = self.fitter.calcLinFitError(self.x_data)
        self.assertTrue(numpy.mean(error) == 0.)    
        line = self.fitter.straightLine(self.x_data)
        self.assertTrue(numpy.all(line - self.x_data) == 0.)    
        
        del coeff
        del error
        del fitPoly
        del popt
        del poptPoly
        del fit
        del line
        
if __name__ == '__main__':
    unittest.main()

