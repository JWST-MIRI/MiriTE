#!/usr/bin/env python

"""

Module test_integrators - Contains the unit tests for the
PoissonIntegrator and ImperfectIntegrator class.

:History:

25 Aug 2010: Created
03 Sep 2010: Test the new leak function.
27 Sep 2010: Python environment for windows verified.
16 Nov 2010: Set the seed of the random number generator, for more
             predictable and repeatable tests.
05 Oct 2011: Enhanced the saturation test to check extreme input.
13 Nov 2012: Major restructuring of the package folder. Import
             statements updated.
23 Apr 2014: Changed the terminology to match common MIRI usage (e.g.
             an "integration" includes an interval between resets, not
             a fragment of that interval). Removed redundant methods.
             ImperfectIntegrator class added to separate effects
             due to other than Poisson statistics.
17 Jun 2014: Removed the old "preflash" latency implementation.
08 Sep 2015: Make compatible with Python 3.
04 Dec 2015: Renamed from poisson_integrator to integrators.
12 Feb 2016: Added tests for extremely small floating point flux values.
17 Feb 2017: The Poisson noise calculation no longer includes a ratchet,
             so successive values can go down as well as up. Test readouts
             converted to signed integers so that negative differences
             can be managed.
05 May 2017: Corrected permission for nosetests.
13 Dec 2017: Added flux and noise tests.
04 Jan 2017: Check the flux is correct when noise is turned on and when there
             is a non-zero pedestal. Also check the flux is correct when
             zeropoint drift and latency effects are included.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

#import os, sys
import copy
import unittest
import numpy as np

from miri.simulators.integrators import PoissonIntegrator, ImperfectIntegrator

class TestPoissonIntegrator(unittest.TestCase):
    
    def setUp(self):
        # Create a very simple 3 x 3 Poisson integrator object
        self.integrator = PoissonIntegrator(3, 3, verbose=0)
        self.integrator.set_seed(42)
        
    def tearDown(self):
        # Tidy up
        del self.integrator

    def test_creation(self):
        # Check for exceptions when creating bad PoissonIntegrator objects.
        # The dimensions of the object must be positive.
        self.assertRaises(ValueError, PoissonIntegrator, -1, -1, verbose=0)

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        title = self.integrator.get_title()
        self.assertIsNotNone(title)
        self.assertIsNotNone(self.integrator.nperiods)
        self.assertIsNotNone(self.integrator.readings)
        self.assertIsNotNone(self.integrator.exposure_time)
        # The shape should be 3x3, as created.
        self.assertEqual(self.integrator.shape[0], 3)
        self.assertEqual(self.integrator.shape[1], 3)
        counts = self.integrator.get_counts()
        del counts
        descr = self.integrator.__str__()
        self.assertIsNotNone(descr)
        
    def test_integrate(self):
        # The readings must never decrease following an integration,
        # no matter how short the integration time.
        # NOTE: This is no longer true. Tests commented out.
        flux = [[1.0, 2.0, 3.0], \
                [4.0, 5.0, 6.0], \
                [7.0, 8.0, 9.0]]
        self.integrator.reset()
        self.integrator.integrate(flux, 100.0)
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 100.0)
        readout2 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 1.0)
        readout3 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 0.0)
        readout4 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
#        diff2 = readout3 - readout2
#        diff3 = readout4 - readout3
        self.assertTrue(diff1.min() >= 0)
#        self.assertTrue(diff2.min() >= 0)
#        self.assertTrue(diff3.min() >= 0)
        # The final exposure time must be the sum of all the
        # integration times.
        expected = 100.0 + 100.0 + 1.0 + 0.0
        actual = self.integrator.exposure_time
        self.assertAlmostEqual(expected, actual)
        # An integration with no flux, or a wait, does not
        # increases the exposure time.
        before = self.integrator.exposure_time
        self.integrator.integrate(None, 100.0)
        after = self.integrator.exposure_time
        self.assertAlmostEqual(before, after)
        before = self.integrator.exposure_time
        self.integrator.wait(100.0)
        after = self.integrator.exposure_time
        self.assertAlmostEqual(before, after)
        # A perfect reset should restore the readings to zero
        self.integrator.reset()
        readout6 = self.integrator.readout()
        self.assertTrue(np.all(readout6 == 0))

    def test_flux(self):
        # Test that, when noise is switched off, the output flux is the
        # input illumination multiplied by the exposure time, at least
        # to within one electron. Defining a pedestal level should make
        # no difference.
        # Create a very simple 3 x 3 Poisson integrator object with
        # Poisson noise turned off.
        intnonoise = PoissonIntegrator(3, 3, simulate_poisson_noise=False,
                                       verbose=0)
        intnonoise.set_pedestal(8000 * np.ones([3, 3]))
        exptime = 100.0
        for fluxlevel in (1.234, 12.345, 123.456, 1234.56):
            flux = fluxlevel * np.ones((3,3), dtype=np.float32)
            intnonoise.reset()
            intnonoise.integrate(flux, exptime)
            readout1 = intnonoise.readout().astype(np.int32)
            intnonoise.integrate(flux, exptime)
            readout2 = intnonoise.readout().astype(np.int32)
            difference = readout2 - readout1
            # The flux difference must be less than 1 electron.
            deviation = difference.mean() - (fluxlevel * exptime)
            self.assertLess(abs(deviation), 1.0)
        del intnonoise
        # Check that the flux level is also approximately
        # correct when noise is turned on
        self.integrator.set_pedestal(8000 * np.ones([3, 3]))
        exptime = 100.0
        for fluxlevel in (1.234, 12.345, 123.456, 1234.56):
            flux = fluxlevel * np.ones((3,3), dtype=np.float32)
            self.integrator.reset()
            self.integrator.integrate(flux, exptime)
            readout1 = self.integrator.readout().astype(np.int32)
            self.integrator.integrate(flux, exptime)
            readout2 = self.integrator.readout().astype(np.int32)
            difference = readout2 - readout1
            # The flux difference must be less than 5%.
            deviation = difference.mean() - (fluxlevel * exptime)
            self.assertLess(abs(deviation), difference.mean()/20.0)

    def test_noise(self):
        # Test that the Poisson noise level is approximately the
        # square root of the input signal (Bug 439).
        # The test is done using a 1024x1024 integrator to
        # eliminate the small number statistics.
        intnoise = PoissonIntegrator(1024, 1024, verbose=0)
        intnoise.set_seed(42)
        exptime = 10.0
        # Try 2 groups with varying flux levels.
        for fluxlevel in (1.234, 12.345, 123.456, 1234.56):
            flux = fluxlevel * np.ones((1024,1024), dtype=np.float32)
            intnoise.reset()
            readout1 = intnoise.readout().astype(np.int32)
            intnoise.integrate(flux, exptime)
            readout2 = intnoise.readout().astype(np.int32)
            difference = readout2 - readout1
            noise = difference.std()
            # The following comparison can only be approximate, since the
            # noise should approximate the square root of expected count,
            # not the actual readings. Make sure it is within 1%.
            deviation = difference.mean() -  (noise * noise)
            self.assertLess(abs(deviation), difference.mean()/100.0)
        # Try 10 groups with a fixed flux level.
        fluxlevel = 4.0
        flux = fluxlevel * np.ones((1024,1024), dtype=np.float32)
        intnoise.reset()
        last_readout = intnoise.readout().astype(np.int32)
        first_readout = copy.deepcopy(last_readout)
        for group in range(1, 10):
            intnoise.integrate(flux, exptime)
            this_readout = intnoise.readout().astype(np.int32)
            difference = this_readout - last_readout
            noise = difference.std()
            # Make sure the noise is within 1% of expectation.
            deviation = difference.mean() -  (noise * noise)
#             print("Mean=",difference.mean(), "Noise^2=", (noise*noise))
#             print("Compare", abs(deviation), "with", difference.mean()/100.0)
            self.assertLess(abs(deviation), difference.mean()/100.0)
            last_readout = this_readout
        # Check the noise level between the last and first reading.
        difference = this_readout - first_readout
        noise = difference.std()
        # Make sure the noise is within 1% of expectation.
        deviation = difference.mean() -  (noise * noise)
#         print("Mean=",difference.mean(), "Noise^2=", (noise*noise))
#         print("Compare", abs(deviation), "with", difference.mean()/100.0)
        self.assertLess(abs(deviation), difference.mean()/100.0)
        del intnoise

    def test_saturation(self):
        # If a bucket size is defined, the integrator will saturate
        # when the count reaches or exceeds this bucket size.
        saturated = 10
        flux = [[100.0, 200.0, 300.0], \
                [400.0, 500.0, 600.0], \
                [700.0, 800.0, 900.0]]
        test = PoissonIntegrator(3, 3, bucket_size=saturated, verbose=0)
        test.reset()
        test.integrate(flux, 1000.0)
        readout = test.readout()
        self.assertTrue(np.all(readout == saturated))
        del test
        
        # If a bucket size is not defined, the integrator must still
        # be able to cope with extreme input without overflowing.
        flux = [[1.0e9, 2.0e9, 3.0e9], \
                [4.0e9, 5.0e9, 6.0e9], \
                [7.0e9, 8.0e9, 9.0e9]]
        test = PoissonIntegrator(3, 3, bucket_size=None, verbose=0)
        test.reset()
        test.integrate(flux, 1000.0)
        readout = test.readout()
        # Overflow will cause spurious negative values.
        self.assertTrue(np.all(readout >= 0))
        del test
        
    def test_extreme(self):
        # Check that the integration and readout functions can accept
        # extremely small floating point values without raising an
        # exception (Bug 16).
        self.integrator.reset()
        flux = [[1.0e-9, 2.0e-9, -3.0e9], \
                [4.0e-9, 5.0e-15, 6.0e-9], \
                [7.0e-9, -8.0e-4, -9.0e-15]]
        self.integrator.integrate(flux, 1.0)
        readout = self.integrator.readout(nsamples=4)
        # Check there are no negative values.
        self.assertTrue(np.all(readout >= 0))
        # Check that the functions can safely accept a stream of random
        # positive and negative values.
        for rep in range(0,10):
            randflux = np.random.randn( 3, 3 ) - 0.5
            self.integrator.integrate(randflux, 10.0)
            readout = self.integrator.readout(nsamples=4)
            # Check there are no negative values.
            self.assertTrue(np.all(readout >= 0))

    def test_wrong_shape(self):
        # Attempting to integrate on a flux array of the wrong shape
        # should raise an exception.
        flux = [1.0, 3.0, 5.0]
        self.assertRaises(TypeError, self.integrator.integrate, flux, 1.0)

    def test_negative_time(self):
        # Attempting to integrate with a negative integration time
        # should raise an exception.
        flux = [[1.0, 2.0, 3.0], \
                [4.0, 5.0, 6.0], \
                [7.0, 8.0, 9.0]]
        self.assertRaises(ValueError, self.integrator.integrate, flux, -1.0)


class TestImperfectIntegrator(unittest.TestCase):
    
    def setUp(self):
        # Create a very simple 3 x 3 Poisson integrator object
        self.integrator = ImperfectIntegrator(3, 3, verbose=0)
        self.integrator.set_seed(42)
        
    def tearDown(self):
        # Tidy up
        del self.integrator

    def test_creation(self):
        # Check for exceptions when creating bad ImperfectIntegrator objects.
        # TBD
        pass

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        title = self.integrator.get_title()
        self.assertIsNotNone(title)
        descr = self.integrator.__str__()
        self.assertIsNotNone(descr)
        
    def test_integrate(self):
        # The readings must never decrease following an integration,
        # no matter how short the integration time.
        # NOTE: This is no longer true. Tests commented out.
        flux = [[1.0, 2.0, 3.0], \
                [4.0, 5.0, 6.0], \
                [7.0, 8.0, 9.0]]
        self.integrator.reset()
        self.integrator.integrate(flux, 100.0)
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 100.0)
        readout2 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 1.0)
        readout3 = self.integrator.readout().astype(np.int32)
        self.integrator.integrate(flux, 0.0)
        readout4 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
#        diff2 = readout3 - readout2
#        diff3 = readout4 - readout3
        self.assertTrue(diff1.min() >= 0)
#        self.assertTrue(diff2.min() >= 0)
#        self.assertTrue(diff3.min() >= 0)
        # The final exposure time must be the sum of all the
        # integration times.
        expected = 100.0 + 100.0 + 1.0 + 0.0
        actual = self.integrator.exposure_time
        self.assertAlmostEqual(expected, actual)
        # An integration with no flux, or a wait, does not
        # increases the exposure time.
        before = self.integrator.exposure_time
        self.integrator.integrate(None, 100.0)
        after = self.integrator.exposure_time
        self.assertAlmostEqual(before, after)
        before = self.integrator.exposure_time
        self.integrator.wait(100.0)
        after = self.integrator.exposure_time
        self.assertAlmostEqual(before, after)
        # Defining a persistence will make the reset
        # leave behind some a residual count.
        self.integrator.set_persistence(0.1)
        self.integrator.reset()
        readout6 = self.integrator.readout()
        self.assertTrue(np.all(readout6 > 0))
        # Turning off the persistence allows a perfect reset.
        self.integrator.set_persistence(0.0)
        self.integrator.reset()
        readout7 = self.integrator.readout()
        self.assertTrue(np.all(readout7 == 0))
 
    def test_flux(self):
        # Check that the flux level is still approximately correct when noise,
        # zeropoint drifts and latency effects are all turned on
        self.integrator.set_pedestal(8000 * np.ones([3, 3]))
        zp_slow = [45000.0, 0.0084]  # Slow zeropoint drift parameters
        zp_fast = [[0.0, -2.917],    # Fast zeropoint drift parameters
                   [0.0, -2.292],
                   [0.0, -2.396],
                   [0.0, -2.408]]
        self.integrator.set_zeropoint(zp_slow, zp_fast)
        slow_latency = [1.67e-9, 136000.0] # Slow latency parameters
        fast_latency = [0.002, 300.0]      # Fast latency parameters
        self.integrator.set_latency(slow_latency, fast_latency)
        sensitivity = [1.1, -0.4]  # Linearity sensitivity coeffs
        self.integrator.set_sensitivity(sensitivity)
        exptime = 100.0
        for fluxlevel in (1.234, 12.345, 123.456, 1234.56):
            flux = fluxlevel * np.ones((3,3), dtype=np.float32)
            self.integrator.reset()
            self.integrator.integrate(flux, exptime)
            readout1 = self.integrator.readout().astype(np.int32)
            self.integrator.integrate(flux, exptime)
            readout2 = self.integrator.readout().astype(np.int32)
            difference = readout2 - readout1
            # The flux difference must be less than 5%.
            deviation = difference.mean() - (fluxlevel * exptime)
            self.assertLess(abs(deviation), difference.mean()/20.0)
        
    def test_saturation(self):
        # If a bucket size is defined, the integrator will saturate
        # when the count reaches or exceeds this bucket size.
        saturated = 10
        flux = [[100.0, 200.0, 300.0], \
                [400.0, 500.0, 600.0], \
                [700.0, 800.0, 900.0]]
        test = ImperfectIntegrator(3, 3, bucket_size=saturated, verbose=0)
        test.reset()
        test.integrate(flux, 1000.0)
        readout = test.readout()
        self.assertTrue(np.all(readout == saturated))
        del test
        
        # If a bucket size is not defined, the integrator must still
        # be able to cope with extreme input without overflowing.
        flux = [[1.0e9, 2.0e9, 3.0e9], \
                [4.0e9, 5.0e9, 6.0e9], \
                [7.0e9, 8.0e9, 9.0e9]]
        test = ImperfectIntegrator(3, 3, bucket_size=None, verbose=0)
        test.reset()
        test.integrate(flux, 1000.0)
        readout = test.readout()
        # Overflow will cause spurious negative values.
        self.assertTrue(np.all(readout >= 0))
        del test

    def test_extreme(self):
        # Check that the integration and readout functions can accept
        # extremely small floating point values without raising an
        # exception.
        self.integrator.reset()
        flux = [[1.0e-9, 2.0e-9, -3.0e9], \
                [4.0e-9, 5.0e-15, 6.0e-9], \
                [7.0e-9, -8.0e-4, -9.0e-15]]
        self.integrator.integrate(flux, 1.0)
        readout = self.integrator.readout(nsamples=4)
        # Check there are no negative values.
        self.assertTrue(np.all(readout >= 0))
        # Check that the functions can safely accept a stream of random
        # positive and negative values.
        for rep in range(0,10):
            randflux = np.random.randn( 3, 3 ) - 0.5
            self.integrator.integrate(randflux, 10.0)
            readout = self.integrator.readout(nsamples=4)
            # Check there are no negative values.
            self.assertTrue(np.all(readout >= 0))

    def test_wrong_shape(self):
        # Attempting to integrate on a flux array of the wrong shape
        # should raise an exception.
        flux = [1.0, 3.0, 5.0]
        self.assertRaises(TypeError, self.integrator.integrate, flux, 1.0)

    def test_negative_time(self):
        # Attempting to integrate with a negative integration time
        # should raise an exception.
        flux = [[1.0, 2.0, 3.0], \
                [4.0, 5.0, 6.0], \
                [7.0, 8.0, 9.0]]
        self.assertRaises(ValueError, self.integrator.integrate, flux, -1.0)
        # Attempting to wait for a negative time should also raise an exception
        self.assertRaises(ValueError, self.integrator.wait, -1.0)

    def test_cosmic_ray_int(self):
        # Test the effect of cosmic ray hits on the integrator.
        # This test defines a single pixel hit and an integer energy value.
        energy = 1000
        self.integrator.reset()
        # A cosmic ray hit should increase the reading if it hits inside
        # the boundaries of the integrator.
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.hit_by_cosmic_ray(energy, 1, 1)
        readout2 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
        self.assertTrue(diff1.max() > 0)
        # A cosmic ray hit outside the boundary is ignored. (The reading
        # might change due to Poisson statistics, but any change should be
        # a lot smaller than the cosmic ray energy).
        self.integrator.hit_by_cosmic_ray(energy, 3, 3)
        readout3 = self.integrator.readout().astype(np.int32)
        diff2 = readout3 - readout2
        self.assertTrue(diff2.max() < energy/2.0)
        # Negative row and column numbers should not cause a crash.
        self.integrator.hit_by_cosmic_ray(energy, -1, -1)
        
    def test_cosmic_ray_map(self):
        # Test the effect of cosmic ray hits defined as an energy map.
        # First a map smaller than the integrator.
        energy = [[200,300], \
                  [400,100]]
        self.integrator.reset()
        # A cosmic ray hit should increase the reading if any non-zero
        # part of the map hits inside the boundaries of the integrator.
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.hit_by_cosmic_ray(energy, 1, 1)
        readout2 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
        self.assertTrue(diff1.max() > 0)
        # A cosmic ray hit outside the boundary is ignored. (The reading
        # might change due to Poisson statistics, but any change should be
        # a lot smaller than any part of the cosmic ray energy).
        self.integrator.hit_by_cosmic_ray(energy, 4, 4)
        readout3 = self.integrator.readout().astype(np.int32)
        diff2 = readout3 - readout2
        self.assertTrue(diff2.max() < 200)
        # Negative row and column numbers should not cause a crash.
        self.integrator.hit_by_cosmic_ray(energy, -1, -1)
        
        # Test a map the same size as the integrator.
        energy = [[10,100,10], \
                  [100,1000,100], \
                  [10,100,10]]
        self.integrator.reset()
        # A cosmic ray hit should increase the reading if any non-zero
        # part of the map hits inside the boundaries of the integrator.
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.hit_by_cosmic_ray(energy, 1, 1)
        readout2 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
        self.assertTrue(diff1.max() > 0)
        # A cosmic ray hit outside the boundary is ignored. (The reading
        # might change due to Poisson statistics, but any change should be
        # a lot smaller than any part of the cosmic ray energy).
        self.integrator.hit_by_cosmic_ray(energy, 5, 5)
        readout3 = self.integrator.readout().astype(np.int32)
        diff2 = readout3 - readout2
        self.assertTrue(diff2.max() < 500)
        
        # Test a map larger than the integrator.
        energy = [[10,100,100,10], \
                  [100,1000,1000,100], \
                  [10,100,100,10]]
        self.integrator.reset()
        # A cosmic ray hit should increase the reading if any non-zero
        # part of the map hits inside the boundaries of the integrator.
        readout1 = self.integrator.readout().astype(np.int32)
        self.integrator.hit_by_cosmic_ray(energy, 1, 1)
        readout2 = self.integrator.readout().astype(np.int32)
        diff1 = readout2 - readout1
        self.assertTrue(diff1.max() > 0)
        # A cosmic ray hit outside the boundary is ignored. (The reading
        # might change due to Poisson statistics, but any change should be
        # a lot smaller than any part of the cosmic ray energy).
        self.integrator.hit_by_cosmic_ray(energy, 10, 10)
        readout3 = self.integrator.readout().astype(np.int32)
        diff2 = readout3 - readout2
        self.assertTrue(diff2.max() < 500)
        
    def test_bad_cosmic_ray_map(self):
        # Attempting to use a cosmic ray map that isn't 2-D
        # should raise an exception.
        energy = [200,300]
        self.integrator.reset()
        self.assertRaises(TypeError, self.integrator.hit_by_cosmic_ray,
                          energy, 0, 0)

    def test_leakage(self):
        # Test the effect of charge leakage from the integrator.
        # First integrate some counts.
        flux = [[1.0, 2.0, 3.0], \
                [4.0, 5.0, 6.0], \
                [7.0, 8.0, 9.0]]
        self.integrator.reset()
        self.integrator.integrate(flux, 1000.0)
        readout1 = self.integrator.readout()
        # Now leak some charge from one of the pixels.
        # The reading at that pixel should have decreased.
        energy = 1000
        self.integrator.leak(energy, 1, 1)
        readout2 = self.integrator.readout()
        self.assertTrue(readout2[1,1] < readout1[1,1])
        
    def test_zeropoint(self):
        # Zeropoint coefficients, must be a tuple, list or None
        self.assertRaises(AssertionError, self.integrator.set_zeropoint,
                          'silly value', 'silly value')
        self.integrator.set_zeropoint(None, None)
        zpslow = [50000.0, 0.0084]
        zpfast = [[0.0, -2.917],
                  [0.0, -2.292],
                  [0.0, -2.396],
                  [0.0, -2.408]]
        self.integrator.set_zeropoint(zpslow, zpfast)
        # Insert more zeropoint tests here

    def test_latency(self):
        # Latency parameters must either be None or tuples of at least 2
        # values
        self.assertRaises(AssertionError, self.integrator.set_latency,
                          'silly value', 'silly value')
        self.assertRaises(AssertionError, self.integrator.set_latency,
                          42.0, 42.0)
        self.assertRaises(AssertionError, self.integrator.set_latency,
                          [42.0], [42.0])
        self.integrator.set_latency(None, None)
        slow_params = [1.0e-8, 0.999]
        fast_params = [0.002, 0.5]
        self.integrator.set_latency(slow_params, fast_params)
        # Insert more latency tests here.

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
