#!/usr/bin/env python

"""

Module test_amplifier - Contains the unit tests for the
Amplifier class.

:History:

26 Aug 2010: Created
01 Sep 2010: Changed cosmic ray simulation algorithm. The bias level
             is no longer affected.
27 Sep 2010: Python environment for windows verified.
16 Nov 2010: Added gain type test.
16 Nov 2010: Set the seed of the random number generator, for more
             predictable and repeatable tests.
13 Nov 2012: Major restructuring of the package folder. Import
             statements updated.
04 Jun 2013: Removed use of noise_mv object.
08 Sep 2015: Make compatible with Python 3
30 Mar 2015: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
05 May 2017: Changed permission for nosetests.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.FATAL) # Turn off most log messages 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import unittest
import numpy as np

from miri.simulators.scasim.amplifier import Amplifier, set_reference_level, \
    random_level_change, check_readout

_NUMBER_OF_AMPS = 4

class TestAmplifier(unittest.TestCase):
    
    def setUp(self):
        # For test purposes assume there is a 4x4 pixel detector and
        # each of the 4 columns on the detector is read out by a
        # separate amplifier. These 4 amplifier objects are created
        # below, with their bias, gain and noise set to arbitrary values.
        bias = 42
        gain = 1.0
        noise = 20.0
        maxdn = 65535.0
        self.amp_list = []
        for ii in range(0,_NUMBER_OF_AMPS):
            amp = Amplifier(ii, None,None,None, ii,None,_NUMBER_OF_AMPS,
                              bias, gain, noise, maxdn,
                              gaintype='LINEAR', verbose=0, logger=LOGGER)
            amp.set_seed(ii)
            self.amp_list.append(amp)
         
    def tearDown(self):
        # Tidy up
        for amp in self.amp_list:
            del amp
        del self.amp_list
        
    def test_creation(self):
        # Check for pathological cases when creating bad Amplifier objects.
        # The bias and noise values must be valid numbers.
        self.assertRaises(ValueError, Amplifier, 0, None,None,None,
                          None,None,None, 'forty_two', 1.0, 20.0, 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        self.assertRaises(ValueError, Amplifier, 0, None,None,None,
                          None,None,None, 42, 1.0, 'twenty', 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        # The maxdn parameter must be either a float or None.
        self.assertRaises(ValueError, Amplifier, 0, None,None,None,
                          None,None,None, 42, 1.0, 20.0, 'not a float',
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        # Should be ok
        amp = Amplifier(0, None,None,None,
                        None,None,None, 42, 1.0, 20.0, None,
                        gaintype='LINEAR', verbose=0, logger=LOGGER)
        del amp
        # The gain type should be recognised
        self.assertRaises(ValueError, Amplifier, 0, None,None,None,
                          None,None,None, 42, 1.0, 20.0, None,
                          gaintype='NoSuchType', verbose=0, logger=LOGGER)
        # When gaintype is LINEAR the gain must be a single valid number.
        self.assertRaises(TypeError, Amplifier, 0, None,None,None,
                          None,None,None, 42, 'silly gain', 20.0, 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        self.assertRaises(TypeError, Amplifier, 0, None,None,None,
                          None,None,None, 42, None, 20.0, 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        self.assertRaises(TypeError, Amplifier, 0, None,None,None,
                          None,None,None, 42, (0.5, 0.1), 20.0, 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        # When gaintype is POLYNOMIAL the gain must be a tuple or list
        # of numbers, not a single value.
        self.assertRaises(TypeError, Amplifier, 0, None,None,None,
                          None,None,None, 42, 1.0, 20.0, 65535.0,
                          gaintype='POLYNOMIAL', verbose=0, logger=LOGGER)
        # Should be ok
        amp = Amplifier(0, None,None,None,
                          None,None,None, 42, 1, 20.0, 65535.0,
                          gaintype='LINEAR', verbose=0, logger=LOGGER)
        del amp
        amp = Amplifier(0, None,None,None,
                          None,None,None, 42, (1,1), 20.0, 65535.0,
                          gaintype='POLYNOMIAL', verbose=0, logger=LOGGER)
        del amp
        amp = Amplifier(0, None,None,None,
                          None,None,None, 42, [1,1], 20.0, 65535.0,
                          gaintype='POLYNOMIAL', verbose=0, logger=LOGGER)
        del amp

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        for amp in self.amp_list:        
            descr = amp.__str__()
            self.assertIsNotNone(descr)
        
    def test_default_zone(self):
        # It should be possible to create an amplifier with all the row
        # and column attributes set to None and use it to read out a whole
        # detector.
        bias = 42
        gain = (1.0, 0.0) # Default polynomial gain
        noise = 20.0
        maxdn = 65535.0
        amp = Amplifier(0, None,None,None, None,None,None,
                        bias, gain, noise, maxdn, verbose=0)
        detector_data = [[1.0, 5.0, 1.0],
                         [5.0, 9.0, 5.0],
                         [1.0, 5.0, 1.0],
                         ]
        readout = amp.readout(detector_data)
#        print "readout=", readout
        # The entire detector surface should be have been read out
        self.assertTrue(check_readout())
        del amp, detector_data, readout

    def test_readout_mode(self):
        # Verify that the set_readout_mode() method doesn't fail.
        for amp in self.amp_list:        
            amp.set_readout_mode(10)
            
    def test_level_changes(self):
        # Verify that the reference level change functions don't fail.
        set_reference_level(42)
        random_level_change(2)

    def test_readout_check(self):
        # Set up a 4x4 set of detector data and deliberately
        # only read it with one amplifier.
        detector_data = [[1.0, 3.0, 5.0, 7.0],
                         [3.0, 5.0, 7.0, 1.0],
                         [5.0, 7.0, 1.0, 3.0],
                         [7.0, 1.0, 3.0, 5.0],
                         ]
        # The check_readout() function must return False
        # when the detector is only partly read out.
        data = self.amp_list[0].readout(detector_data)
        self.assertIsNotNone(data)
        self.assertFalse(check_readout())

    def test_readout(self):
        # Set up a 4x4 set of detector data and attempt to read it out
        # using the 4 amplifiers.
        detector_data = [[1.0, 3.0, 5.0, 7.0],
                         [3.0, 5.0, 7.0, 1.0],
                         [5.0, 7.0, 1.0, 3.0],
                         [7.0, 1.0, 3.0, 5.0],
                         ]
        for amp in self.amp_list:        
            detector_data = amp.readout(detector_data)
#            print "readout=", detector_data
        # After all the amplifiers have been processed, the entire
        # detector surface should have been read out.
        self.assertTrue(check_readout())

    def test_gain(self):
        # The readout tests above will have generated output looking
        # different from the original detector data because of the
        # bias level and readout noise. The amplifier gain can be
        # tested by setting the bias level readout noise to zero.
        test_gain = 42.0
        original_data = np.array([[1.0, 3.0, 5.0, 7.0],
                         [3.0, 5.0, 7.0, 1.0],
                         [5.0, 7.0, 1.0, 3.0],
                         [7.0, 1.0, 3.0, 5.0],
                         ])
        detector_data = np.copy(original_data)
        set_reference_level(0)
        for amp in self.amp_list:        
            amp._set_bias(0)
            amp._set_read_noise(0.0)
            amp._set_gain(test_gain, gaintype='LINEAR')
        for amp in self.amp_list:        
            detector_data = amp.readout(detector_data)
        # The data read out should have all been multiplied by
        # the linear gain set above.
        residual = abs((detector_data / original_data) - test_gain)
        self.assertTrue(np.all(residual < 0.0001))
        
        # Now test a polynomial gain.
        del detector_data
        detector_data = np.copy(original_data)
        test_gain = (-1.0e-6, 10.0, 3.0)
        for amp in self.amp_list:        
            amp._set_gain(test_gain, gaintype='POLYNOMIAL')
        for amp in self.amp_list:        
            detector_data = amp.readout(detector_data)
        # Verify that each element of the readout has had the defined
        # polynomial applied to it. (Values are only tested to 5
        # significant figures because of possible rounding errors.)
        for row in range(0,4):
            for col in range(0,4):
                expected = test_gain[-1] + \
                    (test_gain[-2] * original_data[row,col]) + \
                    (test_gain[-3] * \
                        original_data[row,col] * original_data[row,col])
                self.assertAlmostEqual(detector_data[row,col], expected,5)

    def test_cosmic_ray_hit(self):
        # Verify each amplifier can be hit by cosmic rays of varying
        # energy. The internal noise level should increase after each
        # hit. There may or may not be an increase in glitches,
        # depending on the cosmic ray energy.
        for energy in (100.0, 10000.0, 1.0e6, 1.0e8):
            for amp in self.amp_list:
                glitches_before = amp._cosmic_ray_glitches
                noise_before = amp._cosmic_ray_noise
                amp.hit_by_cosmic_ray(energy)
                glitches_after = amp._cosmic_ray_glitches
                noise_after = amp._cosmic_ray_noise
                self.assertTrue(glitches_after >= glitches_before)
                self.assertTrue(noise_after > noise_before)
                # A reset should remove all cosmic ray effects.
                amp.reset()
                self.assertEqual(amp._cosmic_ray_glitches, 0)
                self.assertAlmostEqual(amp._cosmic_ray_noise, 0.0)


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
