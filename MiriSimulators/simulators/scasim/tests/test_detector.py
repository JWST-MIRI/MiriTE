#!/usr/bin/env python

"""

Module test_detector - Contains the unit tests for the
DetectorArray class.
NOTE: Since this module depends on the amplifier, poisson_integrator
and measured_variable modules, the unit tests belonging to those modules
should be executed first.

:History:

26 Aug 2010: Created
27 Sep 2010: Python environment for windows verified.
12 Oct 2010: No such FPM should raise ValueError not KeyError.
16 Nov 2010: Set the seed of the random number generator, for more
             predictable and repeatable tests.
26 Nov 2010: Get known detectors directly from the DETECTORS_DICT.
05 Oct 2011: Test the get_counts function.
25 Oct 2011: Read detector_properties using ParameterFileManager.
13 Nov 2012: Major restructuring of the teams/miri folder. Import
             statements updated.
05 Jun 2013: Tidied up constructor and moved bad pixel mask and dark map
             access to separate functions. Rationalised attribute names.
23 Apr 2014: Removed redundant methods.
21 May 2015: Amplifier level removed and gain and read noise defined by
             calibration data products.
08 Sep 2015: Make compatible with Python 3
30 Mar 2015: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
05 May 2017: Changed permission for nosetests.
13 Oct 2017: test_frame_time and test_exposure_time updated to assume the new
             frame time calculation from Mike Ressler. SLOW mode now uses
             8 out of 9 samples. Test the new frame_rti function.
             Import more detector parameters from detector_properties.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.FATAL) # Turn off most log messages 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import os
import unittest
import numpy as np

from miri.simulators.scasim.cosmic_ray import CosmicRay
from miri.simulators.scasim.detector import frame_rti, DetectorArray

from miri.tools.filesearching import ParameterFileManager, make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
detector_properties = ParameterFileManager("detector_properties.py",
                            search_path=search_path,
                            description="detector properties",
                            logger=LOGGER)
#from miri.simulators.scasim import detector_properties

# Actual detector parameters
_KNOWN_DETECTORS = list(detector_properties['DETECTORS_DICT'].keys())
_FIRST_DETECTOR = detector_properties['DETECTORS_DICT'][str(_KNOWN_DETECTORS[0])]

# Smaller-sized detector used for testing.
_PIXELS_PER_SIDE = 16
_REF_PIXELS_LEFT = 4
_REF_PIXELS_RIGHT = 4
_REF_PIXELS_BOTTOM = 0
_REF_PIXELS_TOP = 4

class TestDetectorArray(unittest.TestCase):
    
    def setUp(self):
        # For test purposes assume there is a 16x16 pixel detector with
        # 4 extra reference columns at the left and right edge, no
        # reference rows at the bottom but an extra 4 reference rows at
        # the top.
        self.detector = DetectorArray(_KNOWN_DETECTORS[0],
                                      _PIXELS_PER_SIDE, _PIXELS_PER_SIDE,
                                      _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                                      left_columns=_REF_PIXELS_LEFT,
                                      right_columns=_REF_PIXELS_RIGHT,
                                      bottom_rows=_REF_PIXELS_BOTTOM,
                                      top_rows=_REF_PIXELS_TOP,
                                      well_depth=_FIRST_DETECTOR['WELL_DEPTH'],
                                      verbose=0, logger=LOGGER)
        self.detector.set_seed(42)
        samples = detector_properties.get('READOUT_MODE', 'SLOW')
        self.detector.set_readout_mode(samples[0], samples[1])
         
    def tearDown(self):
        # Tidy up
        del self.detector
        
    def test_creation(self):
        # Check for pathological cases when creating bad DetectorArray objects.
        # All focal plane modules are supported and an unknown module id
        # will raise an exception.
        for detid in _KNOWN_DETECTORS:
            det = DetectorArray(detid, _PIXELS_PER_SIDE, _PIXELS_PER_SIDE,
                                _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                                left_columns=_REF_PIXELS_LEFT,
                                right_columns=_REF_PIXELS_RIGHT,
                                bottom_rows=_REF_PIXELS_BOTTOM,
                                top_rows=_REF_PIXELS_TOP,
                                well_depth=_FIRST_DETECTOR['WELL_DEPTH'], verbose=0,
                                logger=LOGGER)
            del det
        self.assertRaises(KeyError, DetectorArray, 'No such SCA',
                          _PIXELS_PER_SIDE, _PIXELS_PER_SIDE,
                          _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                          left_columns=_REF_PIXELS_LEFT,
                          right_columns=_REF_PIXELS_RIGHT,
                          bottom_rows=_REF_PIXELS_BOTTOM,
                          top_rows=_REF_PIXELS_TOP,
                          well_depth=_FIRST_DETECTOR['WELL_DEPTH'], verbose=0, logger=LOGGER)
        # Zero or negative detector sizes should be rejected.
        self.assertRaises(ValueError, DetectorArray, _KNOWN_DETECTORS[0],
                          0, 0, _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                          left_columns=_REF_PIXELS_LEFT,
                          right_columns=_REF_PIXELS_RIGHT,
                          bottom_rows=_REF_PIXELS_BOTTOM,
                          top_rows=_REF_PIXELS_TOP,
                          well_depth=_FIRST_DETECTOR['WELL_DEPTH'], verbose=0, logger=LOGGER)
        self.assertRaises(ValueError, DetectorArray, _KNOWN_DETECTORS[0],
                          -1, -1, _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                          left_columns=_REF_PIXELS_LEFT,
                          right_columns=_REF_PIXELS_RIGHT,
                          bottom_rows=_REF_PIXELS_BOTTOM,
                          top_rows=_REF_PIXELS_TOP,
                          well_depth=_FIRST_DETECTOR['WELL_DEPTH'], verbose=0)
        self.assertRaises(ValueError, DetectorArray, _KNOWN_DETECTORS[0],
                          _PIXELS_PER_SIDE, _PIXELS_PER_SIDE, 6.5,
                          left_columns=-1, right_columns=-1,
                          bottom_rows=-1, top_rows=-1,
                          well_depth=_FIRST_DETECTOR['WELL_DEPTH'], verbose=0, logger=LOGGER)
        # Having no reference pixels at all should be ok
        det = DetectorArray(_KNOWN_DETECTORS[0],
                            _PIXELS_PER_SIDE, _PIXELS_PER_SIDE,
                            _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                            left_columns=0, right_columns=0,
                            bottom_rows=0, top_rows=0,
                            well_depth=_FIRST_DETECTOR['WELL_DEPTH'],
                            verbose=0, logger=LOGGER)
        del det

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        descr = self.detector.__str__()
        self.assertIsNotNone(descr)
        shape = self.detector.illuminated_shape
        self.assertEqual(shape[0], _PIXELS_PER_SIDE)
        self.assertEqual(shape[1], _PIXELS_PER_SIDE)
        shape = self.detector.get_subarray_shape()
        shape = self.detector.get_subarray_shape(subarray=(1,1,2,2))
        counts = self.detector.get_counts()
        del counts

    def test_readout_mode(self):
        # Test the setting of the readout mode
        self.detector.set_readout_mode(1, 0)
        self.detector.set_readout_mode(8, 1)
        
    def test_frame_time(self):
        # First test the global frame_rti function by recalculating some of
        # the key rows in Mike Ressler's table.
        test_parameters = [(  1, 258,   1, 1024, 0, 1, 4, False,  2.77504),
                           (  1, 258,   1, 1024, 0, 1, 6, False,  2.79552),
                           (  1,  64,  65,  320, 0, 1, 4, False,  0.24320),
                           (  2,   5,   1,   16, 0, 1, 4, False,  0.07296),
                           (  2,   5,   1,   16, 0, 1, 4, True,   0.07296),
                           (  2,  65,  65,  320, 0, 1, 4, False,  0.24576),
                           ( 91, 106, 765,  828, 0, 1, 4, False,  0.14144),
                           ( 91, 106, 765,  828, 0, 1, 4, True,   0.09600),
                           ( 96,  99, 787,  802, 0, 1, 4, False,  0.08800),
                           ( 96,  99, 787,  802, 0, 1, 4, True,   0.07600),   
                           (181, 244,  65,  320, 0, 1, 4, False,  0.70400),
                           (181, 244,  65,  320, 0, 1, 4, True,   0.33792),
                           
                           (  1, 258,   1, 1024, 1, 8, 4, False, 23.88992),
                           (  1, 258,   1, 1024, 2, 8, 4, False, 26.51136)
                           ]
        REFPIX_SAMPLESKIP = 3
        for (colstart, colstop, rowstart, rowstop, sampleskip, samplesum,
             resetwidth, burst_mode, expected) in test_parameters:
            frame_clks = frame_rti(1024, 258,
                                   colstart, colstop, rowstart, rowstop,
                                   sampleskip, REFPIX_SAMPLESKIP, samplesum,
                                   resetwidth, _FIRST_DETECTOR['RESET_OVERHEAD'],
                                   burst_mode)
            frame_time = frame_clks * _FIRST_DETECTOR['CLOCK_TIME']
            msg = "Frame time calculation for \'%s\' incorrect " % \
                str( (colstart, colstop, rowstart, rowstop, sampleskip,
                      samplesum, resetwidth, burst_mode) )
            msg += "(%g != %g)." % \
                (frame_time, expected)
            self.assertAlmostEqual(frame_time, expected, places=5, msg=msg)

        # Repeat this test for all known MIRI detectors
        for detid in _KNOWN_DETECTORS:
            # The frame time for a MIRI 1024x1024 pixel detector in
            # full frame mode should be around 2.77504 seconds in FAST
            # mode (samplesum=1, sampleskip=0) and around 23.88992 seconds
            # in SLOW mode (samplesum=8, sampleskip=1).
            miri_detector = DetectorArray(detid,
                                      _FIRST_DETECTOR['ILLUMINATED_ROWS'],
                                      _FIRST_DETECTOR['ILLUMINATED_COLUMNS'],
                                      _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                                      left_columns=_FIRST_DETECTOR['LEFT_COLUMNS'],
                                      right_columns=_FIRST_DETECTOR['RIGHT_COLUMNS'],
                                      bottom_rows=_FIRST_DETECTOR['BOTTOM_ROWS'],
                                      top_rows=_FIRST_DETECTOR['TOP_ROWS'],
                                      well_depth=_FIRST_DETECTOR['WELL_DEPTH'],
                                      verbose=0, logger=LOGGER)
            full_fast_time = miri_detector.frame_time(1, 0)
            self.assertAlmostEqual(full_fast_time, 2.77504, places=5)
            full_slow_time = miri_detector.frame_time(8, 1)
            self.assertAlmostEqual(full_slow_time, 23.88992, places=5)
            
            # The subarray readout time should be less than the full
            # frame readout.
            sub_fast_time = miri_detector.frame_time(1, 0, subarray=(1,1,256,256))
            self.assertTrue(sub_fast_time < full_fast_time)
            sub_slow_time = miri_detector.frame_time(8, 1, subarray=(1,1,256,256))
            self.assertTrue(sub_slow_time < full_slow_time)        
            del miri_detector
        
    def test_exposure_time(self):
        # The exposure time for a particular readout mode should
        # scale with the number of integrations and number of groups.
        time_1_1, elapsed_1_1 = self.detector.exposure_time(1, 1, 8, 1)
        time_1_10, elapsed_1_10 = self.detector.exposure_time(1, 10, 8, 1)
        time_10_1, elapsed_10_1 = self.detector.exposure_time(10, 1, 8, 1)
        time_10_10, elapsed_10_10 = self.detector.exposure_time(10, 10, 8, 1)
        ratio1 = time_1_10 / time_1_1
        ratio2 = time_10_1 / time_1_1
        ratio3 = time_10_10 / time_1_1
        self.assertAlmostEqual(ratio1, 10.0)
        self.assertAlmostEqual(ratio2, 10.0)
        self.assertAlmostEqual(ratio3, 100.0)
        # The elapsed time should always be greater than or equal to
        # the integration time.
        self.assertTrue(elapsed_1_1 >= time_1_1)
        self.assertTrue(elapsed_1_10 >= time_1_10)
        self.assertTrue(elapsed_10_1 >= time_10_1)
        self.assertTrue(elapsed_10_10 >= time_10_10)

    def test_integrate(self):
        # Test integrating on a photon flux which is all 1.0
        # with a short integration and a long integration.
        flux = np.ones([_PIXELS_PER_SIDE,_PIXELS_PER_SIDE])
        self.detector.reset()
        self.detector.integrate(flux, 0.1)
        readout1 = self.detector.readout()
        self.detector.integrate(flux, 100.0)
        readout2 = self.detector.readout()
        # The readouts must be positive everywhere.
        self.assertTrue(np.all(readout1 >= 0.0))
        self.assertTrue(np.all(readout2 >= 0.0))
        
        # Test waiting for an elapsed time and integrating again.
        self.detector.wait(42.0)
        self.detector.integrate(flux, 100.0)
        readout3 = self.detector.readout()
        self.assertTrue(np.all(readout3 >= 0.0))
        
    def test_wrong_shape(self):
        # Attempting to integrate on a flux array of the wrong shape
        # should raise an exception.
        flux = [1.0, 3.0, 5.0]
        self.assertRaises(AttributeError, self.detector.integrate, flux, 1.0)

    def test_subarray(self):
        # Create a DetectorArray object with no reference pixels, so the
        # size of the subarray is more predictable.
        detector = DetectorArray(_KNOWN_DETECTORS[0],
                                  _PIXELS_PER_SIDE, _PIXELS_PER_SIDE,
                                  _FIRST_DETECTOR['TARGET_TEMPERATURE'],
                                  left_columns=0, right_columns=0,
                                  bottom_rows=0, top_rows=0,
                                  well_depth=_FIRST_DETECTOR['WELL_DEPTH'],
                                  verbose=0, logger=LOGGER)
        dshape = detector.detector_shape
        # Create a flux array the same shape as the detector whose values
        # are predictable from their row and column.
        # flux(x,y) will contain x + 100y.
        x, y = np.meshgrid(list(range(0,dshape[0])),list(range(0,dshape[1])))
        flux = x + 100.0 * y
        # Try a subarray at the bottom left corner which fits within
        # the detector data. The subarray must have the shape requested
        # and must start at the pixel expected.
        subarray1 = detector._extract_subarray(flux, subarray=(1,1,8,8))
        self.assertEqual(subarray1.shape[0], 8)
        self.assertEqual(subarray1.shape[1], 8)
        self.assertEqual(subarray1[0,0], flux[0,0])
        # Try a subarray lying completely within the detector data.
        # Again, the subarray must have the size and start point expected.
        subarray2 = detector._extract_subarray(flux, subarray=(5,6,6,7))
        self.assertEqual(subarray2.shape[0], 6)
        self.assertEqual(subarray2.shape[1], 7)
        self.assertEqual(subarray2[0,0], flux[4,5])
        del detector
        
        # Repeat the above tests using the default detector which has
        # reference rows.
        dshape = self.detector.detector_shape
        # Create a flux array the same shape as the detector whose values
        # are predictable from their row and column.
        # flux(x,y) will contain x + 100y.
        x, y = np.meshgrid(list(range(0,dshape[0])),list(range(0,dshape[1])))
        flux = x + 100.0 * y
        # Try a subarray at the bottom left corner which fits within
        # the detector data. The subarray must start at the pixel expected
        # and have the number of columns requested (it may have more rows
        # than requested because of the reference rows).
        subarray1 = self.detector._extract_subarray(flux, subarray=(1,1,8,8))
        self.assertTrue(subarray1.shape[0] >= 8)
        self.assertEqual(subarray1.shape[1], 8)
        self.assertEqual(subarray1[0,0], flux[0,0])
        # Try a subarray lying completely within the detector data.
        # Again, the subarray must have the start point and number of columns
        # expected (and the number of rows might be greater).
        subarray2 = self.detector._extract_subarray(flux, subarray=(5,6,6,7))
        self.assertTrue(subarray2.shape[0] >= 6)
        self.assertEqual(subarray2.shape[1], 7)
        self.assertEqual(subarray2[0,0], flux[4,5])

    def test_bad_subarray(self):
        # A readout with a stupid subarray shape should generate an exception
        self.detector.reset()
        self.assertRaises(ValueError, self.detector.readout,
                          subarray=(1,1,0,0))
        # A readout with a subarray completely outside the detector area
        # should generate an exception.
        self.assertRaises(ValueError, self.detector.readout,
                          subarray=(10000,10000,8,8))        
        # A readout with a subarray partially outside the detector area
        # will generate an exception when there are reference rows.
        self.assertRaises(ValueError, self.detector.readout,
                          subarray=(10,10,100,100))        

    def test_cosmic_ray(self):
        # Create a list of 8 cosmic ray objects of various energies
        # and targets. The 7th cosmic ray hits completely outside
        # the detector area - it should be accepted and ignored (this
        # is already tested in the units tests for the PoissonIntegrator
        # class).
        energies = (1000.0,    300.0,       100.0,      5000.0,
                    10000.0,   42.0,        780.0,      65000.0)
        coords = ((4,4),       (6,6),      (2,10),     (1,1),
                  (14,15),     (3,3),      (100,100),  (13,8) )
        hit_map_base = np.array([[0.0, 0.15, 0.0],
                                 [0.15, 1.0, 0.15],
                                 [0.0, 0.15, 0.0]])
        cosmic_ray_list = []
        for ii in range(0,len(energies)):
            hit_map = hit_map_base * energies[ii]
            cosmic_ray = CosmicRay(energies[ii],coords[ii],  hit_map, verbose=0)
            cosmic_ray_list.append(cosmic_ray)
        
        # Hit the detector with the cosmic rays. The readout after
        # the cosmic ray hits will contain spikes that will increase
        # the maximum.
        self.detector.reset()
        readout1 = self.detector.readout()
        self.detector.hit_by_cosmic_rays(cosmic_ray_list)
        readout2 = self.detector.readout()
        self.assertTrue(readout2.max() > readout1.max())


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
