#!/usr/bin/env python

"""

Module test_cosmic_ray - Contains the unit tests for the
CosmicRayEnvironment and CosmicRay classes.

:History:

26 Aug 2010: Created.
30 Aug 2010: Further tests.
03 Sep 2010: Separated tests on CosmicRay and CosmicRayEnvironment objects.
27 Sep 2010: Python environment for windows verified.
16 Nov 2010: Set the seed of the random number generator, for more
             predictable and repeatable tests.
13 Nov 2012: Major restructuring of the package folder. Import
             statements updated.
08 Sep 2015: Make compatible with Python 3
05 May 2017: Changed permission for nosetests.
28 Sep 2021: Changed np.int to np.int32.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


# Python logging facility
import logging
logging.basicConfig(level=logging.FATAL) # Turn off most log messages 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import unittest
import numpy as np

from miri.simulators.scasim.cosmic_ray import CosmicRayEnvironment, CosmicRay

_NUMBER_OF_MEASUREMENTS = 10

class TestCosmicRay(unittest.TestCase):

    def setUp(self):
        # Create a cosmic ray object.
        im = ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.))
        self.cosmic_ray = CosmicRay(1000.0, (10,10), hit_map=im,
                                    nucleon='C', verbose=0, logger=LOGGER)
                
    def tearDown(self):
        # Tidy up
        del self.cosmic_ray
    
    def test_creation(self):
        # Check for exceptions when creating bad CosmicRay objects.
        im = ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.))
        # If target is DETECTOR coordinates are expected to be a tuple
        # of 2 integers.
        self.assertRaises(TypeError, CosmicRay, 42.0, 3,
                          hit_map=im, nucleon='H', verbose=0, logger=LOGGER)
        self.assertRaises(TypeError, CosmicRay, 42.0, (3,4,5),
                          hit_map=im, nucleon='H', verbose=0, logger=LOGGER)

    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        descr = self.cosmic_ray.__str__()
        self.assertIsNotNone(descr)
        target_coords = self.cosmic_ray.get_target_coords()
        self.assertIsNotNone(target_coords)
        energy = self.cosmic_ray.get_energy()
        self.assertIsNotNone(energy)
        electrons = self.cosmic_ray.get_electrons()
        self.assertIsNotNone(electrons)
        electrons = self.cosmic_ray.get_electrons(scale=10.0)
        self.assertIsNotNone(electrons)
        nucleon = self.cosmic_ray.get_nucleon()
        self.assertIsNotNone(nucleon)
        hit_map = self.cosmic_ray.get_hit_map()
        self.assertIsNotNone(hit_map)
        hit_map = self.cosmic_ray.get_hit_map(scale=10.0)
        self.assertIsNotNone(hit_map)
        

class TestCosmicRayEnvironment(unittest.TestCase):
    
    def setUp(self):
        # Define arrays containing cosmic ray energies and a probability
        # distribution and create a RANDOM mode cosmic ray environment.
        self.energies_rnd = (1.0, 3.0, 10.0, 30.0, 100.0,
                             300.0, 1000.0, 3000.0, 10000.0, 30000.0)
        self.dist_rnd = (0.05, 0.07, 0.13, 0.15, 0.16,
                         0.13, 0.1, 0.09, 0.07, 0.05)
        self.cr_flux = 1.0 # Cosmic ray events per second per square micron
        self.cr_env_rnd = CosmicRayEnvironment(self.cr_flux, 'RANDOM',
                                               energies=self.energies_rnd,
                                               distribution=self.dist_rnd,
                                               images=None, nucleons=None,
                                               verbose=0, logger=LOGGER)
        self.cr_env_rnd.set_seed(42)
        
        # Define a 5 element cosmic ray library and use it to create a
        # LIBRARY mode cosmic ray environment. Each cosmic ray event is
        # described by a 3x3 image.
        self.energies_lib = (1.0, 3.0, 10.0, 30.0, 100.0)
        self.images = (
                  ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
                  ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
                  ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
                  ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
                  ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.))
                  )
        self.nucleons = (0, 1, 2, 3, 4)
        self.cr_env_lib = CosmicRayEnvironment(self.cr_flux, 'LIBRARY',
                                               energies=self.energies_lib,
                                               distribution=None,
                                               images=self.images,
                                               nucleons=self.nucleons,
                                               verbose=0, logger=LOGGER)
        self.cr_env_lib.set_seed(42)
        
    def tearDown(self):
        # Tidy up
        del self.energies_rnd, self.energies_lib, self.dist_rnd
        del self.images, self.nucleons
        del self.cr_env_rnd
        del self.cr_env_lib
        
    def test_creation(self):
        # Check for exceptions when creating bad CosmicRayEnvironment objects.
        # In RANDOM mode, the energies and distribution arrays must be the
        # same size.
        en = [1.0, 10.0, 100.0, 1000.0]
        dist = [0.1, 0.1, 0.1]
        self.assertRaises(TypeError, CosmicRayEnvironment, self.cr_flux,
                          'RANDOM', energies=en, distribution=dist,
                          images=None, nucleons=None, verbose=0, logger=LOGGER)
        # The distribution array must be present in RANDOM mode.
        self.assertRaises(TypeError, CosmicRayEnvironment, self.cr_flux,
                          'RANDOM', energies=en, distribution=None,
                          images=None, nucleons=None, verbose=0, logger=LOGGER)
        # The distribution array cannot be all zero.
        dist = [0.0, 0.0, 0.0, 0.0]
        self.assertRaises(ValueError, CosmicRayEnvironment, self.cr_flux,
                          'RANDOM', energies=en, distribution=dist,
                          images=None, nucleons=None, verbose=0, logger=LOGGER)
        # In LIBRARY mode the energies, images and nucleons arrays must
        # all contain the same number of events.
        en = [1.0, 3.0, 10.0]
        im = (
              ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
              ((0.,1.,0.),(1.,10.,1.),(0.,1.,0.)),
              )
        nu = (0, 1, 2)
        self.assertRaises(TypeError, CosmicRayEnvironment, self.cr_flux,
                          'LIBRARY', energies=en, distribution=None,
                          images=im, nucleons=nu, verbose=0, logger=LOGGER)
        # The images array must be a data cube.
        im = (0.0, 1.0, 0.0)
        self.assertRaises(TypeError, CosmicRayEnvironment, self.cr_flux,
                          'LIBRARY', energies=en, distribution=None,
                          images=im, nucleons=nu, verbose=0, logger=LOGGER)
    
    def test_description(self):
        # Test that the querying and description functions work.
        # For the test to pass these only need to run without error.
        descr1 = self.cr_env_rnd.__str__()
        descr2 = self.cr_env_lib.__str__()
        self.assertIsNotNone(descr1)
        self.assertIsNotNone(descr2)
        
    def test_generate_events_random(self):
        # Define some test cosmic ray environments with a probability
        # distribution designed to favour a particular cosmic ray
        # energy. When random cosmic ray events are generated they should
        # have that same energy.
        energies = [100.0, 200.0, 300.0, 400.0]
        distribution = [0.0, 0.0, 0.0, 0.0]
        
        for ii in range(0,len(energies)):
            distribution[ii] = 1.0
        
            cr_env = CosmicRayEnvironment(self.cr_flux, 'RANDOM',
                                      energies=energies,
                                      distribution=distribution,
                                      images=None, nucleons=None,
                                      verbose=0, logger=LOGGER)
            cr_env.set_seed(ii)
            cosmic_ray1 = cr_env.generate_event(1024, 1024)
            energy1 = cosmic_ray1.get_energy()
            self.assertAlmostEqual(energy1, energies[ii])
            distribution[ii] = 0.0
            del cr_env, cosmic_ray1
        
    def test_generate_events_library(self):
        # Generate a list of random events based on the library
        # provided. There is an even probability of selecting any
        # particular event from the library, so if sufficient
        # events have been generated there should be at least one
        # example of each event in the list.
        cosmic_ray_list = self.cr_env_lib.generate_events(10, 10, 10.0, 1.0)
        nevents = len(cosmic_ray_list)
        npossibilities = len(self.energies_lib)
        if nevents > (npossibilities * 10):
            ecount = np.zeros([len(self.energies_lib)], dtype=np.int32)
            for cosmic_ray in cosmic_ray_list:
                energy = cosmic_ray.get_energy()
                for ii in range(0,npossibilities):
                    if (energy == self.energies_lib[ii]):
                        ecount[ii] += 1
                        break
            # Assert there are counts in every element of the ecount array.
            self.assertTrue(np.all(ecount > 0))
        del cosmic_ray_list
        
    def test_generate_events_quantity(self):
        # Generate a random set of cosmic ray events for detectors only.
        # There should be more events when the pixel size is increased
        # or the integration time is increased.
        cosmic_ray_list1 = self.cr_env_rnd.generate_events(10, 10,
                                                1.0, 1.0)
        cosmic_ray_list2 = self.cr_env_rnd.generate_events(10, 10,
                                                1.0, 10.0)
        cosmic_ray_list3 = self.cr_env_rnd.generate_events(10, 10,
                                                10.0, 1.0)
        num1 = len(cosmic_ray_list1)
        num2 = len(cosmic_ray_list2)
        num3 = len(cosmic_ray_list3)
        self.assertTrue(num2 > num1)
        self.assertTrue(num3 > num1)
        del cosmic_ray_list1, cosmic_ray_list2, cosmic_ray_list3
        

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
