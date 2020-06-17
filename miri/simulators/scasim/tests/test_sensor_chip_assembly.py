#!/usr/bin/env python

"""

Module test_sensor_chip_assembly - Contains the unit tests for the
SensorChipAssembly class.
NOTE: Since this module depends on all the other scasim modules
it is recommended that all the other modules be tested first, with
the primitive, stand-alone modules tested first.
NOTE: Some tests within this module can take several minutes to run.

:History:

31 Aug 2010: Created
01 Sep 2010: Corrected host-specific time formatting problem.
03 Sep 2010: illumination_maps renamed to data_maps
27 Sep 2010: Python environment for windows verified.
16 Nov 2010: Set the seed of the random number generator, for more
             predictable and repeatable tests.
11 Jan 2011: Added saturation test.
25 Oct 2011: Read cosmic_ray_properties, detector_properties and
             amplifier_properties using ParameterFileManager.
31 Oct 2011: find_writable function used to find a writable
             directory for the temporary files.
11 Jan 2012: Renamed symbol to avoid potential name clash:
             format-->fileformat.
13 Nov 2012: Major restructuring of the package folder. Import
             statements updated.
05 Jun 2013: Tidied up constructor and moved bad pixel mask and dark map
             access to separate functions. Rationalised attribute names.
10 Jun 2013: Modified to use the new MIRI exposure model as an alternative.
             Replaced the old MiriMetadata model with a simpler one in
             data_maps. Metadata keywords and FITS header access functions
             changed to match the new data model.
02 Sep 2013: Added option to suppress output.
09 Sep 2013: Added TEST_MODE_COVERAGE option to reduce the execution time
             for automated testing.
03 Mar 2015: Simulator simplified by removing the ability to save data
             to legacy file formats.
08 Sep 2015: Removed dependency on data_maps module. Made compatible with
             Python 3.
23 Sep 2015: Added MINSAMPLES to ensure that at least 2 possibilities per
             case are tested. simulate_sca global function replaced
             by simulate_files method.
23 Mar 2016: Python logger defined and log messages turned off. All
             combinations of readout modes are tested, but not all have
             corresponding CDP files available.
08 Apr 2016: amplifier_properties removed.
06 Jun 2016: Added test_flags, to ensure the simulator still works when
             flags are turned on and off.
03 Aug 2016: Replaced _set_illumination_map() with _set_illumination and
             renamed subarray parameter to subarray_input.
05 Aug 2016: Shuffle unsampled lists to improve the test coverage.
05 May 2017: Changed permission for nosetests.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
19 Jul 2017: Make the DARK simulation optional for situations where the
             test system has limited memory
17 Jun 2020: Work-around to allow the test to work with nosetests after
             installation by pip. Unzip the data file if not found.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


# Python logging facility
import logging
logging.basicConfig(level=logging.FATAL) # Turn off most log messages 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import os, time
import unittest
import warnings
import random
import numpy as np

from zipfile import ZipFile

# Get the MIRI illumination data model and SensorChipAssembly classes
from miri.datamodels.sim import MiriExposureModel, \
    MiriIlluminationModel
from miri.simulators.scasim.sensor_chip_assembly import SensorChipAssembly,\
    SensorChipAssembly1, SensorChipAssembly2, SensorChipAssembly3

# Search for the cosmic ray and detector parameters files and parse them
# into properties dictionaries. The files are searched for in 3 places:
# (a) The current directory
# (b) The directory where this Python file is being executed from
# (c) The miri.simulators.scasim installation directory.
from miri.tools.filesearching import find_writable, ParameterFileManager, \
    make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
cosmic_ray_properties = ParameterFileManager("cosmic_ray_properties.py",
                            search_path=search_path,
                            description="cosmic ray properties",
                            logger=LOGGER)
detector_properties = ParameterFileManager("detector_properties.py",
                            search_path=search_path,
                            description="detector properties",
                            logger=LOGGER)

# Get the path to where the data directory has been installed. This allows
# the data files to be loaded using an absolute file path
from miri.simulators.scasim import data
_datapath = data.__path__[0]

# Query the lists of known detectors, readout modes and cosmic ray modes
# from the properties dictionaries.
_KNOWN_DETECTORS = list(detector_properties['DETECTORS_DICT'].keys())
_DEFAULT_SCA = _KNOWN_DETECTORS[0]
_DEFAULT_READOUT = 'SLOW'
_DEFAULT_CR_MODE = cosmic_ray_properties['DEFAULT_CR_MODE']

# These test input files are part of the scasim installation.
_TEST_INPUT_FILE = os.path.join(_datapath, 'SCATestInput80x64.fits')
_TEST_INPUT_FILE2 = os.path.join(_datapath, 'SCATestHorseHead1024.fits')
_TEST_INPUT_FILE2_ZIP = os.path.join(_datapath, 'SCATestHorseHead1024.zip')

# Find a writable directory and make up an output file stub guaranteed
# not to overwrite any existing files accidentally.
_filestring = 'SCATestOutput%s' % time.strftime('%Y%m%d%H%M%S')
_outpath = find_writable()
_TEST_OUTPUT_STUB = os.path.join(_outpath, _filestring)
REMOVE_FILES = True

# Control verbose output
VERBOSE = True

# Control the test coverage for readout modes, cosmic ray modes and subarray
# modes. Most problems will manifest themselves in any of these modes, so
# most of the time it is only necessary to test a subset. The subset is
# randomized so all the modes will eventually be tested given a sufficient
# number of repeats.
#
# For infrequent, interactive testing during development it is recommended
# that TEST_MODE_COVERAGE be set to 1.0 to ensure every combination of modes
# is tested. The test will take several minutes.
#
# For frequent, automated testing it is recommended that TEST_MODE_COVERAGE
# be reduced to save time. A value of ~0.25 seems the right compromise between
# paranoia and time taken.
#
# MINSAMPLES defines the minimum number of samples. 2 ensures that at least
# two cases are covered per test.
TEST_MODE_COVERAGE = 0.25
MINSAMPLES = 2

# Control whether DARK simulation is included in the test.
# Set to False in situations where the test system has limited memory.
SIMULATE_DARK = False

class TestSensorChipAssembly(unittest.TestCase):
    
    def setUp(self):
        # Begin by randomizing the random number generator.
        random.seed()
        # Nothing else needs to be set up.
         
    def tearDown(self):
        # Nothing needs to be tidied up
        pass
        
    def test_creation(self):
        # Check for pathological cases when creating bad SensorChipAssembly
        # objects.
        sca = SensorChipAssembly(logger=LOGGER)
        # The focal plane module must exist.
        self.assertRaises(KeyError, sca.setup, 'No such SCA',
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=100.0, ngroups=12, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # The readout mode must be recognisable.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode='No such mode', subarray='FULL',
                          inttime=100.0, ngroups=12, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # The subarray mode must be recognisable.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT,
                          subarray='No such area',
                          inttime=100.0, ngroups=12, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # The cosmic ray mode must be recognisable.
        self.assertRaises(KeyError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=100.0, ngroups=12, nints=2,
                          temperature=6.5, cosmic_ray_mode='in a bucket',
                          verbose=0)
        # Integration time must be a valid number when ngroups is not given.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime='a fortnight', ngroups=None, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # Integration time cannot be zero or negative - there must be at
        # least one group.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=0.0, ngroups=None, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=-10.0, ngroups=None, nints=2,
                          temperature=6.5, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # Temperature must be a valid number.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=100.0, ngroups=12, nints=2,
                          temperature='cold', cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        # There must be at least 1 group or integration.
        self.assertRaises(ValueError, sca.setup, _DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray='FULL',
                          inttime=100.0, ngroups=0, nints=0,
                          temperature=100.0, cosmic_ray_mode=_DEFAULT_CR_MODE,
                          verbose=0)
        del sca
       
    def test_simulation(self):
        # Run a simulation on the test input file for all possible focal
        # plane modules and all possible readout modes. The subarray modes
        # cannot all be tested here because the 80x64 test data is too
        # small to accommodate them.
        # For all possible detector modules
        sca = SensorChipAssembly1(logger=LOGGER)
        for detector in list(detector_properties['DETECTORS_DICT'].keys()):
            if VERBOSE:
                print( "\nTesting detector", detector, "with readout modes: ", \
                       end='')
            if TEST_MODE_COVERAGE > 0.99:
                modes_to_test = list(detector_properties['READOUT_MODE'].keys())
                random.shuffle(modes_to_test)
            else:
                nsamples = int(0.5 + \
                               len(list(detector_properties['READOUT_MODE'].keys())) * \
                               TEST_MODE_COVERAGE)
                nsamples = max(nsamples, MINSAMPLES)
                modes_to_test = random.sample( \
                                    list(detector_properties['READOUT_MODE'].keys()),
                                    nsamples)
            with warnings.catch_warnings(): # Suppress FITS header warnings.
                warnings.simplefilter("ignore")
                for rdmode in modes_to_test:
                    if VERBOSE:
                        print( rdmode + ' ', end='')
                    # Use a longer integration time for SLOW mode.
                    if rdmode == 'SLOW':
                        inttime = 10.0
                    else:
                        inttime = 1.0
                    # Case 1: Specify the integration time and let the simulator
                    # calculate the number of groups.
                    test_output_file_name = "%s_%s_%s_TIME.fits" % \
                        (_TEST_OUTPUT_STUB, detector, rdmode)
                    sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                            detector, readout_mode=rdmode, subarray='FULL',
                            inttime=inttime, nints=1, wait_time=10.0,
                            cosmic_ray_mode='NONE',
                            simulate_dark_current=SIMULATE_DARK,
                            overwrite=True,
                            seedvalue=1, verbose=0 )
                    if REMOVE_FILES:
                        os.remove(test_output_file_name)
                    
                    # Case 2: Specify the number of groups and let the simulator
                    # calculate the integration time.
                    test_output_file_name = "%s_%s_%s_GROUPS.fits" % \
                        (_TEST_OUTPUT_STUB, detector, rdmode)
                    sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                            detector, readout_mode=rdmode, subarray='FULL',
                            ngroups=4, nints=2, wait_time=10.0,
                            cosmic_ray_mode='NONE',
                            simulate_dark_current=SIMULATE_DARK,
                            overwrite=True,
                            seedvalue=2, verbose=0 )
                    if REMOVE_FILES:
                        os.remove(test_output_file_name)
            if VERBOSE:
                print( "", end='')
        del sca

    def test_flags(self):
        # Run a simulation on the test input file for a selection of flag
        # combinations.       
        if VERBOSE:
            print( "\nTesting Flags:", end='')
        sca = SensorChipAssembly1(logger=LOGGER)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")

            if VERBOSE:
                print( " QE=", end='')
            for qe_adjust in (True, False):
                if VERBOSE:
                    print( str(qe_adjust)[0], end='')
                test_output_file_name = "%s_QE_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(qe_adjust))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        qe_adjust=qe_adjust, overwrite=True, seedvalue=1,
                        verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " POISSON=", end='')
            for simulate_poisson_noise in (True, False):
                if VERBOSE:
                    print( str(simulate_poisson_noise)[0], end='')
                test_output_file_name = "%s_POISSON_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_poisson_noise))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        simulate_poisson_noise=simulate_poisson_noise,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " READ=", end='')
            for simulate_read_noise in (True, False):
                if VERBOSE:
                    print( str(simulate_read_noise)[0], end='')
                test_output_file_name = "%s_READ_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_read_noise))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        simulate_read_noise=simulate_read_noise,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " REF=", end='')
            for simulate_ref_pixels in (True, False):
                if VERBOSE:
                    print( str(simulate_ref_pixels)[0], end='')
                test_output_file_name = "%s_REF_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_ref_pixels))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_ref_pixels=simulate_ref_pixels,
                        simulate_dark_current=SIMULATE_DARK,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " BAD=", end='')
            for simulate_bad_pixels in (True, False):
                if VERBOSE:
                    print( str(simulate_bad_pixels)[0], end='')
                test_output_file_name = "%s_BAD_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_bad_pixels))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_bad_pixels=simulate_bad_pixels,
                        simulate_dark_current=SIMULATE_DARK,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if SIMULATE_DARK:
                if VERBOSE:
                    print( " DARK=", end='')
                for simulate_dark_current in (True, False):
                    if VERBOSE:
                        print( str(simulate_dark_current)[0], end='')
                    test_output_file_name = "%s_DARK_%s.fits" % \
                        (_TEST_OUTPUT_STUB, str(simulate_dark_current))
                    sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                            _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                            inttime=1.0, nints=1, wait_time=1.0,
                            cosmic_ray_mode='NONE',
                            simulate_dark_current=simulate_dark_current,
                            overwrite=True, seedvalue=1, verbose=0 )
                    if REMOVE_FILES:
                        os.remove(test_output_file_name)

            if VERBOSE:
                print( " FLAT=", end='')
            for simulate_flat_field in (True, False):
                if VERBOSE:
                    print( str(simulate_flat_field)[0], end='')
                test_output_file_name = "%s_FLAT_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_flat_field))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        simulate_flat_field=simulate_flat_field,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " GAIN=", end='')
            for simulate_gain in (True, False):
                if VERBOSE:
                    print( str(simulate_gain)[0], end='')
                test_output_file_name = "%s_GAIN_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_gain))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        simulate_gain=simulate_gain,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)

            if VERBOSE:
                print( " LIN=", end='')
            for simulate_nonlinearity in (True, False):
                if VERBOSE:
                    print( str(simulate_nonlinearity)[0], end='')
                test_output_file_name = "%s_LIN_%s.fits" % \
                    (_TEST_OUTPUT_STUB, str(simulate_nonlinearity))
                sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode='FAST', subarray='FULL',
                        inttime=1.0, nints=1, wait_time=1.0,
                        cosmic_ray_mode='NONE',
                        simulate_dark_current=SIMULATE_DARK,
                        simulate_nonlinearity=simulate_nonlinearity,
                        overwrite=True, seedvalue=1, verbose=0 )
                if REMOVE_FILES:
                    os.remove(test_output_file_name)
        if VERBOSE:
            print( "", end='')            
        del sca
      
    def test_pipeline(self):
        # Make sure the test file has been unzipped.
        if not os.path.isfile(_TEST_INPUT_FILE2):
            if VERBOSE:
                print( "Unzipping", _TEST_INPUT_FILE2_ZIP )
            zf = ZipFile(_TEST_INPUT_FILE2_ZIP)
            zf.extractall(_datapath)

        # Test the alternative pipeline API for the simulator.
        test_map = MiriIlluminationModel(_TEST_INPUT_FILE2)
        test_map.set_instrument_metadata(_DEFAULT_SCA)
        
        sca = SensorChipAssembly3(logger=LOGGER)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")
            exposure_data = sca.simulate_pipe(test_map, scale=1.0,
                    fringemap=None, readout_mode=_DEFAULT_READOUT,
                    subarray=None, nints=1, ngroups=10, temperature=6.5,
                    cosmic_ray_mode='NONE',
                    simulate_dark_current=SIMULATE_DARK, verbose=0)
            # The simulated data must have a valid shape and contain a
            # range of values.
            self.assertTrue(exposure_data.data.shape[0] > 0)
            self.assertTrue(exposure_data.data.shape[1] > 0)
            self.assertNotAlmostEqual(exposure_data.data.min(),
                                      exposure_data.data.max())
        del test_map, exposure_data, sca
        
    def test_saturation(self):
        # Test that the simulator handles saturated data correctly.
        sca = SensorChipAssembly2(logger=LOGGER)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")
            test_output_file_name = "%s_SATURATED.fits" % _TEST_OUTPUT_STUB
            sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                _DEFAULT_SCA, readout_mode='SLOW', subarray='FULL',
                nints=1, ngroups=20, scale=10000.0, cosmic_ray_mode='NONE',
                simulate_dark_current=SIMULATE_DARK,
                overwrite=True, seedvalue=1, verbose=0 )
            # TODO: Check the data?
            if REMOVE_FILES:
                os.remove(test_output_file_name)
        del sca
                                 
    def test_subarray_modes(self):
        # Test the extraction of subarrays from full frame data.
        test_map = MiriIlluminationModel(_TEST_INPUT_FILE2)
        map_data = test_map.get_illumination()
        # The illumination data is assumed to have a valid shape and be
        # filled with a range of values (the data should have a slope).
        self.assertTrue(map_data.shape[0] > 0)
        self.assertTrue(map_data.shape[1] > 0)
        self.assertNotAlmostEqual(map_data.min(), map_data.max())
        if VERBOSE:
            print( "\nTesting FULL to subarray modes: ", end='')
        if TEST_MODE_COVERAGE > 0.99:
            modes_to_test = list(detector_properties['SUBARRAY'].keys())
            random.shuffle(modes_to_test)
        else:
            nsamples = int(0.5 + len(list(detector_properties['SUBARRAY'].keys())) * \
                           TEST_MODE_COVERAGE)
            nsamples = max(nsamples, MINSAMPLES)
            modes_to_test = random.sample( \
                                list(detector_properties['SUBARRAY'].keys()),
                                nsamples)
        sca = SensorChipAssembly3(logger=LOGGER)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")
            for submode in modes_to_test:
                if VERBOSE:
                    print( submode + ' ', end='')
                sca.setup(_DEFAULT_SCA,
                          readout_mode=_DEFAULT_READOUT, subarray=submode,
                          inttime=None, ngroups=2, nints=1, temperature=6.5,
                          cosmic_ray_mode='NONE',
                          simulate_dark_current=SIMULATE_DARK, verbose=0)
                sca.set_illumination(test_map)
                simulated_data = sca.exposure()
                # The simulated data must have a valid shape and contain a
                # range of values. (A more specific test is difficult because
                # read noise is included. The accuracy of subarray extraction
                # should already have been tested in the test_detector unit
                # tests.)
                self.assertTrue(simulated_data.shape[0] > 0)
                self.assertTrue(simulated_data.shape[1] > 0)
                self.assertNotAlmostEqual(simulated_data.min(),
                                          simulated_data.max())
        if VERBOSE:
            print( "", end='')
        del test_map, sca

    def test_subarray_to_subarray(self):
        # Test the extraction of a subarray from input data containing
        # the same subarray.
        if VERBOSE:
            print( "\nTesting subarray to subarray modes: ", end='')
        if TEST_MODE_COVERAGE > 0.99:
            modes_to_test = list(detector_properties['SUBARRAY'].keys())
            random.shuffle(modes_to_test)
        else:
            nsamples = int(0.5 + len(list(detector_properties['SUBARRAY'].keys())) * \
                           TEST_MODE_COVERAGE)
            nsamples = max(nsamples, MINSAMPLES)
            modes_to_test = random.sample( \
                                list(detector_properties['SUBARRAY'].keys()),
                                nsamples)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")
            for submode in modes_to_test:
                print( submode + ' ', end='')
            
                subarray = detector_properties.get('SUBARRAY',submode)
                if subarray is not None:
                    
                    testvalues = np.fromfunction(
                            lambda i,j: 2.0 + i*1.0 + j*1.0,
                            [subarray[2],subarray[3]])
                    wavelength = np.empty_like(testvalues)
                    wav = np.linspace(5.0, 25.0, wavelength.shape[1])
                    wavelength[:,:] = wav
                    
                    test_map = MiriIlluminationModel(intensity=testvalues,
                                                     wavelength=wavelength)

                    map_data = test_map.get_illumination()
                    # The illumination data is assumed to have a valid shape and be
                    # filled with a range of values (the data should have a slope).
                    self.assertTrue(map_data.shape[0] > 0)
                    self.assertTrue(map_data.shape[1] > 0)
                    self.assertNotAlmostEqual(map_data.min(),map_data.max())
            
                    sca = SensorChipAssembly1(logger=LOGGER)
                    sca.setup(_DEFAULT_SCA,
                              readout_mode=_DEFAULT_READOUT, subarray=submode,
                              inttime=None, ngroups=2, nints=1, temperature=6.5,
                              cosmic_ray_mode='NONE',
                              simulate_dark_current=SIMULATE_DARK, verbose=0)
                    sca.set_illumination(test_map, subarray_input=submode)
                    simulated_data = sca.exposure()
                    # The simulated data must have a valid shape and contain a
                    # range of values. (A more specific test is difficult because
                    # read noise is included. The accuracy of subarray extraction
                    # should already have been tested in the test_detector unit
                    # tests.)
                    self.assertTrue(simulated_data.shape[0] > 0)
                    self.assertTrue(simulated_data.shape[1] > 0)
                    self.assertNotAlmostEqual(simulated_data.min(),
                                          simulated_data.max())
                    del sca
                    del test_map
        if VERBOSE:
            print( "", end='' )

    def test_cosmic_ray_modes(self):
        # Run a simulation in all possible cosmic ray modes.
        if VERBOSE:
            print( "\nTesting cosmic ray modes: ", end='')
        if TEST_MODE_COVERAGE > 0.99:
            modes_to_test = list(cosmic_ray_properties['CR_FLUX'].keys())
            random.shuffle(modes_to_test)
        else:
            nsamples = int(0.5 + \
                           len(list(cosmic_ray_properties['CR_FLUX'].keys())) * \
                           TEST_MODE_COVERAGE)
            nsamples = max(nsamples, MINSAMPLES)
            modes_to_test = random.sample( \
                                list(cosmic_ray_properties['CR_FLUX'].keys()),
                                nsamples)
        variants_to_test = cosmic_ray_properties['CR_LIBRARY_FILES']['VARIANTS']
        sca = SensorChipAssembly2(logger=LOGGER)
        with warnings.catch_warnings(): # Suppress FITS header warnings.
            warnings.simplefilter("ignore")
            for crmode in modes_to_test:
                for variant in variants_to_test:
                    crmode_variant = crmode + variant
                    if VERBOSE:
                        print( crmode_variant + ' ', end='')
                    test_output_file_name = "%s_%s_TIME.fits" % \
                        (_TEST_OUTPUT_STUB, crmode)
                    sca.simulate_files(_TEST_INPUT_FILE, test_output_file_name,
                        _DEFAULT_SCA, readout_mode=_DEFAULT_READOUT,
                        subarray='FULL', inttime=20.0, nints=1,
                        cosmic_ray_mode=crmode_variant,
                        simulate_dark_current=SIMULATE_DARK,
                        overwrite=True, seedvalue=42, verbose=0 )
                    if REMOVE_FILES:
                        os.remove(test_output_file_name)
        if VERBOSE:
            print( "", end='' )
        del sca
            

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
