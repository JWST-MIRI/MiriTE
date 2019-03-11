#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""

Simple tests of the data model utility functions in datamodels/util.py.

:History:

19 Jan 2018: created
27 Apr 2018: Specify overwrite=True when calling util.verify_cdp_file
05 Sep 2018: Flatfield CDPs must now include FILTER, CHANNEL and BAND.
             Rewritten to use Python context manager.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import numpy as np

import astropy.io.fits as pyfits

import miri.datamodels.util as util
import miri.datamodels.cdp as cdp
import miri.datamodels.sim as sim

KEEPFILES = False

class TestFileIO(unittest.TestCase):
    
    # Test the file I/O utilities within datamodels.util.
    
    def setUp(self):
        # Create some simple data models which adhere to the MIRI and JWST
        # data model standards and save them to temporary FITS files.
        sci1d = [1.0, 2.0, 3.0, 4.0, 5.0]
        sci2d = [sci1d, sci1d, sci1d, sci1d, sci1d]
        sci3d = [sci2d, sci2d, sci2d]
        sci4d = [sci3d, sci3d]
        self.cdp_models_to_test = \
            [('TestUtilFlatfield.fits', cdp.MiriFlatfieldModel),
             ('TestUtilGain.fits', cdp.MiriGainModel)]
        self.sim_models_to_test = \
            [('TestUtilRamp.fits', sim.MiriRampModel)]
        self.models_to_test = self.cdp_models_to_test + self.sim_models_to_test
        
        with cdp.MiriFlatfieldModel( data=sci2d ) as flatfield:
            flatfield.set_telescope()
            flatfield.set_housekeeping_metadata(origin='UKATC', author='Joe Bloggs',
                                                pedigree='GROUND', version='V1.0',
                                                description='Dummy flat-field data')        
            flatfield.set_instrument_metadata(detector='MIRIMAGE', filt='F560W',
                                              channel='N/A', band='N/A')
            flatfield.set_exposure_metadata(readpatt='SLOW', nints=2, ngroups=3)
            flatfield.set_subarray_metadata((1,1,5,5))
            flatfield.set_exposure_type()
            flatfield.add_history('DOCUMENT: Marvel Comics')
            flatfield.add_history('SOFTWARE: test_util.py')
            flatfield.add_history('DATA USED: Lots')
            flatfield.add_history('DIFFERENCES: Not a lot')
            flatfield.save( self.cdp_models_to_test[0][0], overwrite=True )
        del flatfield

        with cdp.MiriGainModel( data=sci2d ) as gain:
            gain.set_telescope()
            gain.set_housekeeping_metadata(origin='UKATC', author='Joe Bloggs',
                                           pedigree='GROUND', version='V1.0',
                                           description='Dummy gain data')
            gain.set_instrument_metadata(detector='MIRIMAGE', filt='F560W')
            gain.set_exposure_metadata(readpatt='SLOW', nints=2, ngroups=3)
            gain.set_subarray_metadata((1,1,5,5))
            gain.set_exposure_type()
            gain.add_history('DOCUMENT: Marvel Comics')
            gain.add_history('SOFTWARE: test_util.py')
            gain.add_history('DESCRIP: Dummy gain data')
            gain.add_history('DATA USED: Lots')
            gain.add_history('DIFFERENCES: Not a lot')
            gain.save( self.cdp_models_to_test[1][0], overwrite=True )
        del gain

        with sim.MiriRampModel( data=sci4d ) as ramp:
            ramp.set_telescope()
            ramp.set_housekeeping_metadata('UKATC', 'Joe Bloggs', 'GROUND', 'V1.0')
            ramp.set_instrument_metadata(detector='MIRIMAGE', filt='F560W')
            ramp.set_exposure_metadata(readpatt='SLOW', nints=2, ngroups=3)
            ramp.set_subarray_metadata((1,1,5,5))
            ramp.set_exposure_type()
            ramp.save( self.sim_models_to_test[0][0], overwrite=True )
        del ramp
       
    def tearDown(self):
        # Tidy up. Remove the temporary FITS files.
        for (filename, datamodel) in self.models_to_test:
            if os.path.isfile(filename) and not KEEPFILES:
                os.remove(filename)

    def test_open_from_filename(self):
        # Check that each data model can be opened successfully
        # by specifying its filename.
        for (filename, datamodel) in self.models_to_test:
            #print("Opening %s using filename" % filename)
            with util.open( filename ) as model:
                #print(model)
                self.assertIsNotNone(model)
                self.assertIsNotNone(model.data)
                self.assertTrue(isinstance(model, datamodel))
            del model

    def test_open_from_hdulist(self):
        # Check that each data model can be opened successfully
        # by specifying an hdulist.
        for (filename, datamodel) in self.models_to_test:
            #print("Opening %s using hdulist" % filename)
            with pyfits.open( filename ) as hdulist:
                model = util.open( hdulist )
                #print(model)
                self.assertIsNotNone(model)
                self.assertIsNotNone(model.data)
                self.assertTrue(isinstance(model, datamodel))
                hdulist.close()
            del hdulist
            del model

    def test_verify_fits_file(self):
        # Check that the FITS verification function passes the simple
        # files created by this test.
        for (filename, datamodel) in self.models_to_test:
            util.verify_fits_file(filename, cdp_checks=False)

    def test_verify_cdp_file(self):
        # Check that the CDP verification function passes the simple
        # files created by this test.
        for (filename, datamodel) in self.cdp_models_to_test:
            util.verify_cdp_file(filename, overwrite=True)

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
