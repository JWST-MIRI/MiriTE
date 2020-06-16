#!/usr/bin/env python

"""

Module test_spectral_models - Contains the unit tests for the data model
classes in the miriteam.spectroscopy module.

:History:

05 Jul 2013: Created.
12 Sep 2013: Swapped the MRS CHANNEL and BAND keywords.
21 Jul 2014: SW detector changed to MIRIFUSHORT.
11 Mar 2015: group_integration_time changed to group_time.
30 Sep 2016: Make all arrays with swapped axes contiguous.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Oct 2019: Moved from MiriTeam to MiriTE to overcome an
             issue with AsdfExtension.

@author: Steven Beard (UKATC)

"""

import os
import unittest
import warnings

import numpy as np
import numpy.ma as ma

from miri.datamodels.spectraldata import Spectrum1D, Spectrum2D, Spectrum3D

class TestSpectrum1D(unittest.TestCase):
    
    # Test the Spectrum1D class.
    # Most of the relevant tests are already done in MiriMeasuredImageModel.
    # Only the additional tests specific to Spectrum1D are included here.

    def setUp(self):
        # Create a typical 1-D spectrum product.
        a = [10,15,10,20,10,30,90,30,20,50,10,12]
        b = [5,4,3,2,1,2,3,4,5,6,7,8]
        c = [0,0,0,0,0,0,1,0,0,0,0,0]
        w = [5,7,9,11,13,15,17,19,21,23,25,27]
        self.dataproduct = Spectrum1D(data=a, err=b, dq=c, wavelength=w,
                                      unit='photons', wunit='microns')

        # Add some example metadata.
        self.dataproduct.set_observation_metadata()
        self.dataproduct.set_target_metadata(0.0, 0.0)
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriSpectrum1D_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        # Remove temporary file, if able to.
        if os.path.isfile(self.testfile):
            try:
                os.remove(self.testfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.testfile + \
                    "\n   " + str(e)
                warnings.warn(strg)

    def test_creation(self):
        # Test the creation of the data product. 
        a = [10,15,10,20,10,30,90,30,20,50,10,12]
        b = [5,4,3,2,1,2,3,4,5,6,7,8]
        c = [0,0,0,0,0,0,1,0,0,0,0,0]
        w = [5,7,9,11,13,15,17,19,21,23,25,27]
        
        # Wavelength, error and data quality arrays are optional
        spectrum = Spectrum1D(data=a)
        # The wavelength array could be added later
        spectrum.wavelength = w
        self.assertTrue( np.allclose(spectrum.wavelength, np.asarray(w)) )
        del spectrum        

        # It must be possible to create a wavelength array from the
        # start and step values.
        wstart = 10.0
        wstep = 1.0
        wstop = wstart + wstep * (len(a)-1)
        spectrum = Spectrum1D(data=a, wstart=wstart, wstep=wstep)
        expected = np.linspace(wstart, wstop, len(a))
        self.assertTrue( np.allclose(spectrum.wavelength, expected) )
        del spectrum, expected
        
        # Specifying a data array that isn't 1-D should generate an error.
        bada = 42.0
        self.assertRaises(ValueError, Spectrum1D, data=bada)
        bada = [[10,15,20],[10,20,30]]
        self.assertRaises(ValueError, Spectrum1D, data=bada)

        # Specifying a wavelength array the wrong size or the wrong
        # shape should generate an error.
        badw = [1,2,3]
        self.assertRaises(ValueError, Spectrum1D, data=a, err=b, dq=c,
                          wavelength=badw)
        badw = [[1,2,3],[1,2,3]]
        self.assertRaises(ValueError, Spectrum1D, data=a, err=b, dq=c,
                          wavelength=badw)

        # It should be possible to set up an empty data product with
        # a specified shape. All the arrays should be initialised to
        # the same shape.
        # FIXME: This doesn't work. The jwst_lib DataModel fails at the
        # extract_extra_elements stage!
#         emptyspec = Spectrum1D( init=[4] )
#         self.assertIsNotNone(emptyspec.data)
#         self.assertEqual(emptyspec.data.shape, [4])
#         self.assertIsNotNone(emptyspec.err)
#         self.assertEqual(emptyspec.err.shape, [4])
#         self.assertIsNotNone(emptyspec.dq)
#         self.assertEqual(emptyspec.dq.shape, [4])
#         self.assertIsNotNone(emptyspec.wavelength)
#         self.assertEqual(emptyspec.wavelength.shape, [4])
#         descr = str(emptyspec)
#         del emptyspec, descr
 
        # It must be possible to create and display a null product
        spectrum = Spectrum1D()
        descr1 = str(spectrum)
        self.assertIsNotNone(descr1)
        # The data arrays could be added later
        spectrum.data = a
        spectrum.err = b
        spectrum.dq = c
        w = [5,7,9,11,13,15,17,19,21,23,25,27]
        spectrum.wavelength = w
        descr2 = str(spectrum)
        self.assertIsNotNone(descr2)
        del spectrum, descr1, descr2
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with Spectrum1D(self.testfile) as readback:
                self.assertTrue( ma.allclose(self.dataproduct.data,
                                             readback.data) )
                self.assertTrue( np.allclose(self.dataproduct.err,
                                             readback.err) )
                self.assertTrue( np.allclose(self.dataproduct.dq,
                                             readback.dq) )
                self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                             readback.wavelength) )
                del readback
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For this test to pass these only need to run without error.
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr


class TestSpectrum2D(unittest.TestCase):
    
    # Test the Spectrum2D class.
    # Most of the relevant tests are already done in Spectrum1D.
    # Only the additional tests specific to Spectrum2D are included here.

    def setUp(self):
        # Create a typical 2-D spectrum product.
        a = np.array([10,15,10,20,10,30,90,30,20,50,10,12])
        b = np.array([5,4,3,2,1,2,3,4,5,6,7,8])
        c = np.array([0,0,0,0,0,0,1,0,0,0,0,0])
        w = np.array([5,7,9,11,13,15,17,19,21,23,25,27])
        s = np.array([1,2,3])
        a2 = np.ascontiguousarray(np.swapaxes(np.array([a, a+10, a+2]), 0, 1))
        b2 = np.ascontiguousarray(np.swapaxes(np.array([b*2, b, b*2]), 0, 1))
        c2 = np.ascontiguousarray(np.swapaxes(np.array([c, c*0, c]), 0, 1))
        self.dataproduct = Spectrum2D(data=a2, err=b2, dq=c2, wavelength=w,
                                      slit=s, unit='photons', wunit='microns',
                                      sunit='mm')

        # Add some example metadata.
        self.dataproduct.set_observation_metadata()
        self.dataproduct.set_target_metadata(0.0, 0.0)
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriSpectrum2D_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        # Remove temporary file, if able to.
        if os.path.isfile(self.testfile):
            try:
                os.remove(self.testfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.testfile + \
                    "\n   " + str(e)
                warnings.warn(strg)

    def test_creation(self):
        # Test the creation of the data product. 
        a = np.array([10,15,10,20,10,30,90,30,20,50,10,12])
        b = np.array([5,4,3,2,1,2,3,4,5,6,7,8])
        c = np.array([0,0,0,0,0,0,1,0,0,0,0,0])
        w = np.array([5,7,9,11,13,15,17,19,21,23,25,27])
        s = np.array([1,2,3])
        a2 = np.ascontiguousarray(np.swapaxes(np.array([a, a+10, a+2]), 0, 1))
        b2 = np.ascontiguousarray(np.swapaxes(np.array([b*2, b, b*2]), 0, 1))
        c2 = np.ascontiguousarray(np.swapaxes(np.array([c, c*0, c]), 0, 1))

        # Slit, Wavelength, error and data quality arrays are optional
        spectrum = Spectrum2D(data=a2)
        # The wavelength and slit arrays could be added later
        spectrum.wavelength = w
        spectrum.slit = s
        self.assertTrue( np.allclose(spectrum.wavelength, w) )
        self.assertTrue( np.allclose(spectrum.slit, s) )
        del spectrum        

        # It must be possible to create a slit array from the
        # start and step values.
        sstart = 5.0
        sstep = 1.0
        sstop = sstart + sstep * (len(s)-1)
        spectrum = Spectrum2D(data=a2, sstart=sstart, sstep=sstep)
        expected = np.linspace(sstart, sstop, len(s))
        self.assertTrue(np.all( spectrum.slit == expected ))
        del spectrum        

        # Specifying a data array that isn't 2-D should generate an error.
        bada = 42.0
        self.assertRaises(ValueError, Spectrum2D, data=bada)
        bada = [10,15,20]
        self.assertRaises(ValueError, Spectrum2D, data=bada)

        # Specifying a slit array the wrong size or the wrong
        # shape should generate an error.
        bads = [1,2,3,4]
        self.assertRaises(ValueError, Spectrum2D, data=a2, err=b2, dq=c2,
                          slit=bads)
        bads = [[1,2,3,4],[1,2,3,4]]
        self.assertRaises(ValueError, Spectrum2D, data=a2, err=b2, dq=c2,
                          slit=bads)

        # It should be possible to set up an empty data product with
        # a specified shape. All the arrays should be initialised to
        # the same shape.
        # FIXME: This doesn't work. The jwst_lib DataModel fails at the
        # extract_extra_elements stage!
#         emptyspec = Spectrum2D( init=[4,5] )
#         self.assertIsNotNone(emptyspec.data)
#         self.assertEqual(emptyspec.data.shape, [4,5])
#         self.assertIsNotNone(emptyspec.err)
#         self.assertEqual(emptyspec.err.shape, [4,5])
#         self.assertIsNotNone(emptyspec.dq)
#         self.assertEqual(emptyspec.dq.shape, [4,5])
#         self.assertIsNotNone(emptyspec.wavelength)
#         self.assertEqual(emptyspec.wavelength.shape, [4])
#         self.assertIsNotNone(emptyspec.slit)
#         self.assertEqual(emptyspec.slit.shape, [5])
#         descr = str(emptyspec)
#         del emptyspec, descr
        
        # It must be possible to create and display a null product
        spectrum = Spectrum2D()
        descr1 = str(spectrum)
        self.assertIsNotNone(descr1)
        # The data arrays could be added later
        spectrum.data = a2
        spectrum.err = b2
        spectrum.dq = c2
        spectrum.wavelength = w
        s = np.array([1,2,3])
        spectrum.slit = s
        self.assertTrue( np.allclose(spectrum.data, np.asarray(a2)) )
        self.assertTrue( np.allclose(spectrum.err, np.asarray(b2)) )
        self.assertTrue( np.allclose(spectrum.dq, np.asarray(c2)) )
        self.assertTrue( np.allclose(spectrum.wavelength, np.asarray(w)) )
        self.assertTrue( np.allclose(spectrum.slit, np.asarray(s)) )
        descr2 = str(spectrum)
        self.assertIsNotNone(descr2)
        del spectrum, descr1, descr2
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with Spectrum2D(self.testfile) as readback:
                self.assertTrue( ma.allclose(self.dataproduct.data,
                                             readback.data) )
                self.assertTrue( np.allclose(self.dataproduct.err,
                                             readback.err) )
                self.assertTrue( np.allclose(self.dataproduct.dq,
                                             readback.dq) )
                self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                             readback.wavelength) )
                self.assertTrue( np.allclose(self.dataproduct.slit,
                                             readback.slit) )
                del readback
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For this test to pass these only need to run without error.
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        
    def test_flatten_spectrum(self):
        # Test the function which extracts a flat 1-D spectrum from
        # a 2-D spectrum. All the various permutations of the option
        # flags must work.
        flat1 = self.dataproduct.flatten_to_spectrum()
        descr1 = str(flat1)
        self.assertIsNotNone(descr1)
        del flat1,descr1

        flat2 = self.dataproduct.flatten_to_spectrum(weighted=True,
                                                     perfect=True)
        descr2 = str(flat2)
        self.assertIsNotNone(descr2)
        del flat2,descr2

        flat3 = self.dataproduct.flatten_to_spectrum(weighted=False,
                                                     perfect=False)
        descr3 = str(flat3)
        self.assertIsNotNone(descr3)
        del flat3,descr3

        flat4 = self.dataproduct.flatten_to_spectrum(weighted=False,
                                                     perfect=True)
        descr4 = str(flat4)
        self.assertIsNotNone(descr4)
        del flat4,descr4

        flat5 = self.dataproduct.flatten_to_spectrum(weighted=True,
                                                     perfect=False)
        descr5 = str(flat5)
        self.assertIsNotNone(descr5)
        
        # The resulting product should be 1-D and have the same wavelength scale
        self.assertEqual(flat5.data.ndim, 1)
        self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                    flat5.wavelength) )
        del flat5,descr5


class TestSpectrum3D(unittest.TestCase):
    
    # Test the Spectrum3D class.
    # Most of the relevant tests are already done in Spectrum1D and Spectrum2D.
    # Only the additional tests specific to Spectrum3D are included here.

    def setUp(self):
        # Create a typical 3-D spectral cube product.
        a = np.array([10,15,10,20,10,30,90,30,20,50,10,12])
        b = np.array([5,4,3,2,1,2,3,4,5,6,7,8])
        c = np.array([0,0,0,0,0,0,1,0,0,0,0,0])
        w = np.array([5,7,9,11,13,15,17,19,21,23,25,27])
        alpha = np.array([1,2,3])
        beta = np.array([2,4,6,8])
        a2 = np.ascontiguousarray(np.swapaxes(np.array([a, a+10, a+2]), 0, 1))
        b2 = np.ascontiguousarray(np.swapaxes(np.array([b*2, b, b*2]), 0, 1))
        c2 = np.ascontiguousarray(np.swapaxes(np.array([c, c*0, c]), 0, 1))
        a3 = np.ascontiguousarray(np.swapaxes(np.array([a2, a2+2, a2+4, a2+6]), 0, 1))
        b3 = np.ascontiguousarray(np.swapaxes(np.array([b2, b2+1, b2+3, b2+4]), 0, 1))
        c3 = np.ascontiguousarray(np.swapaxes(np.array([c2, c2, c2, c2]), 0, 1))

        self.dataproduct = Spectrum3D(data=a3, err=b3, dq=c3, wavelength=w,
                                      alpha=alpha, beta=beta, unit='photons',
                                      wunit='microns', aunit='arcsec',
                                      bunit='arcsec')

        # Add some example metadata.
        self.dataproduct.set_observation_metadata()
        self.dataproduct.set_target_metadata(0.0, 0.0)
        self.dataproduct.set_instrument_metadata(detector='MIRIFUSHORT',
                                                 channel='1',
                                                 ccc_pos='OPEN',
                                                 deck_temperature=11.0,
                                                 detector_temperature=6.0)
        self.dataproduct.set_exposure_metadata(readpatt='FAST',
                                               nints=1, ngroups=1,
                                               frame_time=1.0,
                                               integration_time=10.0,
                                               group_time=10.0,
                                               reset_time=0, frame_resets=3)
        self.testfile = "MiriSpectrum3D_test.fits"
        
    def tearDown(self):
        # Tidy up
        del self.dataproduct
        # Remove temporary file, if able to.
        if os.path.isfile(self.testfile):
            try:
                os.remove(self.testfile)
            except Exception as e:
                strg = "Could not remove temporary file, " + self.testfile + \
                    "\n   " + str(e)
                warnings.warn(strg)

    def test_creation(self):
        # Test the creation of the data product. 
        a = np.array([10,15,10,20,10,30,90,30,20,50,10,12])
        b = np.array([5,4,3,2,1,2,3,4,5,6,7,8])
        c = np.array([0,0,0,0,0,0,1,0,0,0,0,0])
        w = np.array([5,7,9,11,13,15,17,19,21,23,25,27])
        alpha = np.array([1,2,3])
        beta = np.array([2,4,6])
        a2 = np.ascontiguousarray(np.swapaxes(np.array([a, a+10, a+2]), 0, 1))
        b2 = np.ascontiguousarray(np.swapaxes(np.array([b*2, b, b*2]), 0, 1))
        c2 = np.ascontiguousarray(np.swapaxes(np.array([c, c*0, c]), 0, 1))
        a3 = np.ascontiguousarray(np.swapaxes(np.array([a2, a2+2, a2+4]), 0, 1))
        b3 = np.ascontiguousarray(np.swapaxes(np.array([b2, b2+1, b2+3]), 0, 1))
        c3 = np.ascontiguousarray(np.swapaxes(np.array([c2, c2, c2]), 0, 1))

        # Alpha, Beta, Wavelength, error and data quality arrays are optional
        cube = Spectrum3D(data=a3)
        # The wavelength, alpha and beta arrays could be added later
        cube.wavelength = w
        cube.alpha = alpha
        cube.beta = beta
        self.assertTrue( np.allclose(cube.wavelength, w) )
        self.assertTrue( np.allclose(cube.alpha, alpha) )
        self.assertTrue( np.allclose(cube.beta, beta) )
        del cube        

        # It must be possible to create alpha and beta arrays from the
        # start and step values.
        astart = 5.0
        astep = 1.0
        astop = astart + astep * (len(alpha)-1)
        bstart = 0.0
        bstep = 0.1
        bstop = bstart + bstep * (len(beta)-1)
        cube = Spectrum3D(data=a3, astart=astart, astep=astep, bstart=bstart,
                          bstep=bstep)
        aexpected = np.linspace(astart, astop, len(alpha))
        self.assertTrue( np.allclose( cube.alpha, aexpected ) )
        bexpected = np.linspace(bstart, bstop, len(beta))
        self.assertTrue( np.allclose( cube.beta, bexpected ) )
        del cube, aexpected, bexpected

        # Specifying a data array that isn't 3-D should generate an error.
        bada = 42.0
        self.assertRaises(ValueError, Spectrum3D, data=bada)
        bada = [10,15,20]
        self.assertRaises(ValueError, Spectrum3D, data=bada)
        bada = [[10,15,20],[10,20,30]]
        self.assertRaises(ValueError, Spectrum3D, data=bada)

        # Specifying alpha or beta arrays the wrong size or the wrong
        # shape should generate an error.
        alpha = [1,2,3,4]
        self.assertRaises(ValueError, Spectrum3D, data=a2, err=b2, dq=c2,
                          alpha=alpha)
        alpha = [[1,2,3,4],[1,2,3,4]]
        self.assertRaises(ValueError, Spectrum3D, data=a2, err=b2, dq=c2,
                          alpha=alpha)
        beta = [1,2,3,4]
        self.assertRaises(ValueError, Spectrum3D, data=a2, err=b2, dq=c2,
                          beta=beta)
        beta = [[1,2,3,4],[1,2,3,4]]
        self.assertRaises(ValueError, Spectrum3D, data=a2, err=b2, dq=c2,
                          beta=beta)

        # It should be possible to set up an empty data product with
        # a specified shape. All the arrays should be initialised to
        # the same shape.
        # FIXME: This doesn't work. The jwst_lib DataModel fails at the
        # extract_extra_elements stage!
#         emptycube = Spectrum3D( init=[2,3,4] )
#         self.assertIsNotNone(emptycube.data)
#         self.assertEqual(emptycube.data.shape, [2,3,4])
#         self.assertIsNotNone(emptycube.err)
#         self.assertEqual(emptycube.err.shape, [2,3,4])
#         self.assertIsNotNone(emptycube.dq)
#         self.assertEqual(emptycube.dq.shape, [2,3,4])
#         self.assertIsNotNone(emptycube.wavelength)
#         self.assertEqual(emptycube.wavelength.shape, [2])
#         self.assertIsNotNone(emptycube.alpha)
#         self.assertEqual(emptycube.alpha.shape, [3])
#         self.assertIsNotNone(emptycube.beta)
#         self.assertEqual(emptycube.beta.shape, [4])
#         descr = str(emptycube)
#         del emptycube, descr
       
        # It must be possible to create and display a null product
        cube = Spectrum3D()
        descr1 = str(cube)
        self.assertIsNotNone(descr1)
        # The data arrays could be added later
        cube.data = a3
        cube.err = b3
        cube.dq = c3
        cube.wavelength = w
        alpha = np.array([1,2,3])
        beta = np.array([2,4,6])
        cube.alpha = alpha
        cube.beta = beta
        self.assertTrue( np.allclose(cube.data, np.asarray(a3)) )
        self.assertTrue( np.allclose(cube.err, np.asarray(b3)) )
        self.assertTrue( np.allclose(cube.dq, np.asarray(c3)) )
        self.assertTrue( np.allclose(cube.wavelength, np.asarray(w)) )
        self.assertTrue( np.allclose(cube.alpha, np.asarray(alpha)) )
        self.assertTrue( np.allclose(cube.beta, np.asarray(beta)) )
        descr2 = str(cube)
        self.assertIsNotNone(descr2)
        del cube, descr1, descr2
        
    def test_fitsio(self):
        # Suppress metadata warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Check that the data product can be written to a FITS
            # file and read back again without changing the data.
            self.dataproduct.save(self.testfile, overwrite=True)
            with Spectrum3D(self.testfile) as readback:
                self.assertTrue( ma.allclose(self.dataproduct.data,
                                             readback.data) )
                self.assertTrue( np.allclose(self.dataproduct.err,
                                             readback.err) )
                self.assertTrue( np.allclose(self.dataproduct.dq,
                                             readback.dq) )
                self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                             readback.wavelength) )
                self.assertTrue( np.allclose(self.dataproduct.alpha,
                                             readback.alpha) )
                self.assertTrue( np.allclose(self.dataproduct.beta,
                                             readback.beta) )
                del readback
        
    def test_description(self):
        # Test that the querying and description functions work.
        # For this test to pass these only need to run without error.
        descr = str(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.dataproduct)
        self.assertIsNotNone(descr)
        del descr

    def test_flatten_spectrum(self):
        # Test the function which extracts a flat 1-D spectrum from
        # a 3-D spectral cube. All the various permutations of the option
        # flags must work.
        flat1 = self.dataproduct.flatten_to_spectrum()
        descr1 = str(flat1)
        self.assertIsNotNone(descr1)
        del flat1,descr1

        flat2 = self.dataproduct.flatten_to_spectrum(weighted=True,
                                                     perfect=True)
        descr2 = str(flat2)
        self.assertIsNotNone(descr2)
        del flat2,descr2

        flat3 = self.dataproduct.flatten_to_spectrum(weighted=False,
                                                     perfect=False)
        descr3 = str(flat3)
        self.assertIsNotNone(descr3)
        del flat3,descr3

        flat4 = self.dataproduct.flatten_to_spectrum(weighted=False,
                                                     perfect=True)
        descr4 = str(flat4)
        self.assertIsNotNone(descr4)
        del flat4,descr4

        flat5 = self.dataproduct.flatten_to_spectrum(weighted=True,
                                                     perfect=False)
        descr5 = str(flat5)
        self.assertIsNotNone(descr5)
        
        # The resulting product should be 1-D and have the same wavelength scale
        self.assertEqual(flat5.data.ndim, 1)
        self.assertTrue( np.allclose(self.dataproduct.wavelength,
                                    flat5.wavelength) )
        del flat5,descr5

    def test_flatten_image(self):
        # Test the function which extracts a 2-D image from
        # a 3-D spectral cube. All the various permutations of the option
        # flags must work.
        flat1 = self.dataproduct.flatten_to_image()
        descr1 = str(flat1)
        self.assertIsNotNone(descr1)
        del flat1,descr1

        flat2 = self.dataproduct.flatten_to_image(weighted=True, perfect=True)
        descr2 = str(flat2)
        self.assertIsNotNone(descr2)
        del flat2,descr2

        flat3 = self.dataproduct.flatten_to_image(weighted=False, perfect=False)
        descr3 = str(flat3)
        self.assertIsNotNone(descr3)
        del flat3,descr3

        flat4 = self.dataproduct.flatten_to_image(weighted=False, perfect=True)
        descr4 = str(flat4)
        self.assertIsNotNone(descr4)
        del flat4,descr4

        flat5 = self.dataproduct.flatten_to_image(weighted=True, perfect=False)
        descr5 = str(flat5)
        self.assertIsNotNone(descr5)
        
        # The resulting image should be 2-D of dimensions (alpha,beta)
        self.assertEqual(flat5.data.ndim, 2)
        self.assertEqual(flat5.data.shape[1], len(self.dataproduct.alpha))
        self.assertEqual(flat5.data.shape[0], len(self.dataproduct.beta))
        del flat5,descr5


# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
