#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Some extensions to the standard STScI data model which define
MIRI photometric flux conversion models.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

21 Aug 2015: Created from miri_fluxconversion_models.
10 Sep 2015: Only implement the imaging and generic models. The format of the
             MRS model is TBD.
11 Sep 2015: Added MiriPixelAreaModel class.
06 Nov 2015: Corrected typo in definition of ARCSEC2_PER_STERADIAN.
17 Nov 2015: SUBARRAY is GENERIC when photometric factors do not depend
             on subarray.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
             Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import math
import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasData

# List all classes and global functions here.
__all__ = ['MiriPhotometricModel', 'MiriImagingPhotometricModel',
           'MiriPixelAreaModel']

# Useful constants
ARCSEC2_PER_STERADIAN = ((180.0/math.pi) * 3600.0)**2 # Square arcsec per steradian
MAX_NLEM = 500 # Maximum size of the wavelength and relresponse arrays

class MiriPhotometricModel(MiriDataModel):
    """
    
    A data model for MIRI photometric flux conversion data.
    Flux factors are looked up by filter name and subarray.
    A wavelength and response array can also be provided for
    each filter and subarray combination.
    
    The model is valid for MIRI post-CDP-4 imaging photometry data.
    It may also be used for MIRI post-CDP-4 LRS photometry data,
    although the format for spectroscopic data might change.
    The data format for MIRI MRS photometry data is TBD.
    
    :Parameters:
    
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the phot_table parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    phot_table. list of tuples or numpy record array (optional)
        A list of tuples, or a numpy record array, where each record
        contains the following information:
        
        * FILTER: A string containing a filter name. Compulsory.
        * SUBARRAY: A string containing a subarray name. The string can
          be 'GENERIC' to indicate that the conversion factor does not
          depend on the subarray.
        * PHOTMJSR: A conversion factor for the given filter and subarray
          combination, converting from DN/s/pixel to MJy/Sr, assuming that
          each pixel has the same area on the sky.
        * NELM: An integer (range 0 to 500) giving the number of elements
          in the wavelength and relresponse arrays to be used. If NELM=0,
          the wavelength and response arrays are ignored.
        * WAVELENGTH: An array of 500 floats containing wavelengths in microns.
          Unused parts of the array are padded with zero.
        * RELRESPONSE: An array of 500 floats containing the relative
          response of the instrument at each wavelength.
          Unused parts of the array are padded with zero.

        A phot_table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
        
    pixar_sr: number (optional)
        The nominal pixel area for the detector in steradians.
        If provided, the value is written to the PIXAR_SR keyword.
    pixel_a2: number (optional)
        The nominal pixel area for the detector in square arcseconds
        If provided, the value is written to the PIXAR_A2 keyword.
    
    """
    schema_url = "miri_photom.schema.yaml"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse')

    def __init__(self, init=None, phot_table=None, pixar_sr=None,
                 pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriPhotometricModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPhotometricModel, self).__init__(init=init, **kwargs)
        # Data type is photometric flux conversion.
        self.meta.model_type = 'PHOTOM'
        self.meta.reftype = 'PHOTOM'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'

        if phot_table is not None:
            try:
                #phot_table = np.recarray(phot_table)
                self.phot_table = phot_table
            except (ValueError, TypeError) as e:
                strg = "phot_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # If provided, define the pixel area metadata.
        if pixar_sr is not None:
            self.meta.photometry.pixelarea_steradians = pixar_sr
            # Should both keywords be written?
#             self.meta.photometry.arcsecsq = pixar_s2 * ARCSEC2_PER_STERADIAN
        if pixar_a2 is not None:
            self.meta.photometry.pixelarea_arcsecsq = pixar_a2
            # Should both keywords be written?
#             self.meta.photometry.pixelarea_steradians = \
#                 pixar_s2 / ARCSEC2_PER_STERADIAN


class MiriImagingPhotometricModel(MiriPhotometricModel):
    """
    
    A data model for MIRI imaging mode photometric flux conversion data.
    Flux factors are looked up by filter name and subarray.
    
    This alternative call to the photometric model initialises the
    wavelength and response arrays to zero, so they don't need to
    be defined by the caller.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the phot_table parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    phot_table. list of tuples or numpy record array (optional)
        A list of tuples, or a numpy record array, where each record
        contains the following information:
        
        * FILTER: A string containing a filter name. Compulsory.
        * SUBARRAY: A string containing a subarray name. The string can
          be 'GENERIC' to indicate that the conversion factor does not
          depend on the subarray. Empty strings are converted to 'GENERIC'.
        * PHOTMJSR: A conversion factor for the given filter and subarray
          combination, converting from DN/s/pixel to MJy/Sr, assuming that
          each pixel has the same area on the sky.

        All the NELM, WAVELENGTH and RELESPONSE columns are automatically
        initialised to zero.
        A phot_table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
        
    pixar_sr: number (optional)
        The nominal pixel area for the detector in steradians.
        If provided, the value is written to the PIXAR_SR keyword.
    pixel_a2: number (optional)
        The nominal pixel area for the detector in square arcseconds
        If provided, the value is written to the PIXAR_A2 keyword.
    
    """
    # Both models use exactly the same schema.
    schema_url = "miri_photom.schema.yaml"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse')

    def __init__(self, init=None, phot_table=None, pixar_sr=None,
                 pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriImagingPhotometricModel class.
        
        Parameters: See class doc string.

        """
        if phot_table is not None:
            # Construct a phot_table from the given phot_table_img and
            # some dummy wavelength and relreponse arrays.
            wavelength = [0.0] * MAX_NLEM
            relresponse = [0.0] * MAX_NLEM
            new_phot_table = []
            for (filter, subarray, photmjsr, uncertainty) in phot_table:
                if not subarray:
                    # Convert an empty SUBARRAY string to 'GENERIC'
                    subarray = 'GENERIC'
                new_phot_table.append( (filter, subarray, photmjsr, uncertainty,
                                    0, wavelength, relresponse) )
        else:
            new_phot_table = None

        # Pass the phot_table to the generic model.
        super(MiriImagingPhotometricModel, self).__init__(init=init,
                                                          phot_table=new_phot_table,
                                                          pixar_sr=pixar_sr,
                                                          pixar_a2=pixar_a2,
                                                          **kwargs)
        self.meta.model_type = 'PHOTOM (Imaging)'
        self.add_comment("WAVELENGTH and RELRESPONSE arrays are all zero.")


class MiriPixelAreaModel(MiriDataModel, HasData):
    """
    
    A data model for MIRI pixel area data.
        
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the mask parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    data: numpy array (optional)
        An array containing the pixel area data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    pixar_sr: number (optional)
        The nominal pixel area for the detector in steradians.
        If provided, the value is written to the PIXAR_SR keyword.
    pixel_a2: number (optional)
        The nominal pixel area for the detector in square arcseconds
        If provided, the value is written to the PIXAR_A2 keyword.
    
    """
    schema_url = "miri_pixelarea.schema.yaml"

    def __init__(self, init=None, data=None, pixar_sr=None, pixar_a2=None,
                 **kwargs):
        """
        
        Initialises the MiriPixelAreaModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPixelAreaModel, self).__init__(init=init, **kwargs)


        # Data type is AREA.
        self.meta.model_type = 'AREA (Pixel Area)'
        self.meta.reftype = 'AREA'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'

        # If provided, define the pixel area metadata.
        if pixar_sr is not None:
            self.meta.photometry.pixelarea_steradians = pixar_sr
            # Should both keywords be written?
#             self.meta.photometry.arcsecsq = pixar_s2 * ARCSEC2_PER_STERADIAN
        if pixar_a2 is not None:
            self.meta.photometry.pixelarea_arcsecsq = pixar_a2
            # Should both keywords be written?
#             self.meta.photometry.pixelarea_steradians = \
#                 pixar_s2 / ARCSEC2_PER_STERADIAN

        # Update the data array if it has been specifically provided.
        HasData.__init__(self, data)


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MIRI photometric flux conversion module.")
    
    import sys
    PLOTTING = False
    SAVE_FILES = False
    
    pixar_a2 = 0.11 * 0.11 # MIRI imager pixels are 0.1 arcsec1 x 0.11 arcsec
    pixar_sr = pixar_a2 / ARCSEC2_PER_STERADIAN
    
    # Imager photometric model. Dummy response values and blank (zero-filled)
    # wavelength and response arrays
    wavelength = [0.0] * MAX_NLEM
    relresponse = [0.0] * MAX_NLEM
    # There is a separate response record for each filter, valid for all
    # subarrays.
    phot_im1 = [('F560W',  'GENERIC', 2.41,  0.26,  0, wavelength, relresponse),
               ('F770W',   'GENERIC', 1.32,  0.013, 0, wavelength, relresponse),
               ('F1000W',  'GENERIC', 1.76,  0.12,  0, wavelength, relresponse),
               ('F1130W',  'GENERIC', 5.76,  0.43,  0, wavelength, relresponse),
               ('F1280W',  'GENERIC', 2.11,  0.16,  0, wavelength, relresponse),
               ('F1500W',  'GENERIC', 1.84,  0.01,  0, wavelength, relresponse),
               ('F1800W',  'GENERIC', 2.68,  0.23,  0, wavelength, relresponse),
               ('F2100W',  'GENERIC', 2.04,  0.15,  0, wavelength, relresponse),
               ('F2550W',  'GENERIC', 4.25,  0.4,   0, wavelength, relresponse),
               ('F2550WR', 'GENERIC', 4.60,  0.24,  0, wavelength, relresponse),
               ('F1065C',  'GENERIC', 1.37,  0.1,   0, wavelength, relresponse),
               ('F1140C',  'GENERIC', 1.43,  0.11,  0, wavelength, relresponse),
               ('F1550C',  'GENERIC', 1.81,  0.13,  0, wavelength, relresponse),
               ('F2300C',  'GENERIC', 3.65,  0.23,  0, wavelength, relresponse),
               ]

    # An alternative way of defining an imaging photometric model, without
    # needing to fill in the dummy wavelength and relresponse arrays.
    # Empty subarray fields are converted to 'GENERIC'.
    phot_im2 = [('F560W',  '', 2.41,  0.26),
               ('F770W',   '', 1.32,  0.013),
               ('F1000W',  '', 1.76,  0.12),
               ('F1130W',  '', 5.76,  0.43),
               ('F1280W',  '', 2.11,  0.16),
               ('F1500W',  '', 1.84,  0.01),
               ('F1800W',  '', 2.68,  0.23),
               ('F2100W',  '', 2.04,  0.15),
               ('F2550W',  '', 4.25,  0.4),
               ('F2550WR', '', 4.60,  0.24),
               ('F1065C',  '', 1.37,  0.1),
               ('F1140C',  '', 1.43,  0.11),
               ('F1550C',  '', 1.81,  0.13),
               ('F2300C',  '', 3.65,  0.23),
               ]

    print("\nImaging PHOTOM model (with provided blank wavelength array):")
    with MiriPhotometricModel( phot_table=phot_im1, pixar_a2=pixar_a2,
                               pixar_sr=pixar_sr ) as testphot1:
        testphot1.set_instrument_metadata(detector='MIRIMAGE',
                                          ccc_pos='OPEN',
                                          deck_temperature=11.0,
                                          detector_temperature=6.0)
        testphot1.set_exposure_metadata(readpatt='FAST',
                                        nints=1, ngroups=1,
                                        frame_time=1.0,
                                        integration_time=10.0,
                                        group_time=10.0,
                                        reset_time=0, frame_resets=3)
        testphot1.set_subarray_metadata('FULL')
        testphot1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testphot1)
        if PLOTTING:
            testphot1.plot(description="testphot1")
        if SAVE_FILES:
            testphot1.save("test_img_photom_model1.fits", overwrite=True)
        del testphot1

    print("\nImaging PHOTOM model (without needing a wavelength array):")
    with MiriImagingPhotometricModel( phot_table=phot_im2, pixar_a2=pixar_a2,
                                      pixar_sr=pixar_sr ) as testphot2:
        testphot2.set_instrument_metadata(detector='MIRIMAGE',
                                          ccc_pos='OPEN',
                                          deck_temperature=11.0,
                                          detector_temperature=6.0)
        testphot2.set_exposure_metadata(readpatt='FAST',
                                        nints=1, ngroups=1,
                                        frame_time=1.0,
                                        integration_time=10.0,
                                        group_time=10.0,
                                        reset_time=0, frame_resets=3)
        testphot2.set_subarray_metadata('FULL')
        testphot2.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testphot2)
        if PLOTTING:
            testphot2.plot(description="testphot2")
        if SAVE_FILES:
            testphot2.save("test_img_photom_model2.fits", overwrite=True)
        del testphot2

    # LRS photometric model. Dummy response values and dummy wavelength and response arrays.
    wavelength = []
    relresponse = []
    resp = 0
    nelm = 0
    for ii in range(300):
        wav = float(ii)/12.0
        resp = (resp + 1) % 30
        r10 = 0.1 + resp/35.0
        wavelength.append(wav)
        relresponse.append(r10)
        nelm += 1
    for ii in range(300,MAX_NLEM):
        # Pad unused elements with zero.
        wavelength.append(0.0)
        relresponse.append(0.0)
    # The LRS data is for the P750L filter only, but there is separate
    # record for the SLITLESSPRISM subarray.
    phot_lrs = [('P750L',  '',              1.0,  0.0,  nelm, wavelength, relresponse),
                ('P750L',  'SLITLESSPRISM', 0.9,  0.0,  nelm, wavelength, relresponse)
               ]

    print("\nLRS PHOTOM with defined wavelength and relresponse arrays:")
    with MiriPhotometricModel( phot_table=phot_lrs, pixar_a2=pixar_a2,
                               pixar_sr=pixar_sr ) as testphot3:
        testphot3.set_instrument_metadata(detector='MIRIMAGE',
                                          ccc_pos='OPEN', filt='P750L',
                                          deck_temperature=11.0,
                                          detector_temperature=6.0)
        testphot3.set_exposure_metadata(readpatt='FAST',
                                        nints=1, ngroups=1,
                                        frame_time=1.0,
                                        integration_time=10.0,
                                        group_time=10.0,
                                        reset_time=0, frame_resets=3)
        testphot3.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testphot3)
        if PLOTTING:
            testphot3.plot(description="testphot3")
        if SAVE_FILES:
            testphot3.save("test_lrs_photom_model3.fits", overwrite=True)
        del testphot3

    print("Pixel area model:")
    
    # Generate a test pattern
    data5x5 = np.ones([5,5])
    for ii in range(0,5):
        for jj in range(0,5):
            d = abs(float(ii)-2.0) + abs(float(jj)-2.0)
            data5x5[ii,jj] += d/10.0
    # Ensure the mean value is 1.0
    dmean = np.mean(data5x5)
    data5x5 = data5x5 / dmean
    
    with MiriPixelAreaModel( data=data5x5, pixar_a2=pixar_a2,
                             pixar_sr=pixar_sr ) as testdata:
        testdata.set_instrument_metadata(detector='MIRIMAGE',
                                        ccc_pos='OPEN',
                                        deck_temperature=11.0,
                                        detector_temperature=6.0)
        testdata.set_exposure_metadata(readpatt='FAST',
                                       nints=1, ngroups=1,
                                       frame_time=1.0,
                                       integration_time=10.0,
                                       group_time=10.0,
                                       reset_time=0, frame_resets=3)
        testdata.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_pixelarea_model1.fits", overwrite=True)
        del testdata
       
    print("Test finished.")
