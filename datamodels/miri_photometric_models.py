#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Some extensions to the standard STScI data model which define
MIRI photometric flux conversion models.

NOTE: There is no MRS version of this data model.
See MiriMrsFluxconversionModel within miri_fluxconversion_models.py

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

21 Aug 2015: Created from miri_fluxconversion_models.
10 Sep 2015: Only implement the imaging and generic models. The format of the
             MRS model is TBD.
11 Sep 2015: Added MiriPixelAreaModel class.
06 Nov 2015: Corrected typo in definition of ARCSEC2_PER_STERADIAN.
17 Nov 2015: SUBARRAY is GENERIC when photometric factors do not depend on
             subarray.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the JWST
             build 7.1 data models release. meta.reffile.type also changed to
             meta.reftype. TYPE keyword replaced by DATAMODL. Do not set
             observation or target metadata. Neither are appropriate for a
             reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
04 Oct 2018: Define exposure type.
10 Oct 2018: Added MiriLrsPhotometricModel class, append and append_srf methods.
             Corrected __all__.
17 Oct 2018: Added relresperror column to MiriPhotometricModel. Old data model
             preserved as MiriPhotometricModel_CDP5.
26 Oct 2018: Added get_srf() function.
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Steven Beard (UKATC)

"""

import warnings
import math
import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasData

# List all classes and global functions here.
__all__ = ['MiriPhotometricModel', 'MiriImagingPhotometricModel', 'MiriPhotometricModel_CDP5',
           'MiriLrsPhotometricModel', 'MiriPixelAreaModel', 'MAX_NELEM']

# Useful constants
ARCSEC2_PER_STERADIAN = ((180.0/math.pi) * 3600.0)**2 # Square arcsec per steradian
MAX_NELEM = 500 # Maximum size of the wavelength and relresponse arrays

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
    schema_url = "miri_photom.schema"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse', 'relresperror')

    def __init__(self, init=None, phot_table=None, pixar_sr=None,
                 pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriPhotometricModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPhotometricModel, self).__init__(init=init, **kwargs)

        # Data type is photometric flux conversion.
        self.meta.reftype = 'PHOTOM'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()

        if phot_table is not None:
            try:
                #phot_table = np.recarray(phot_table)
                self.phot_table = phot_table
            except (ValueError, TypeError) as e:
                strg = "phot_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
            
        # Copy the table column units from the schema, if defined.
        phot_table_units = self.set_table_units('phot_table')

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

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def append(self, other):
        """
        
        Append the phot_table from another MiriPhotometricModel to the
        phot_table of this model. New rows are added at the end of the
        existing table.
         
        :Parameters:
        
        other: MiriPhotometricModel
            Other data model to be appended to this one.
            
        :Attributes:
        
        phot_table: numpy recarray
            New rows appended.
       
        """
        #import numpy.lib.recfunctions as rfn
        assert isinstance(other, MiriPhotometricModel,
                "append: must provide a MiriPhotometricModel object")
        # TODO: A slow, brute force method. There is probably a faster way,
        # but the function is not time critical for small CDPs..
        new_phot_table = []
        for this_row in self.phot_table:
            new_phot_table.append( tuple(this_row) )
        for other_row in other.phot_table:
            new_phot_table.append( tuple(other_row) )
        self.phot_table = new_phot_table

    def append_lrs(self, *args):
        """
        
        Append the SRF tables contained in one or more MiriLrsFluxconversionModel
        data models to the phot_table of this model. Each data mode given
        results in a new row added at the end of the existing table.
        
        NOTE: This function has a variable number of parameters. Each
        parameter is assumed to be a MiriLrsFluxconversionModel to be
        appended.
        
        :Parameters:
            
        lrs_model1: MiriLrsFluxconversionModel
            A MiriLrsFluxconversionModel data model to be appended
            to this one.
        lrs_model2: MiriLrsFluxconversionModel
            A MiriLrsFluxconversionModel data model to be appended
            to this one.
        etc...
            
        :Attributes:
        
        phot_table: numpy recarray
            New row appended.
        
        """
        # The function only works is at least one MiriLrsFluxconversionModel
        # has been specified
        if len(args) > 0:
            # Initialise a new phot_table.
            new_phot_table = []
            for this_row in self.phot_table:
                new_phot_table.append( tuple(this_row) )
                
            # Append each LRS data model one at a time.
            # All LRS models use the same filter and define an absolute
            # response, so PHOTMJSR is 1.0 and UNCERTAINTY is 0.0.
            mirifilter = 'P750L'
            photmjsr = 1.0
            uncertainty = 0.0
            argnum = 0
            for lrs_model in args:
                argnum += 1
                if hasattr( lrs_model, 'flux_table'):
                    if hasattr(lrs_model, 'meta') and hasattr(lrs_model.meta, 'subarray'):
                        subarray = lrs_model.meta.subarray.name
                    else:
                        subarray = 'GENERIC'
                    # Add the SRF flux table to the data model.
                    # TODO: A slow but readable method. There is probably a
                    # faster way, but the function is not time critical for
                    # small CDPs.
                    wavelength = [0.0] * MAX_NELEM
                    relresponse = [0.0] * MAX_NELEM
                    relresperror = [0.0] * MAX_NELEM
                    ii = 0
                    for (wav, srf, unc) in lrs_model.flux_table:
                        wavelength[ii] = wav
                        relresponse[ii] = srf
                        relresperror[ii] = unc
                        ii += 1
                    nelem = len(lrs_model.flux_table)
                    new_phot_table.append( (mirifilter, subarray, photmjsr, uncertainty,
                                      nelem, tuple(wavelength), tuple(relresponse),
                                      tuple(relresperror)))
                else:
                    strg = "Function argument %d is not a MiriFluxconversionModel" % argnum
                    raise TypeError(strg)
            self.phot_table = new_phot_table
            
    def get_srf(self, filt, subarray):
        """
        
        Return the SRF arrays associated with the given filter and subarray.
        
        This function is mainly designed for LRS data where filt='P750L'.
        
        :Parameters:
            
        filt: str
            The MIRI filter whose SRF arrays are to be returned.
            Use file='P750L' to obtain the LRS SRF arrays.
        subarray: str
            The subarray whose SRF arrays are to be returned.
            Use subarray='FULL' or 'SLITLESSPRISM' to obtain the LRS SRF arrays.
            
        :Returned:
        
        srf_arrays: tuple of (wavelength, srf, uncertainty)
            The wavelength, srf and uncertainty arrays
        
        """
        # Find the relevant row in the PHOTOM table
        for row in self.phot_table:
            if str(row[0]).strip() == filt and str(row[1]).strip() == subarray:
                # A matching row has been found
                nelem = int(row[4])
                if nelem > 0:
                    return( row[5][:nelem], row[6][:nelem], row[7][:nelem] )
                else:
                    strg = "Table row for filter=\'%s\' subarray=\'%s\'" % (filt,subarray)
                    strg += " has no spectral response data "
                    warnings.warn(strg)
                    return ([], [], [])
        # No matching filter
        strg = "No table row matching filter=\'%s\' subarray=\'%s\'" % (filt,subarray)
        warnings.warn(strg)
        return (None, None, None)

class MiriPhotometricModel_CDP5(MiriDataModel):
    """
    
    This class can be used to access the old CDP-6 version of the
    MiriPhotometricModel data model.
    
    See the MiriPhotometricModel class for full documentation.
    
    """
    schema_url = "miri_photom_CDP5.schema"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse')

    def __init__(self, init=None, phot_table=None, pixar_sr=None,
                 pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriPhotometricModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPhotometricModel_CDP5, self).__init__(init=init, **kwargs)

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

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def append(self, other):
        """
        
        Append the phot_table from another MiriPhotometricModel to the
        phot_table of this model. New rows are added at the end of the
        existing table.
         
        :Parameters:
        
        other: MiriPhotometricModel
            Other data model to be appended to this one.
            
        :Attributes:
        
        phot_table: numpy recarray
            New rows appended.
       
        """
        #import numpy.lib.recfunctions as rfn
        assert isinstance(other, MiriPhotometricModel,
                "append: must provide a MiriPhotometricModel object")
        # TODO: A slow, brute force method. There is probably a faster way,
        # but the function is not time critical for small CDPs..
        new_phot_table = []
        for this_row in self.phot_table:
            new_phot_table.append( tuple(this_row) )
        for other_row in other.phot_table:
            new_phot_table.append( tuple(other_row) )
        self.phot_table = new_phot_table

    def append_lrs(self, *args):
        """
        
        Append the SRF tables contained in one or more MiriLrsFluxconversionModel
        data models to the phot_table of this model. Each data mode given
        results in a new row added at the end of the existing table.
        
        NOTE: This function has a variable number of parameters. Each
        parameter is assumed to be a MiriLrsFluxconversionModel to be
        appended.
        
        :Parameters:
            
        lrs_model1: MiriLrsFluxconversionModel
            A MiriLrsFluxconversionModel data model to be appended
            to this one.
        lrs_model2: MiriLrsFluxconversionModel
            A MiriLrsFluxconversionModel data model to be appended
            to this one.
        etc...
            
        :Attributes:
        
        phot_table: numpy recarray
            New row appended.
        
        """
        # The function only works is at least one MiriLrsFluxconversionModel
        # has been specified
        if len(args) > 0:
            # Initialise a new phot_table.
            new_phot_table = []
            for this_row in self.phot_table:
                new_phot_table.append( tuple(this_row) )
                
            # Append each LRS data model one at a time.
            # All LRS models use the same filter and define an absolute
            # response, so PHOTMJSR is 1.0 and UNCERTAINTY is 0.0.
            mirifilter = 'P750L'
            photmjsr = 1.0
            uncertainty = 0.0
            argnum = 0
            for lrs_model in args:
                argnum += 1
                if hasattr( lrs_model, 'flux_table'):
                    if hasattr(lrs_model, 'meta') and hasattr(lrs_model.meta, 'subarray'):
                        subarray = lrs_model.meta.subarray.name
                    else:
                        subarray = 'GENERIC'
                    # Add the SRF flux table to the data model.
                    # TODO: A slow but readable method. There is probably a
                    # faster way, but the function is not time critical for
                    # small CDPs.
                    wavelength = [0.0] * MAX_NELEM
                    relresponse = [0.0] * MAX_NELEM
                    ii = 0
                    for (wav, srf, unc) in lrs_model.flux_table:
                        wavelength[ii] = wav
                        relresponse[ii] = srf
                        ii += 1
                    nelem = len(lrs_model.flux_table)
                    new_phot_table.append( (mirifilter, subarray, photmjsr, uncertainty,
                                      nelem, tuple(wavelength), tuple(relresponse)))
                else:
                    strg = "Function argument %d is not a MiriFluxconversionModel" % argnum
                    raise TypeError(strg)
            self.phot_table = new_phot_table


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
    schema_url = "miri_photom.schema"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse', 'relresperror')

    def __init__(self, init=None, phot_table=None, pixar_sr=None,
                 pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriImagingPhotometricModel class.
        
        Parameters: See class doc string.

        """
        if phot_table is not None:
            # Construct a new phot_table from the given phot_table and
            # some dummy wavelength and relreponse arrays.
            wavelength = [0.0] * MAX_NELEM
            relresponse = [0.0] * MAX_NELEM
            relresperror = [0.0] * MAX_NELEM
            new_phot_table = []
            for phot_row in phot_table:
                mirifilter = phot_row[0]
                subarray = phot_row[1]
                photmjsr = phot_row[2]
                uncertainty = phot_row[3]
                if not subarray:
                    # Convert an empty SUBARRAY string to 'GENERIC'
                    subarray = 'GENERIC'
                new_phot_table.append( (mirifilter, subarray, photmjsr, uncertainty,
                                    0, wavelength, relresponse, relresperror) )
        else:
            new_phot_table = None

        # Pass the phot_table to the generic model.
        super(MiriImagingPhotometricModel, self).__init__(init=init,
                                                          phot_table=new_phot_table,
                                                          pixar_sr=pixar_sr,
                                                          pixar_a2=pixar_a2,
                                                          **kwargs)
        #self.add_comment("WAVELENGTH and RELRESPONSE arrays are all zero for imager.")
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type

class MiriLrsPhotometricModel(MiriPhotometricModel):
    """
    
    A data model for MIRI LRS mode photometric flux conversion data.
    Flux factors are looked up by filter name and subarray.
    
    This alternative call to the photometric model initialises the
    filter, photmjsr and uncertainty arrays, so they don't
    need to be defined by the caller.
    
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
        
    subarray: str
        Subarray for which the following srf_table is valid. Compulsory
        if srf_table provided. Expected to be 'FULL' or 'SLITLESSPRISM'
        
    srf_table. list of tuples or numpy record array (optional)
        A list of tuples, or a numpy record array, where each record
        contains the following information:
        
        * WAVELENGTH: Wavelength.
        * SRF: Spectral response function.
        * UNCERTAINTY: Uncertainty in the spectral response function.

        The UNCERTAINTY data provided in the srf_table is ignored.
        All the FILTER, PHOTMJSR and UNCERTAINTY columns in the phot_table
        are automatically initialised.

    pixar_sr: number (optional)
        The nominal pixel area for the detector in steradians.
        If provided, the value is written to the PIXAR_SR keyword.
    pixel_a2: number (optional)
        The nominal pixel area for the detector in square arcseconds
        If provided, the value is written to the PIXAR_A2 keyword.
    
    """
    # Both models use exactly the same schema.
    schema_url = "miri_photom.schema"
    fieldnames = ('filter', 'subarray', 'photmjsr', 'uncertainty', 'nelem',
                  'wavelength', 'relresponse', 'relresperror')

    def __init__(self, init=None, subarray=None, srf_table=None,
                 pixar_sr=None, pixar_a2=None, **kwargs):
        """
        
        Initialises the MiriLrsPhotometricModel class.
        
        Parameters: See class doc string.

        """
        if srf_table is not None:
            # Construct a phot_table from the given subarray and srf_table and
            # some dummy photmjsr and uncertainty arrays.
            # TODO: A slow, brute force method. There is probably a faster way,
            # but the function is not time critical for small CDPs..
            mirifilter = 'P750L'
            if subarray is None or not subarray:
                strg = "If an srf_table is parameter provided, "
                strg += "a subarray parameter must also be provided."
                raise AttributeError(strg)
            photmjsr = 1.0
            uncertainty = 0.0
            wavelength = [0.0] * MAX_NELEM
            relresponse = [0.0] * MAX_NELEM
            relresperror = [0.0] * MAX_NELEM
            ii = 0
            for (wav, srf, unc) in srf_table:
                wavelength[ii] = wav
                relresponse[ii] = srf
                relresperror[ii] = unc
                ii += 1
            nelem = len(srf_table)
            new_phot_table = [(mirifilter, subarray, photmjsr, uncertainty,
                              nelem, tuple(wavelength), tuple(relresponse),
                              tuple(relresperror))]
        else:
            new_phot_table = None

        # Pass the phot_table to the generic model.
        super(MiriLrsPhotometricModel, self).__init__(init=init,
                                                      phot_table=new_phot_table,
                                                      pixar_sr=pixar_sr,
                                                      pixar_a2=pixar_a2,
                                                      **kwargs)
        #self.add_comment("RELRESPONSE is absolute response so PHOTMJSR is 1.0.")
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type:
            self.meta.model_type = model_type
            
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
    schema_url = "miri_pixelarea.schema"

    def __init__(self, init=None, data=None, pixar_sr=None, pixar_a2=None,
                 **kwargs):
        """
        
        Initialises the MiriPixelAreaModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPixelAreaModel, self).__init__(init=init, **kwargs)

        # Data type is AREA.
        self.meta.reftype = 'AREA'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        
        # This is a reference data model.
        self._reference_model()

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
    wavelength = [0.0] * MAX_NELEM
    relresponse = [0.0] * MAX_NELEM
    relresperror = [0.0] * MAX_NELEM
    # There is a separate response record for each filter, valid for all
    # subarrays.
    phot_im1 = [('F560W',  'GENERIC', 2.41,  0.26,  0, wavelength, relresponse, relresperror),
               ('F770W',   'GENERIC', 1.32,  0.013, 0, wavelength, relresponse, relresperror),
               ('F1000W',  'GENERIC', 1.76,  0.12,  0, wavelength, relresponse, relresperror),
               ('F1130W',  'GENERIC', 5.76,  0.43,  0, wavelength, relresponse, relresperror),
               ('F1280W',  'GENERIC', 2.11,  0.16,  0, wavelength, relresponse, relresperror),
               ('F1500W',  'GENERIC', 1.84,  0.01,  0, wavelength, relresponse, relresperror),
               ('F1800W',  'GENERIC', 2.68,  0.23,  0, wavelength, relresponse, relresperror),
               ('F2100W',  'GENERIC', 2.04,  0.15,  0, wavelength, relresponse, relresperror),
               ('F2550W',  'GENERIC', 4.25,  0.4,   0, wavelength, relresponse, relresperror),
               ('F2550WR', 'GENERIC', 4.60,  0.24,  0, wavelength, relresponse, relresperror),
               ('F1065C',  'GENERIC', 1.37,  0.1,   0, wavelength, relresponse, relresperror),
               ('F1140C',  'GENERIC', 1.43,  0.11,  0, wavelength, relresponse, relresperror),
               ('F1550C',  'GENERIC', 1.81,  0.13,  0, wavelength, relresponse, relresperror),
               ('F2300C',  'GENERIC', 3.65,  0.23,  0, wavelength, relresponse, relresperror),
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
                                           useafter='TODAY',
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
                                           useafter='TODAY',
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
    relresperror = []
    resp = 0
    nelm = 0
    for ii in range(300):
        wav = float(ii)/12.0
        resp = (resp + 1) % 30
        r10 = 0.1 + resp/35.0
        wavelength.append(wav)
        relresponse.append(r10)
        relresperror.append(0.001)
        nelm += 1
    for ii in range(300,MAX_NELEM):
        # Pad unused elements with zero.
        wavelength.append(0.0)
        relresponse.append(0.0)
        relresperror.append(0.0)
    # The LRS data is for the P750L filter only, but there is separate
    # record for the SLITLESSPRISM subarray.
    phot_lrs = [('P750L',  'FULL',          1.0,  0.0,  nelm, wavelength, relresponse, relresperror),
                ('P750L',  'SLITLESSPRISM', 0.9,  0.0,  nelm, wavelength, relresponse, relresperror)
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
                                           useafter='TODAY',
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
                                           useafter='TODAY',
                                           description='Test data')
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_pixelarea_model1.fits", overwrite=True)
        del testdata
       
    print("Test finished.")
