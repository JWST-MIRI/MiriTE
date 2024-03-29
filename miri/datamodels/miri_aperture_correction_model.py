#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines the MIRI
aperture correction model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

20 Nov 2015: Created (using the transmission correction table as a template).
10 Dec 2015: APPCORR changed to APERCORR throughout.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema. Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
26 Mar 2020: Ensure the model_type remains as originally defined when saving
             to a file.
22 Nov 2021: MiriLrsPathlossCorrectionModel added.
             
@author: Steven Beard (UKATC), Vincent Geers (UKATC), Juergen Schreiber (MPIA)

"""

import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel
import numpy as np

# List all classes and global functions here.
__all__ = ['MiriMrsApertureCorrectionModel',\
           'MiriLrsApertureCorrectionModel',\
           'MiriLrsThroughputCorrectionModel',\
           'MiriLrsPositionCorrectionModel',\
           'MiriLrsPathlossCorrectionModel']



class MiriMrsApertureCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI aperture correction table,
    based on the STScI base model, DataModel.

    See MIRI-DS-00011-NLC for a detailed description of the content
    of the data model.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    apercorr_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the aperture
        correction table;
        Or: A numpy record array containing the same information as above.
        A aperture correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_aperture_correction_mrs.schema"
    fieldnames = ('wavelength', 'solution', 'a_aperture',
                  'a_annulus_in', 'a_annulus_out',
                  'a_over_b', 'pos_angle', 'aper_corr', 'aper_corr_err')
    
    def __init__(self, init=None, apercorr_table=None, **kwargs):
        """
        
        Initialises the MiriMrsApertureCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsApertureCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is aperture correction.
        self.meta.reftype = 'APERCORR'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()
        
        if apercorr_table is not None:
            try:
                self.apercorr_table = apercorr_table
            except (ValueError, TypeError) as e:
                strg = "apercorr_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
            
        # Copy the table column units from the schema, if defined.
        apercorr_units = self.set_table_units('apercorr_table')
 
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriMrsApertureCorrectionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()


class MiriLrsThroughputCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS throughput correction table,
    based on the STScI base model, DataModel.
    The throughput is calculated as the flux loss of the centered PSF due to the
    finite extension of the slit.

    See "MIRI-TN-10023-MPI" for a detailed description of the data model.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    throughcorr_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the throughput
        correction table;
        Or: A numpy record array containing the same information as above.
        A throughput correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_throughput_correction_lrs.schema"
    fieldnames = ('wavelength', 'through_corr', 'through_corr_err')
    
    def __init__(self, init=None, throughcorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsThroughputCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsThroughputCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is throughput correction.
        self.meta.reftype = 'THROUGHCORR'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()
        
        if throughcorr_table is not None:
            try:
                self.throughcorr_table = throughcorr_table
            except (ValueError, TypeError) as e:
                strg = "throughcorr_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
            
        # Copy the table column units from the schema, if defined.
        throughcorr_units = self.set_table_units('throughcorr_table')
        
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriLrsThroughputCorrectionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

class MiriLrsApertureCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS aperture correction table,
    based on the STScI base model, DataModel.
    The aperture correction is calculated as the flux loss ratio of the 
    centered PSF in a user given aperture
    to the extension of the complete slit.

    See "MIRI-TN-10021-MPI" for a detailed description of the data model.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    apcorr_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the aperture
        correction table;
        Or: A numpy record array containing the same information as above.
        A aperture correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
        The table contains the columns wavelength, nelem (number of apertures),
        width (width of aperture in pixels) and apcorr which is the correction
        factor for each width to an infinite aperture
        
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_aperture_correction_lrs.schema"
    fieldnames = ('subarray','wavelength', 'nelem_wl', 'size', 'nelem_size', 'apcorr', 'apcorr_err')
    
    def __init__(self, init=None, apcorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsApertureCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsApertureCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is aperture correction.
        self.meta.reftype = 'APCORR'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()
        
        if apcorr_table is not None:
            try:
                self.apcorr_table = apcorr_table
            except (ValueError, TypeError) as e:
                strg = "apcorr_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
            
        # Copy the table column units from the schema, if defined.
        apcorr_units = self.set_table_units('apcorr_table')
        
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriLrsApertureCorrectionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()


class MiriLrsPositionCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS position correction table,
    based on the STScI base model, DataModel.
    It calculates the flux loss for offsets of the webbPsf in dispersion direction from the slit center in steps 
    of a quarter of a pixel.

    See "MIRI-TN-10024-MPI" for a detailed description of the data model.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    poscorr_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the position
        correction table;
        Or: A numpy record array containing the same information as above.
        A position correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_position_correction_lrs.schema"
    fieldnames = ('wavelength', 'pos_corr')
    
    def __init__(self, init=None, poscorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsPositionCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsPositionCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is position correction.
        self.meta.reftype = 'POSCORR'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()
        
        if poscorr_table is not None:
            try:
                self.poscorr_table = poscorr_table
            except (ValueError, TypeError) as e:
                strg = "poscorr_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Copy the table column units from the schema, if defined.
        poscorr_units = self.set_table_units('poscorr_table')        
        
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriLrsPositionCorrectionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

class MiriLrsPathlossCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS pathloss correction table,
    based on the STScI base model, DataModel.
    It calculates the flux loss for offsets of the webbPsf from the slit center in all directions.

    See "MIRI-TN-10025-MPI" for a detailed description of the data model.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape.
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    pathloss_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing columns in the position
        correction table;
        Or: A numpy record array containing the same information as above.
        A position correction table must either be defined in the
        initializer or in this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_pathloss_correction_lrs.schema"
    fieldnames = ('wavelength', 'pathloss', 'pathloss_err')
    
    def __init__(self, init=None, pathloss_table=None, **kwargs):
        """
        
        Initialises the MiriLrsPathlossCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsPathlossCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is position correction.
        self.meta.reftype = 'PATHLOSS'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()
        
        if pathloss_table is not None:
            try:
                self.pathloss_table = pathloss_table
            except (ValueError, TypeError) as e:
                strg = "pathloss_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Copy the table column units from the schema, if defined.
        poscorr_units = self.set_table_units('pathloss_table')        
        
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriLrsPathlossCorrectionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the miri_aperture_correction_model module.")
    
    PLOTTING = False
    SAVE_FILES = False

    apercorrdata = [(5.0, 'nominal', 4.87,  7.76,  8.57, 1.0,  0.0, 42.0, 0.1),
                   (10.0, 'nominal', 5.87,  8.76,  9.57, 1.0,  0.0, 32.0, 0.1),
                   (15.0, 'nominal', 6.87,  9.76, 10.57, 1.0,  0.0, 32.0, 0.1),
                   (15.0, 'nominal', 7.87, 10.76, 11.57, 0.7, 45.0, 12.0, 0.1)]

    print("\nAperture correction with factors derived from list of tuples:")
    with MiriMrsApertureCorrectionModel( apercorr_table=apercorrdata ) as testapercorr1:
        testapercorr1.set_instrument_metadata(detector='MIRIFUSHORT',
                                         channel='1', band='SHORT',
                                         ccc_pos='OPEN')
        testapercorr1.set_subarray_metadata('FULL')
        testapercorr1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2015-11-20',
                                           description='Test data')
        print(testapercorr1)
        if PLOTTING:
            testapercorr1.plot(description="testapercorr1")
        if SAVE_FILES:
            testapercorr1.save("test_apercorr_model1.fits", overwrite=True)
        del testapercorr1
    
    throughcorrdata_lrs = [( 4.0, 0.9,  0.01),
                        ( 6.0, 0.85,  0.01),
                        ( 8.0, 0.8,  0.01),
                        (10.0, 0.75,  0.01),
                        (12.0, 0.7,  0.01),
                        (14.0, 0.65,  0.01),
                        (16.0, 0.6,  0.01)]

    print("\nThroughput correction with factors derived from list of tuples:")
    with MiriLrsThroughputCorrectionModel( throughcorr_table=throughcorrdata_lrs ) as testthroughcorr_lrs:
        testthroughcorr_lrs.set_instrument_metadata(detector='MIRIMAGE', filt='P750L')
        testthroughcorr_lrs.set_subarray_metadata('FULL')
        testthroughcorr_lrs.set_housekeeping_metadata('MPIA', author='Juergen Schreiber',
                                           version='1.0', useafter='2019-07-19',
                                           description='Test data')
        print(testthroughcorr_lrs)
        if PLOTTING:
            testthroughcorr_lrs.plot(description="testthroughcorr_lrs")
        if SAVE_FILES:
            testthroughcorr_lrs.save("test_throughcorr_model_lrs.fits", overwrite=True)
        del testthroughcorr_lrs
    
    
    wave = [5.,8.,11.]
    poscorr = np.zeros([81])
    poscorr = np.append(np.linspace(0.5,1,num = 41), np.linspace(0.99,0.5, num = 40))
     
    
    poscorr_table=[]
    poscorr_table.append((wave[0], poscorr.tolist()))
    poscorr_table.append((wave[1], poscorr.tolist()))
    poscorr_table.append((wave[2], poscorr.tolist()))
    
    with MiriLrsPositionCorrectionModel( poscorr_table=poscorr_table) as testposcorr_lrs:
        testposcorr_lrs.set_instrument_metadata(detector='MIRIMAGE', filt='P750L')
        testposcorr_lrs.set_subarray_metadata('FULL')
        testposcorr_lrs.set_housekeeping_metadata('MPIA', author='Juergen Schreiber',
                                           version='1.0', useafter='2019-07-19',
                                           description='Test data')
        print(testposcorr_lrs)
        if PLOTTING:
            testposcorr_lrs.plot(description="testposcorr_lrs")
        if SAVE_FILES:
            testposcorr_lrs.save("test_poscorr_model_lrs.fits", overwrite=True)
        del poscorr_table
        del testposcorr_lrs
   
    apcorr = np.zeros([388,40])
    apcorr_x = np.linspace(1.5,1,num = 40)
    apcorr[:,:] = apcorr_x
    apcorr_err = np.zeros([388, 40]) 
    width = np.arange(40) + 1
    wav = np.arange(388)/100.
    print(width.tolist())
    print(apcorr.tolist())
    apcorr_table=[]
    apcorr_table.append(('FULL' ,wav.tolist(), 388, width.tolist(), 40, apcorr.tolist(), apcorr_err.tolist()))
    apcorr_table.append(('SLITLESSPRISM' ,wav.tolist(), 388, width.tolist(), 40, apcorr.tolist(), apcorr_err.tolist()))
    print(apcorr_table)
    
    with MiriLrsApertureCorrectionModel( apcorr_table=apcorr_table) as testapcorr_lrs:
        testapcorr_lrs.set_instrument_metadata(detector='MIRIMAGE', filt='P750L')
        testapcorr_lrs.set_housekeeping_metadata('MPIA', author='Juergen Schreiber',
                                           version='1.0', useafter='2019-07-19',
                                           description='Test data')
        print(testapcorr_lrs)
        if PLOTTING:
            testapcorr_lrs.plot(description="testapcorr_lrs")
        if SAVE_FILES:
            testapcorr_lrs.save("test_apcorr_model_lrs.fits", overwrite=True)
        del apcorr_table
        del testapcorr_lrs
    
    
    
    pathloss = np.zeros([3,40,4])
    pathloss_err = np.zeros([3,40,4])
    pathloss_table=[]
    pathloss_table.append((wave[0], pathloss[0].tolist(), pathloss_err[0].tolist()))
    pathloss_table.append((wave[1], pathloss[1].tolist(), pathloss_err[1].tolist()))
    pathloss_table.append((wave[2], pathloss[2].tolist(), pathloss_err[2].tolist()))  
    print(pathloss_table)
    
    with MiriLrsPathlossCorrectionModel( pathloss_table=pathloss_table) as testpathloss_lrs:
        testpathloss_lrs.set_instrument_metadata(detector='MIRIMAGE', filt='P750L')
        testpathloss_lrs.set_housekeeping_metadata('MPIA', author='Juergen Schreiber',
                                           version='1.0', useafter='2019-07-19',
                                           description='Test data')
        testpathloss_lrs.set_wcs_metadata(wcsaxes = 2, crpix = [1,1], cdelt = [1,1], crval = [-1, -10])
        print(testpathloss_lrs)
        if PLOTTING:
            testpathloss_lrs.plot(description="testpathloss_lrs")
        if SAVE_FILES:
            testpathloss_lrs.save("test_pathloss_model_lrs.fits", overwrite=True)
        del pathloss_table
        del testpathloss_lrs
          
    print("Test finished.")
