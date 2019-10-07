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

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

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
           'MiriLrsPositionCorrectionModel']



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
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type
        
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


class MiriLrsThroughputCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS throughput correction table,
    based on the STScI base model, DataModel.
    The throughput is calculated as the flux loss of the centered PSF due to the
    finite extension of the slit.

    See " " for a detailed description of the data model.
    
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
    schema_url = "miri_throughput_correction_lrs.schema.yaml"
    fieldnames = ('wavelength', 'through_corr', 'through_corr_err')
    
    def __init__(self, init=None, throughcorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsThroughputCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsThroughputCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is throughput correction.
        self.meta.reftype = 'THROUGHCORR'
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type
        
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

class MiriLrsApertureCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS aperture correction table,
    based on the STScI base model, DataModel.
    The aperture correction is calculated as the flux loss ratio of the 
    centered PSF in a user given aperture
    to the extension of the complete slit.

    See " " for a detailed description of the data model.
    
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
    schema_url = "miri_aperture_correction_lrs.schema.yaml"
    fieldnames = ('wavelength', 'aper_corr')
    
    def __init__(self, init=None, apercorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsApertureCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsApertureCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is aperture correction.
        self.meta.reftype = 'APERCORR'
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type
        
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
        
class MiriLrsPositionCorrectionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS position correction table,
    based on the STScI base model, DataModel.
    It calculates the flux loss for offsets of the webbPsf in dispersion direction from the slit center in steps 
    of a quarter of a pixel.

    See " " for a detailed description of the data model.
    
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
    schema_url = "miri_position_correction_lrs.schema.yaml"
    fieldnames = ('wavelength', 'pos_corr')
    
    def __init__(self, init=None, poscorr_table=None, **kwargs):
        """
        
        Initialises the MiriLrsPositionCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsPositionCorrectionModel, self).__init__(init=init, **kwargs)

        # Data type is position correction.
        self.meta.reftype = 'POS_CORR'
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type
        
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


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsApertureCorrectionModel module.")
    
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
   
    apercorr = np.zeros([40])
    apercorr = np.linspace(0.5,1,num = 40) 
    
    apercorr_table=[]
    apercorr_table.append((wave[0], apercorr.tolist()))
    apercorr_table.append((wave[1], apercorr.tolist()))
    apercorr_table.append((wave[2], apercorr.tolist()))
    
    with MiriLrsApertureCorrectionModel( apercorr_table=apercorr_table) as testapercorr_lrs:
        testapercorr_lrs.set_instrument_metadata(detector='MIRIMAGE', filt='P750L')
        testapercorr_lrs.set_subarray_metadata('FULL')
        testapercorr_lrs.set_housekeeping_metadata('MPIA', author='Juergen Schreiber',
                                           version='1.0', useafter='2019-07-19',
                                           description='Test data')
        print(testapercorr_lrs)
        if PLOTTING:
            testapercorr_lrs.plot(description="testapercorr_lrs")
        if SAVE_FILES:
            testapercorr_lrs.save("test_apercorr_model_lrs.fits", overwrite=True)
        del apercorr_table
        del testapercorr_lrs
    
    print("Test finished.")
