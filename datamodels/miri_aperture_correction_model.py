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

# List all classes and global functions here.
__all__ = ['MiriMrsApertureCorrectionModel']


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
    schema_url = "miri_aperture_correction_mrs.schema.yaml"
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
        if model_type is not None:
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
        
    print("Test finished.")
