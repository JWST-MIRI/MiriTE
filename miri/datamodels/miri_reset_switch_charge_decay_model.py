#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An experimental extension to the standard STScI data model, which 
defines MIRI reset switch charge decay model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

19 Nov 2015: Created as a skeleton. Three pieces of information needed:
             (1) Names of decay constants; (2) units of decay constants;
             and (3) any additional metadata.
27 Nov 2015: Updated to match actual data file delivered for CDP-5.
01 Jun 2016: "ROWS" column added to RSCD table.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
15 Nov 2018: Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
04 Oct 2019: Updated to match build 7.3 data model.
26 Mar 2020: Ensure the model_type remains as originally defined when saving
             to a file.
11 May 2020: Removed the CDP-6 version of the data model.
24 Jul 2020: Remove the dependency on the actual structure of the data model,
             since it changes frequently.

@author: Steven Beard (UKATC)

"""

#import warnings
#import numpy as np

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriResetSwitchChargeDecayModel']


class MiriResetSwitchChargeDecayModel(MiriDataModel):
    """
    
    A generic data model for a MIRI RSCD table, based on the STScI
    base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the flux parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.

    """
    schema_url = "miri_reset_switch_charge_decay.schema"
    fieldnames_group = ('subarray', 'readpatt', 'group_skip')
    fieldnames_gen = ('subarray', 'readpatt', 'lower_cuttoff', 'alpha_even', 'alpha_odd')
   
    def __init__(self, init=None, **kwargs):
        """
        
        Initialises the MiriResetSwitchChargeDecayModel class.
        
        Parameters: See class doc string.

        NOTE: This data model changes sufficiently frequently that the MIRI
        data model does not attempt to process it.

        """
        super(MiriResetSwitchChargeDecayModel, self).__init__(init=init, **kwargs)

        # Data type is RSCD.
        self.meta.reftype = 'RSCD'
        # Initialise the model type
        self._init_data_type()       
        # This is a reference data model.
        self._reference_model()
        
    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriResetSwitchChargeDecayModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriResetSwitchChargeDecayModel module.")
    
    PLOTTING = False
    SAVE_FILES = False

    # Test an empty data model.
    with MiriResetSwitchChargeDecayModel( ) as testrscd1:
        testrscd1.set_instrument_metadata(detector='MIRIFUSHORT',
                                          ccc_pos='OPEN',
                                          deck_temperature=11.0,
                                          detector_temperature=6.0)
        testrscd1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', useafter='2020-05-20',
                                           description='Test data')
        print(testrscd1)
        if PLOTTING:
            testrscd1.plot(description="testrscd1")
        if SAVE_FILES:
            testrscd1.save("test_rscd_model1.fits", overwrite=True)
        del testrscd1
        
    print("Test finished.")
