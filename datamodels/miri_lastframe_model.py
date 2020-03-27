#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI last frame 
correction data. Essentially the same as the MIRI measured image data model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

27 Mar 2014: Created
09 Jul 2014: field_def changed to dq_def.
16 Jul 2014: Removed unnecessary arguments from __init__.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
07 Nov 2014: Modified to use miri_lastframe.schema, which defines
             the expected data and err units.
20 Aug 2015: Duplicated parts of schema now reference STScI model.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
26 Mar 2020: Ensure the model_type remains as originally defined when saving
             to a file.

@author: Vincent Geers (DIAS)

"""

import numpy as np

# Import the MIRI image model.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['MiriLastFrameModel']

lastframe_reference_setup = \
            [(0, 'DO_NOT_USE',         'Bad pixel. Do not use.'),
             (1, 'NO_LASTFRAME_CORR',  'No last frame correction')]
#              (2, 'CDP_LOW_QUAL',       'Data of low quality')]
lastframe_reference_flags = insert_value_column( lastframe_reference_setup )


class MiriLastFrameModel(MiriMeasuredModel):
    """
    
    A data model for MIRI last frame correction data.
        
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
        An array containing the last frame data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    err: numpy array (optional)
        An array containing the error data.
        Must be broadcastable onto the data array.
    dq: numpy array (optional)
        An array containing the quality data.
        Must be broadcastable onto the data array.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, comment:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    # TODO: Could the lastframe and reset models share the same schema?
    schema_url = "miri_lastframe.schema"
    _default_dq_def = lastframe_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None,
                 dq_def=None, **kwargs):
        """
        
        Initialises the MiriLastFrameModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLastFrameModel, self).__init__(init=init, data=data,
                                                 dq=dq, err=err,
                                                 dq_def=dq_def, **kwargs)

        # Data type is last frame.
        self.meta.reftype = 'LASTFRAME'
        # Initialise the model type
        self._init_data_type()      
        # This is a reference data model.
        self._reference_model()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriLastFrameModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriLastFrameModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])

    print("Last frame data with data + err + dq:")
    with MiriLastFrameModel(data=data3x3, err=err3x3, dq=dq3x3,
                            dq_def=lastframe_reference_flags ) as testdata:
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_lastframe_model1.fits", overwrite=True)
        del testdata

    print("Test finished.")
