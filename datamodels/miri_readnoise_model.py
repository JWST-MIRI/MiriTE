#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI readnoise data, based
on the base MIRI data model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

16 Jul 2014: Created
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
17 Oct 2014: Updated MiriGainModel which is now based on MiriDataModel, and 
             no longer has an ERR nor DQ extension.
20 May 2015: Corrected mistake in which data array was ignored. The model
             now inherits from the HasData base class.
19 Aug 2015: Duplicated parts of schema now reference STScI model.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Vincent Geers (DIAS), Steven Beard (UKATC)

"""

import numpy as np

# Import the MIRI data model.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasData

# List all classes and global functions here.
__all__ = ['MiriReadnoiseModel']


class MiriReadnoiseModel(MiriDataModel, HasData):
    """
    
    A data model for MIRI readnoise data.
        
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
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_readnoise.schema.yaml"

    def __init__(self, init=None, data=None, **kwargs):
        """
        
        Initialises the MiriReadnoiseModel class.
        
        Parameters: See class doc string.

        """
        super(MiriReadnoiseModel, self).__init__(init=init, **kwargs)

        # Data type is last frame.
        self.meta.reftype = 'READNOISE'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()

        # Update the data array if it has been specifically provided.
        HasData.__init__(self, data)

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriReadnoiseModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])

    print("Readnoise data with data + err + dq:")
    with MiriReadnoiseModel(data=data3x3 ) as testdata:
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_readnoise_model1.fits", overwrite=True)
        del testdata

    print("Test finished.")
