#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI pixel saturation 
data. Essentially the same as the MIRI measured image data model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

26 Feb 2013: Created
04 Mar 2013: Median saturation keyword, SATMED, added to metadata.
04 Jun 2013: Shortened the names of the ramp, slope and image models.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
09 Jul 2014: JWST reference flags table added.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
19 Aug 2015: Use MiriMeasuredModel, not MiriImageModel.
20 Aug 2015: Duplicated parts of schema now reference STScI model.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Steven Beard (UKATC)

"""

import numpy as np
#import numpy.ma as ma

# Import the MIRI image model.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['saturation_reference_flags', 'MiriPixelSaturationModel']

# The new JWST saturation reference flags
saturation_reference_setup = \
            [(0, 'DO_NOT_USE',   'Bad pixel. Do not use.'),
             (1, 'NO_SAT_CHECK', 'Saturation check not available')]
saturation_reference_flags = insert_value_column( saturation_reference_setup )


class MiriPixelSaturationModel(MiriMeasuredModel):
    """
    
    A data model for MIRI pixel saturation data.
        
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
        An array containing the pixel saturation data.
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
    tolerance: number, optional
        Pixel saturation tolerance for a 2pt difference. If specified,
        this value is written to the SATT2PT keyword.
    lowlimit: number, optional
        Lower limit on pixel saturation. If specified, this value is
        written to the SATL keyword.
    median: number, optional
        Median pixel saturation. If specified, this value is
        written to the SATMED keyword.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    # There is no separate pixel saturation schema as yet.
    schema_url = "miri_pixel_saturation.schema"
    _default_dq_def = saturation_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None,
                 dq_def=None, tolerance=None, lowlimit=None, median=None,
                 **kwargs):
        """
        
        Initialises the MiriPixelSaturationModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPixelSaturationModel, self).__init__(init=init, data=data,
                                                 dq=dq, err=err,
                                                 dq_def=dq_def, **kwargs)

        # Data type is pixel saturation.
        self.meta.reftype = 'SATURATION'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()
        
        # Set the metadata, if provided.
        if tolerance is not None:
            self.meta.sat_tolerance = tolerance
        if lowlimit is not None:
            self.meta.sat_lower = lowlimit
        if median is not None:
            self.meta.sat_median = median
#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriPixelSaturationModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])

    print("Pixel saturation data with data + err + dq:")
    with MiriPixelSaturationModel(data=data3x3, err=err3x3, dq=dq3x3,
                            dq_def=saturation_reference_flags, tolerance=3,
                            lowlimit=65535, median=59047.0) as testdata:
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_saturation_model1.fits", overwrite=True)
        del testdata

    print("Test finished.")
