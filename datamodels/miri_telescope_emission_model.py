#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI telescope 
emission data. Essentially the same as the MIRI measured image data 
model, with additional metadata.

NOTE: This is the simple way of describing MIRI telescope emission data. 
There is a separate data object for each filter and temperature. An 
alternative model would be to have one large table consisting of the 
following columns:

   filter, temperature, emission-map

This would be a more complex model.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

10 Jan 2013: Created
21 Jan 2013: Included data filter and telescope temperature in metadata.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
26 Feb 2013: Only update the filter and temperature if explicitly
             provided.
04 Jun 2013: Shortened the names of the ramp, slope and image models.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
09 Jul 2014: field_def changed to dq_def.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
19 Aug 2015: Removed MiriImageModel.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jul 2016: Removed obsolete flags.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import numpy as np
#import numpy.ma as ma

# Import the MIRI measured data model.
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['MiriTelescopeEmissionModel']

telem_reference_setup = \
            [(0, 'DO_NOT_USE',         'Bad pixel. Do not use.'),
             (1, 'UNRELIABLE_SLOPE',   'Slope variance large')]
telem_reference_flags = insert_value_column( telem_reference_setup )


class MiriTelescopeEmissionModel(MiriMeasuredModel):
    """
    
    A data model for MIRI telescope emission data.
        
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the data parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    data: numpy array (optional)
        An array containing the telescope emission data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    err: numpy array (optional)
        An array containing the error data.
        Must be broadcastable onto the data array.
    dq: numpy array (optional)
        An array containing the quality data.
        Must be broadcastable onto the data array.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    filt: string, optional
        The MIRI filter for which these telescope emission data are valid
    temperature: float, optional
        The temperature for which these telescope emission data are valid.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_telescope_emission.schema.yaml"
    _default_dq_def = telem_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 filt=None, temperature=None, **kwargs):
        """
        
        Initialises the MiriTelescopeEmissionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriTelescopeEmissionModel, self).__init__(init=init, data=data,
                                                    dq=dq, err=err,
                                                    dq_def=dq_def,
                                                    **kwargs)
        
        # Data type is telescope emission map.
        self.meta.model_type = 'TEL_EMISSION'
        self.meta.reftype = 'TEL_EMISSION'
        
        # This is a reference data model.
        self._reference_model()

        # Add filter and temperature to the metadata
        if filt is not None:
            self.meta.instrument.filter = filt
        if temperature is not None:
            self.meta.telescope_temperature = temperature

    def __str__(self):
        """
        
        Return the contents of the telescope emission object as a readable
        string.
        
        """
        # First obtain a string describing the underlying data model.
        strg = super(MiriTelescopeEmissionModel, self).__str__()
        
        # Add the extras
        if self.meta.instrument.filter is not None:
            strg += "Data valid for filter=\'%s\' " % \
                self.meta.instrument.filter
        else:
            strg += "Data valid for UNKNOWN filter "
        if self.meta.telescope_temperature is not None:
            strg += "and telescope temperature=%.2fK" % \
                self.meta.telescope_temperature
        else:
            strg += "and UNKNOWN telescope temperature"
        return strg

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriTelescopeEmissionModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])

    print("Telescope emission data with data + err + dq:")
    with MiriTelescopeEmissionModel(data=data3x3, err=err3x3, dq=dq3x3,
                                    dq_def=telem_reference_flags, filt='F560W',
                                    temperature=14.6) \
            as testdata:
        print(testdata)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_emission_model1.fits", overwrite=True)
        del testdata

    print("Test finished.")
