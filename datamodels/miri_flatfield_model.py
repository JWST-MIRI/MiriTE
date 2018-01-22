#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI pixel flat-field 
data. Essentially the same as the MIRI measured data model.

:Reference:

The STScI jwst.datamodels documentation. See

http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

10 Jan 2013: Created
21 Jan 2013: Included data type in metadata.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
19 Feb 2013: Changed MiriImageModel to MiriMeasuredModel to allow
             flat field cubes (Vincent Geers, DIAS).
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
09 Jul 2014: JWST reference flags table added.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
13 Oct 2014: Re-instated SKYFLAT as valid type, for MIRI Imager.
09 Dec 2014: Metadata added to example data product.
11 Mar 2015: group_integration_time changed to group_time.
20 Aug 2015: Duplicated parts of schema now reference STScI model.
20 Nov 2015: Resurrected the PIXELFLAT option.
16 Feb 2016: Default pedigree for SKYFLAT changed from 'GROUND' to 'DUMMY'.
23 Mar 2016: Set default fill value to 1.0.
08 Jun 2016: Corrected typo in set_exposure_type. Added "detector" to
             constructor arguments.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
             Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Nov 2017: Added more DQ flags.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import numpy as np
#import numpy.ma as ma

# Import the MIRI measured data model.
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel, HasDataErrAndDq

# List all classes and global functions here.
__all__ = ['flat_reference_flags', 'MiriFlatfieldModel']

# The new JWST flat-field reference flags
flat_reference_setup = \
            [(0, 'DO_NOT_USE',          'Bad pixel. Do not use.'),
             (1, 'NON_SCIENCE',         'Pixel not on science portion of detector'),
             (2, 'UNRELIABLE_FLAT',     'Flat variance large'),
             (3, 'CDP_PARTIAL_DATA',    'Data derived from incomplete input'),
             (4, 'CDP_LOW_QUAL',        'Data of low quality'),
             (5, 'CDP_UNRELIABLE_ERROR','Data without reliable error estimate'),
             (6, 'NO_FLAT_FIELD',       'No flat-field data available'),
             (7, 'DIFF_PATTERN',        'Diffraction pattern')]
flat_reference_flags = insert_value_column( flat_reference_setup )


class MiriFlatfieldModel(MiriMeasuredModel):
    """
    
    A data model for MIRI flat-field data.
        
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
        An array containing the flat-field data.
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
    flattype: str (optional)
        A string giving the specific kind of flat-field contained
        in the data object: 'SkyFlat', 'PixelFlat' or 'FringeFlat'.
        If not given, the type defaults to the generic term 'FLAT'.
    detector: str (optional)
        The name of the detector associated with this fiat-field data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    # There is no separate flat-field schema as yet.
    schema_url = "miri_flatfield.schema.yaml"
    _default_dq_def = flat_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None, fiterr=None,
                 dq_def=None, flattype='', detector=None, **kwargs):
        """
        
        Initialises the MiriFlatfieldModel class.
        
        Parameters: See class doc string.

        """
        super(MiriFlatfieldModel, self).__init__(init=init, data=data,
                                                 dq=dq, err=err,
                                                 dq_def=dq_def, **kwargs)
        # Missing sections of the flat-field should be filled with 1.0
        HasDataErrAndDq.set_data_fill( self, 1.0 )

        # Data type is flat-field.
        if not flattype:
            # Set to the default type, if not already defined
            if not self.meta.reftype:
                self.meta.model_type = 'FLAT'
                self.meta.reftype = 'FLAT'
        else:
            self.meta.model_type = flattype
            ftupper = flattype.upper()
            if "FRINGE" in ftupper:
                self.meta.reftype = "FRINGE"
            elif "SKY" in ftupper:
                self.meta.reftype = "SKYFLAT"
            if "PIX" in ftupper:
                self.meta.reftype = "PIXELFLAT"
            else:
                # Pixel flat is just FLAT. 
                self.meta.reftype = "FLAT"
       
        # The default pedigree is 'DUMMY' for a sky flat and 'GROUND'
        # for everything else.
        if not self.meta.pedigree:
            if "SKY" in self.meta.reftype:
                self.meta.pedigree = 'DUMMY'
            else:
                self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'

        # Define the detector identifier, if specified.
        if detector is not None:
            self.meta.instrument.detector = detector
 
        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type( datatype='FLAT' )
                
        # The fill value for a flat-field is 1.0
        self._data_fill = 1.0
        self._data_fill_value = 1.0


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriFlatfieldModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.0,1.2,1.1],[1.3,1.2,1.0],[1.1,0.8,0.9]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])

    print("\nFlat-field data with data + err + dq:")
    with MiriFlatfieldModel(data=data3x3, err=err3x3, dq=dq3x3,
                            dq_def=flat_reference_flags) as testdata:
        # Add some example metadata.
        testdata.set_instrument_metadata(detector='MIRIFUSHORT',
                                         channel='1',
                                         ccc_pos='OPEN',
                                         deck_temperature=11.0,
                                         detector_temperature=6.0)
        testdata.set_exposure_metadata(readpatt='FAST',
                                       nints=1, ngroups=1,
                                       frame_time=1.0,
                                       integration_time=10.0,
                                       group_time=10.0,
                                       reset_time=0, frame_resets=3)
        testdata.set_subarray_metadata('FULL')
        testdata.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testdata)
        print("Filled data:\n", testdata.data_filled)
        if PLOTTING:
            testdata.plot(description="testdata")
        if SAVE_FILES:
            testdata.save("test_flatfield_model1.fits", overwrite=True)
        del testdata

    print("\nPIXEL Flat-field data with data + err + dq:")
    with MiriFlatfieldModel(data=data3x3, err=err3x3, dq=dq3x3,
                            dq_def=flat_reference_flags, flattype='PIXELFLAT',
                            detector='MIRIFUSHORT') as testdata2:
        # Add some example metadata.
        testdata2.set_instrument_metadata(detector='MIRIFUSHORT',
                                         channel='1',
                                         ccc_pos='OPEN',
                                         deck_temperature=11.0,
                                         detector_temperature=6.0)
        testdata2.set_exposure_metadata(readpatt='FAST',
                                       nints=1, ngroups=1,
                                       frame_time=1.0,
                                       integration_time=10.0,
                                       group_time=10.0,
                                       reset_time=0, frame_resets=3)
        testdata2.set_subarray_metadata('FULL')
        testdata2.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='Test data')
        print(testdata2)
        print("Filled data:\n", testdata2.data_filled)
        if PLOTTING:
            testdata2.plot(description="testdata2")
        if SAVE_FILES:
            testdata2.save("test_flatfield_model2.fits", overwrite=True)
        del testdata2

    print("Test finished.")
