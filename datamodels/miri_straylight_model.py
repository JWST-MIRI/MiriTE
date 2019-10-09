#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines a MIRI Straylight 
pixel mask model for MRS.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

11 Mar 2014: Created
19 Mar 2014: Minor update to tests in main program.
10 Apr 2014: Removed code that duplicated the MiriBadPixelModel class.
09 Jul 2014: field_def changed to dq_def.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
23 Jul 2014: DQ definitions removed from metadata.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
19 Aug 2015: Mask data are now stored in a "dq" array rather than a "data"
             array.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
06 Sep 2016: Operations work again. Test code restored.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
15 Nov 2018: New data model which uses the JWST schema, saturation.schema.
             Previous data model renamed to MiriMrsStraylightModel_CDP3.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Vincent Geers (DIAS), Steven Beard (UKATC)

"""

#import warnings
import numpy as np

# Import the MIRI reserved data quality flags and flags table class
from miri.datamodels.dqflags import master_flags, \
    FlagsTable

# Import the MIRI base data model and utilities.
# from miri.datamodels.plotting import DataModelPlotVisitor
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_badpixel_model import MiriBadPixelMaskModel
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasData

# List all classes and global functions here.
__all__ = ['MiriMrsStraylightModel', 'MiriMrsStraylightModel_CDP3']

straylight_reference_setup = \
            [(0, 'DO_NOT_USE',  'Bad pixel. Do not use.'),
             (1, 'NON_SCIENCE', 'Pixel not on science portion of detector')]
straylight_reference_flags = insert_value_column( straylight_reference_setup )


class MiriMrsStraylightModel(MiriDataModel, HasData):
    """
    
    A data model for a MIRI MRS straylight mask.
    
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
        An array containing the straylight pixel mask data. Must be 2-D.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    dq: numpy array (optional)
        FOR BACKWARDS COMPATIBILITY.
        An alias for the "data" parameter.
    detector: str (optional)
        FOR BACKWARDS COMPATIBILITY.
        The name of the detector associated with this straylight pixel data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_straylight_mrs.schema"
                                                 
    def __init__(self, init=None, data=None, dq=None, dq_def=None, detector=None,
                 **kwargs):
        """
        
        Initialises the MiriMrsStraylightModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsStraylightModel, self).__init__(init=init, **kwargs)

        # Data type is Straylight mask.
        self.meta.reftype = 'STRAYMASK'
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

        # Set the instrument detector, if provided (backwards compatibility).
        if detector is not None:
            self.meta.instrument.detector = detector
        
        # This is a reference data model.
        self._reference_model()

        # The data array is provided either in the data parameter or
        # the dq parameter (backwards compatibility).
        if data is not None:
            HasData.__init__(self, data)
        else:
            HasData.__init__(self, dq)

    def __str__(self):
        """
        
        Return the contents of the MRS straylight pixel mask object
        as a readable string.
        
        """
        strg = super(MiriMrsStraylightModel, self).__str__()
        # Add any straylight specific code here.
        return strg


class MiriMrsStraylightModel_CDP3(MiriBadPixelMaskModel):
    """
    
    This class can be used to access the old CDP-6 version of the
    MiriBadPixelMaskModel data model, which was based on the MIRI
    bad pixel model, MiriBadPixelMaskModel.
    
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
        
    dq: numpy array (optional)
        An array containing the straylight pixel mask data. Must be 2-D.
        If a dq parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the straylight pixel mask. For
        example: [(0, 'good','Good data'), (1, 'non_science', 'Non Science 
        data')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    detector: str (optional)
        The name of the detector associated with this straylight pixel data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_straylight_mrs_CDP3.schema"
    
    # Set the default dq_def table to the JWST master flags
    # TODO: Can the default declared in the schema be used?
    _default_dq_def = straylight_reference_flags
                                             
    def __init__(self, init=None, dq=None, dq_def=None, detector=None,
                 **kwargs):
        """
        
        Initialises the MiriMrsStraylightModel_CDP3 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsStraylightModel_CDP3, self).__init__(init=init, dq=dq,
                                                          dq_def=dq_def,
                                                          detector=detector,
                                                          **kwargs)

        # Data type is Straylight mask.
        self.meta.reftype = 'STRAY'
        
        # This is a reference data model.
        self._reference_model()

    def __str__(self):
        """
        
        Return the contents of the MRS straylight pixel mask object
        as a readable string.
        
        """
        strg = super(MiriMrsStraylightModel_CDP3, self).__str__()
        # Add any straylight specific code here.
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMrsStraylightModel module.")

    PLOTTING = False  
    SAVE_FILES = False
    
    dqdata1 = np.array([[0,0,1,0],
                        [0,0,1,0],
                        [0,1,0,0],
                        [0,1,0,0]])
    dqdata2 = np.array([[0,0,1,0],
                        [0,0,1,0],
                        [0,1,1,0],
                        [0,1,0,0]])

    # ------------- NEW DATA MODEL
    print("Testing the new data model: MiriMrsStraylightModel.")
    # Define a specific mask definition which is different from the default.
    print("\nFirst straylight pixel mask:")
    with MiriMrsStraylightModel( data=dqdata1 ) as testmask1:
        print(testmask1)
        if PLOTTING:
            testmask1.plot(description="testmask1")
        if SAVE_FILES:
            testmask1.save("test_mrsstraypixel_model1.fits", overwrite=True)

        print("\nSecond straylight pixel mask (use dq alias):")
        with MiriMrsStraylightModel( dq=dqdata2 ) \
                as testmask2:
            print(testmask2)
            if SAVE_FILES:
                testmask2.save("test_mrsstraypixel_model2.fits", overwrite=True)

            print("Combine the two masks:")
            result1 = testmask1.data | testmask2.data
            result2 = testmask1.data ^ testmask2.data
            result3 = testmask1.data & testmask2.data
            print(result3)
            if SAVE_FILES:
                result3.save("test_mrsstraypixel_result3.fits", overwrite=True)
            del result1, result2, result3
            del testmask2
        del testmask1

    print("\nEmpty 3x3 mask:")
    with MiriMrsStraylightModel( (3,3) ) as testmask4:
        print(testmask4)
        if PLOTTING:
            testmask4.plot(description="testmask4")
        if SAVE_FILES:
            testmask4.save("test_mrsstraypixel_model4.fits", overwrite=True)
        del testmask4


    # ------------- OLD DATA MODEL
    print("Testing the old data model: MiriMrsStraylightModel_CDP3.")
    # Define a specific mask definition which is different from the default.
    print("\nFirst straylight pixel mask - explicitly defined dq_def values:")
    with MiriMrsStraylightModel_CDP3( dq=dqdata1, dq_def=straylight_reference_flags,
                                 detector='MIRIFUSHORT' ) as testmask1:
        print(testmask1)
        print("FlagsTable:\n" + str(testmask1.flags_table))
        print("No straylight estimate values:\n" + str(testmask1.get_dq_for_field('DO_NOT_USE')))
        if PLOTTING:
            testmask1.plot(description="testmask1_cdp3")
        if SAVE_FILES:
            testmask1.save("test_mrsstraypixel_model1_cdp3.fits", overwrite=True)

        print("\nSecond straylight pixel mask:")
        with MiriMrsStraylightModel_CDP3( dq=dqdata2, dq_def=straylight_reference_flags ) \
                as testmask2:
            print(testmask2)
            print("FlagsTable:\n" + str(testmask2.flags_table))
            if SAVE_FILES:
                testmask2.save("test_mrsstraypixel_model2_cdp3.fits", overwrite=True)

            print("Combine the two masks:")
            result1 = testmask1 | testmask2
            result2 = testmask1 ^ testmask2
            result3 = testmask1 & testmask2
            print(result3)
            if SAVE_FILES:
                result3.save("test_mrsstraypixel_result3_cdp3.fits", overwrite=True)
            del result1, result2, result3
            del testmask2
        del testmask1

    print("\nMask with no dq_def values:")
    with MiriMrsStraylightModel_CDP3( dq=dqdata1, dq_def=None,
                                detector='MIRIFUSHORT' ) as testmask3:
        print(testmask3)
        print("FlagsTable:\n" + str(testmask3.flags_table))
        if PLOTTING:
            testmask3.plot(description="testmask3_cdp3")
        if SAVE_FILES:
            testmask3.save("test_mrsstraypixel_model3_cdp3.fits", overwrite=True)
        del testmask3

    print("\nEmpty 3x3 mask:")
    with MiriMrsStraylightModel_CDP3( (3,3) ) as testmask4:
        print(testmask4)
        print("FlagsTable:\n" + str(testmask4.flags_table))
        if PLOTTING:
            testmask4.plot(description="testmask4_cdp3")
        if SAVE_FILES:
            testmask4.save("test_mrsstraypixel_model4_cdp3.fits", overwrite=True)
        del testmask4

    print("Test finished.")
