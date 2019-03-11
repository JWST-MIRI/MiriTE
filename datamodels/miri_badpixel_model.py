#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines a MIRI bad 
pixel mask.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

15 Nov 2012: Created
19 Dec 2012: MIRI schema names prefixed with "miri.".
10 Jan 2013: Test file names changed.
16 Jan 2013: columnnames converted to class variable.
21 Jan 2013: Included data type in metadata.
22 Jan 2013: Bitwise operations included. Added plotting.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
13 Feb 2013: Update the .good_value attribute on the fly.
01 Mar 2013: Updated to match Jane Morrison's latest DHAS flags.
23 Apr 2013: mask element renamed data to keep up with changes
             in jwst.datamodel data model (where inclusion of a data
             element is now compulsory). mask added as an alias.
17 May 2013: Create a default table object explicitly to stop
             the underlying data model filling the table with
             junk by default.
01 Jul 2013: HDU name changed from 'MASK' to 'DQ'.
22 Aug 2013: columnnames renamed to dq_def_names.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
03 Sep 2013: Check there is a DETECTOR identification in the header.
04 Oct 2013: Data quality array changed to 32-bit. Standard bad pixel
             flags declared here. Modified to use FlagsTable object.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
11 Apr 2014: Removed redundant coordinate system setting code.
25 Jun 2014: Added dq_def alongside field_def. field_def is deprecated.
             Automatically create DQ_n header keywords from the dq_def
             binary table. 
09 Jul 2014: JWST reference flags table added.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
23 Jul 2014: DQ definitions removed from metadata.
29 Aug 2014: JWST reference flags table updated.
             Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
29 Sep 2014: Corrected a mistake in the logic of the creation of a default
             DQ_DEF table.
30 Sep 2014: Added the convert_flags method, to aid CDP-3 delivery.
             Superflous flags commented out.
27 Oct 2014: Disabled automatic creation / update of the BPCOUNT keyword.
09 Dec 2014: Obsolete field_def table removed.
18 Aug 2015: Bad pixel mask schema merged with STScI model. Obsolete
             dq_metadata schema and associated test code removed.
             Duplicated dq_def schema also removed. Mask data are now stored
             in a "dq" array rather than a "data" array.
             Removed bad pixel count metadata.
02 Sep 2015: Added subarray extraction function.
11 Sep 2015: Removed duplicated plot method.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
28 Jan 2016: Moved the .mask alias to the HasMask class.
08 Jun 2016: Removed obsolete comment.
15 Jul 2016: Added NON_SCIENCE flag.
06 Sep 2016: Operations work again. Test code restored.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
28 Jun 2018: Switch to using get_title_and_metadata() to display data model
             information.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
@author: Steven Beard (UKATC), Michael Droettboom (STScI), Vincent Geers (UKATC)

"""

import numpy as np

# Import the MIRI reserved data quality flags and flags table class
from miri.datamodels.dqflags import FlagsTable, insert_value_column, convert_dq

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasMask

# List all classes and global functions here.
__all__ = ['mask_reference_flags', 'MiriBadPixelMaskModel']

# The new JWST mask reference flags. Note that these flags are a subset
# of the JWST master flags table. The names and descriptions are the same
# but the bit values can be different.
mask_reference_setup = \
            [(0, 'DO_NOT_USE',  'Bad pixel. Do not use.'),
             (1, 'DEAD',        'Dead pixel'),
             (2, 'HOT',         'Hot pixel'),
             (3, 'UNRELIABLE_SLOPE', 'Slope variance large (=noisy pixel)'), # Was NOISY.
             (4, 'RC',          'RC pixel'),
             (9, 'NON_SCIENCE', 'Pixel not on science portion of detector')]
mask_reference_flags = insert_value_column( mask_reference_setup )


class MiriBadPixelMaskModel(MiriDataModel, HasMask):
    """
    
    A data model for a MIRI bad pixel mask, based on the STScI base model,
    DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the dq parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    dq: numpy array (optional)
        An array containing the bad pixel mask data. Must be 2-D.
        If a dq parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (bit:int, value:int, name:str, title:str),
        giving the meaning of values stored in the bad pixel mask. For
        example: [(0, 1, 'good','Good data'), (1, 2, 'dead', 'Dead Pixel'),
        (2, 4, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags plus
        standard bad pixel flags.
    detector: str (optional)
        The name of the detector associated with this bad pixel data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_bad_pixel_mask.schema.yaml"
    dq_def_names = ('BIT', 'VALUE', 'NAME', 'DESCRIPTION')
    
    # Define a default value for the dq_def table.
    # TODO: Can the default declared in the schema be used?
    _default_dq_def = mask_reference_flags
                                             
    def __init__(self, init=None, dq=None, dq_def=None, detector=None,
                 **kwargs):
        """
        
        Initialises the MiriBadPixelMaskModel class.
        
        Parameters: See class doc string.

        """
        super(MiriBadPixelMaskModel, self).__init__(init=init, **kwargs)

        # Data type is bad pixel mask.
        self.meta.reftype = 'MASK'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type:
            self.meta.model_type = model_type

        # This is a reference data model.
        self._reference_model()

        # Define the detector identifier, if specified.
        if detector is not None:
            self.meta.instrument.detector = detector
        
        # Update the dq array if it has been specifically provided.
        HasMask.__init__(self, dq)

        # Set the bit field definitions table, if provided.
        if dq_def is not None:
            try:
                self.dq_def = dq_def
            except (ValueError, TypeError) as e:
                strg = "dq_def must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        elif self.dq_def is None or len(self.dq_def) < 1:
            # No dq_def is provided.
            # Explicitly create a DQ_DEF table with default values.
            # TODO: Can the default be declared in the schema?
            self.dq_def = self._default_dq_def

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model.
        
        """
        return 'dq'

    def get_value_for_field(self, name):
        """
        
        Returns the dq value associated with a particular
        field name.

        :Parameters:

        name : str
            The name of the bit field to retrieve

        :Returns:

        field_value : int
            The dq value associated with the given field name.
            Returns None if the field name is not recognised.
            
        """
        if self.flags_table is not None and name in self.flags_table:
            return self.flags_table[name]
        else:
            return None

    def get_dq_for_field(self, name):
        """
        
        Returns an array that is `True` everywhere the given flag(s) are
        True in the dq array.

        :Parameters:

        name : str or list of str
            The name or names of the flag(s) to retrieve

        :Returns:

        array : boolean numpy array
            `True` everywhere the requested bitfield is `True`.  This
            is the same shape as the dq array.  This array is a copy
            and changes to it will not affect the underlying model.
            Returns None if the dq_def table is empty.
            
        """
        if self.flags_table is not None:
            return self.flags_table.test_flags_any( self.dq, name )
        else:
            return None
    
    def get_subarray(self, subarray):
        """
        
        Extract a subarray from the bad pixel mask.
        
        :Parameters:
        
        subarray: tuple of 4 ints
            (bottom_row, left_column, rows, columns)
            Must lie within the bad pixel mask.
            
        :Returns:
        
        subarray_dq: array-like, int
            Subarray data of size rows x columns
        
        """
        assert isinstance(subarray, (list,tuple))
        assert len(subarray) >= 4
        # Note that subarray rows and columns start at 1 but
        # numpy array indices start at 0.
        r1 = int(subarray[0]) - 1
        c1 = int(subarray[1]) - 1
        r2 = r1 + int(subarray[2])
        c2 = c1 + int(subarray[3])
        if r1 >= 0 and c1 >= 0 and r2 <= self.dq.shape[0] and \
           c2 <= self.dq.shape[1]:
           return self.dq[r1:r2, c1:c2]
        else:
            strg = "Subarray (%d, %d, %d, %d) " % subarray
            strg += "must be within bounds of bad pixel mask (%d, %d)" % \
                self.dq.shape
            raise IndexError(strg)
        
    def __str__(self):
        """
        
        Return the contents of the bad pixel mask object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()

        # Describe the mask array and dq_def table
        strg += self.get_data_str('dq', underline=True, underchar="-")
        strg += "\nColumn names: " + str(self.dq_def_names) + "\n"
        strg += self.get_data_str('dq_def', underline=True, underchar="-")
        #strg += str(self.flags_table)
        return strg

    # "flags_table" is a FlagsTable object created on the fly
    # from the contents of the dq_def table
    @property
    def flags_table(self):
        if hasattr(self, 'dq_def') and self.dq_def is not None:
            # Convert the dq_def table into a FlagsTable object
            # and return it.
            return FlagsTable( self.dq_def )
        else:
            return None

    @flags_table.setter
    def flags_table(self, data):
        raise AttributeError("The flags_table object is read-only")

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriBadPixelMaskModel module.")

    PLOTTING = False  
    SAVE_FILES = False
    
    dqdata1 = np.array([[0,8,8,0],
                        [1,3,3,1],
                        [0,5,5,0],
                        [2,0,0,2]])
    dqdata2 = np.array([[0,1,0,1],
                        [1,0,1,0],
                        [0,1,0,1],
                        [1,0,1,0]])

    # Define a specific dq definition which is different from the default.
    dqdef = [(0, 'bad',   'Bad data'),
             (1, 'dead',  'Dead Pixel'),
             (2, 'hot',   'Hot Pixel'),
             (3, 'noisy', 'Noisy pixel'),
             (4, 'flag4', 'General flag 4'),
             (5, 'flag5', 'Meaning of flag 5'),
             (6, 'flag6', 'Meaning of flag 6')
            ]
    dqdef_4 = insert_value_column( dqdef )

    print("\nFirst bad pixel mask - explicitly defined dq_def values:")
    with MiriBadPixelMaskModel( dq=dqdata1, dq_def=dqdef_4,
                                detector='MIRIMAGE' ) as testmask1:
        print(testmask1)
        print("FlagsTable:\n" + str(testmask1.flags_table))
        print("Dead values:\n" + str(testmask1.get_dq_for_field('dead')))
        print("Noisy values:\n" + str(testmask1.get_dq_for_field('noisy')))
        if PLOTTING:
            testmask1.plot(description="testmask1")
        if SAVE_FILES:
            testmask1.save("test_badpixel_model1.fits", overwrite=True)

        print("\nSecond bad pixel mask:")
        with MiriBadPixelMaskModel( dq=dqdata2, dq_def=dqdef_4 ) \
                as testmask2:
            print(testmask2)
            print("FlagsTable:\n" + str(testmask2.flags_table))
            if SAVE_FILES:
                testmask2.save("test_badpixel_model2.fits", overwrite=True)

            print("Combine the two masks:")
            result1 = testmask1 | testmask2
            result2 = testmask1 ^ testmask2
            result3 = testmask1 & testmask2
            print(result3)
            if SAVE_FILES:
                result3.save("test_badpixel_result3.fits", overwrite=True)
            del result1, result2, result3
            del testmask2
        del testmask1

    print("\nMask with no dq_def values:")
    with MiriBadPixelMaskModel( dq=dqdata1, dq_def=None,
                                detector='MIRIFUSHORT' ) as testmask3:
        print(testmask3)
        print("FlagsTable:\n" + str(testmask3.flags_table))
        if PLOTTING:
            testmask3.plot(description="testmask3")
        if SAVE_FILES:
            testmask3.save("test_badpixel_model3.fits", overwrite=True)
        del testmask3

    print("\nEmpty 3x3 mask:")
    with MiriBadPixelMaskModel( (3,3) ) as testmask4:
        print(testmask4)
        print("FlagsTable:\n" + str(testmask4.flags_table))
        if PLOTTING:
            testmask4.plot(description="testmask4")
        if SAVE_FILES:
            testmask4.save("test_badpixel_model4.fits", overwrite=True)
        field_names = testmask4.get_field_names('dq_def')
        print("The mask contains these field names:", str(field_names))
        del testmask4

    print("Test finished.")
