#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

This module contains flag definitions and binary operator functions for
the DQ data array. The module contains three different kinds of object:

1) Global constants

   * bit_to_value: A lookup table to make 2**bit calculations faster.

   * reserved_flags: A table containing the flag definitions common to
     all MIRI data.
     
2) Global functions

   The functions work with bit masks. They can be used directly or (if you
   prefer to work with flag names) using the class described in part (3).

    * make_mask: A function which makes a bit mask when given a lit of
      bits to set.
      
    * format_mask: A function which shows how to format a bit mask into
      binary (for display purposes).
      
    * raise_mask and lower_mask: Functions which raise or lower the
      specified bits in a data quality array.
      
    * test_mask_all, test_mask_any, test_mask_none: Functions which test
      whether a data quality array has elements which match all the bits
      in a bit mask, match any of the bits in a bit mask, or match none of
      the bits in a bit mask
      
    * mask_array_all, mask_array_any and mask_array_none: Functions which
      returns a numpy masked array containing science data masked in the
      places where all, any or no bits in a given bit mask are matched.
      
    * combine_quality: A function a pipeline process can use to combine
      two data quality array and make a third.
      
3) The FlagsTable class

   A class which defines an object capable of storing and using the
   flag definitions contained in a FIELD_DEF table. The class is created
   from a table
   
       flags_table = FlagsTable( dataobject.dq_def )
       
   The bit value associated with any named flag can be looked up from
   the table
   
       value = flags_table[ name ]
       
   Particular flags can be raised or lowered within a data quality array
   by calling the raiseflags or lowerflags methods. For example
   
       newdqdata = flags_table.raiseflags( olddqdata, [flag1, flag2, flag3] )
       
   A data quality array can be tested to see where particular flags are
   raised using the test_flags_all, test_flags_any or test_flags_none
   methods.
   
   See the doc strings of the individual methods and functions for more
   details.

:Reference:

http://miri.ster.kuleuven.be/bin/view/Internal/DataQualityFlags

:History:

25 Sep 2013: Created
01 Oct 2013: Converted to FlagsTable class and reserved_flags table added.
             Modified so that flags are converted to a bitmask before being
             used, which allows multiple flags to be raised, lowered or
             tested at once.
10 Apr 2014: Included some "assert" statements to check for problems.
26 Jun 2014: Added the global flags_table_to_metadata function, which
             converts the contents of a FlagsTable object into a collection
             of DQ_n metadata keywords.
09 Jul 2014: JWST master flags table added.
23 Jul 2014: Modified the flags table to match the decisions made at the
             22 July JWST pipeline working group meeting.
             master_flags table updated and group_flags table added.
29 Aug 2014: master_flags flags table updated again.
11 Sep 2014: Allow for 3 column and 4 column variations on the flags table.
             reserved_flags removed.
29 Sep 2014: Removed the assertion in insert_value_column function which
             prevented the function working with recarray and FITS_record
             objects.
30 Sep 2014: Added the convert_dq function, which converts a data quality
             array defined using one table into an array defined using a
             different table. Useful for converting CDP-2 format files to
             CDP-3. Ensure superflous blanks are stripped from keywords.
27 May 2015: Replaced pyfits with astropy.io.fits
06 Aug 2015: Updated list of PIXELDQ flags to match the JWST build 5
             pipeline specification.
17 Mar 2017: Corrected some documentation typos.
14 Nov 2017: Added some new DQ flags.
17 May 2018: Python 3: Converted dictionary keys return into a list.
12 Mar 2019: Removed use of astropy.extern.six (since Python 2 no longer used).
29 Apr 2019: Added REFERENCE_PIXEL to pixeldq_setup.
14 May 2019: Added Christophe's masking functions.
20 Jun 2019: Allow a FlagsTable to be created from a FITS_rec object.

@author: Ruyman Azzollini (DIAS), Steven Beard (UKATC), Christophe Cossou (CEA)

"""

import warnings
import copy

import numpy as np
import numpy.ma as ma

from astropy.io.fits import FITS_rec

# 1) Global constants
#
# This global list may be used a short cut to converting bit numbers
# to mask values. Looking up bit_to_value[bit] is faster than calculating
# 2**bit every time.
#
MAXBITS = 32  # The maximum number of bits in a data quality flag

bit_to_value = []
for bit in range(0, MAXBITS):
    bit_to_value.append( 2 ** bit )
#
# The above code makes creating the list more reliable against typos, but
# it hides what the list actually looks like. This is what the list would
# look like if defined directly using numbers.
#
# bit_to_value = [1, 2, 4, 8,
#                 16, 32, 64, 128,
#                 256, 512, 1024, 2048,
#                 4096, 8192, 16384, 32768,
#                 65536, 131072, 262144, 524288,
#                 1048576, 2097152, 4194304, 8388608,
#                 16777216, 33554432, 67108864, 134217728,
#                 268435456, 536870912, 1073741824]

# This global variable contains the original data quality flag definitions
# used by MIRI CDPs. NOW DEPRECATED. For details see
# http://miri.ster.kuleuven.be/bin/view/Internal/DataQualityFlags
#                  Bit  Name            Description
#                  ---  ----            -----------
# reserved_flags = [(0,  'unusable',     'Unusable data'),
#                   (1,  'non_science',  'Non science data'),
#                   (2,  'low_qual',     'Data quality below threshold'),
#                   (3,  'partial_data', 'Partial data only'),
#                   (4,  'unrel_err',    'ERR information unreliable'),
#                   (5,  'bad_sol',      'Bad Solution')
#                   ]

def insert_value_column(flagtable):
    """
    
    Convert a 3-column table of the form (BIT,NAME,DESCRIPTION) into a
    4-column table of the form ((BIT,VALUE,NAME,DESCRIPTION).
    
    :Parameters:
    
    flagtable: record array or numpy array
        A data structure containing a data quality flag table.
        A 3-column table is converted to a 4-column table.
        If the table already has 4 columns it is returned unchanged.
        
    :Returns:
    
    newflagtable: record array or numpy array
        A data structure containing a data quality flag table.
        The returned table will always have 4 columns.
            
    """
    # FIXME: The following assertions do not work, even supplemented by recarray and FITS_record.
    # Is there a generic "list-like" and "record-like" test that
    # could be used instead?
#     assert isinstance(flagtable, (tuple,list,np.recarray,pyfits.FITS_rec))
#     assert isinstance(flagtable[0], (tuple,list,np.record,pyfits.FITS_record))
    if len(flagtable[0]) == 4:
        # The table already contains a value column
        return flagtable
    elif len(flagtable[0]) == 3:
        # Insert a value column to make 4 columns.
        newflagtable = []
        for (bit, name, description) in flagtable:
            value = bit_to_value[bit]
            newflagtable.append( (bit,value,name,description) )
        return newflagtable
    else:
        strg = "Flag table should contain 3 or 4 columns"
        strg += " (it actually contains %d)" % len(flagtable[0])
        raise TypeError(strg)

def convert_to_recarray(flagtable):
    """
    
    Convert a 4-column table of the form (BIT,VALUE,NAME,DESCRIPTION)
    into a numpy record array with column types (int32,int32,str,str)
    
    :Parameters:
    
    flagtable: tuple containing record array
        A data structure containing a data quality flag table.
        
    :Returns:
    
    newflagtable: record array
        The same structure converted to a numpy record array.
            
    """
    assert isinstance(flagtable, (tuple,list))
    assert isinstance(flagtable[0], (tuple,list))
    assert (len(flagtable[0]) == 4)

    newflagtable = np.array( flagtable,
                             dtype=([('BIT', '<i4'),
                                     ('VALUE', '<u4'),
                                     ('NAME','<U40'),
                                     ('DESCRIPTION','<U80')])
                            )
    return newflagtable

# This global variable contains the master set of JWST pipeline flags
# used for the PIXELDQ and DQ arrays. The first 8 bits are also used
# for the GROUPDQ array.For details see
# http://jwst-reffiles.stsci.edu/source/data_quality.html
#
#                Bit Name                Description
#                --- ----                -----------
groupdq_setup = [(0,  'DO_NOT_USE',       'Bad pixel. Do not use.'),
                 (1,  'SATURATED',        'Pixel saturated during exposure'),
                 (2,  'JUMP_DET',         'Jump detected during exposure'),
                 (3,  'DROPOUT',          'Data lost in transmission'),
                 (4,  'RESERVED4',        'Reserved for future use'),
                 (5,  'RESERVED5',        'Reserved for future use'),
                 (6,  'RESERVED6',        'Reserved for future use'),
                 (7,  'RESERVED7',        'Reserved for future use')]
groupdq_flags = insert_value_column( groupdq_setup )
                 
pixeldq_setup = groupdq_setup + \
                [(8,  'UNRELIABLE_ERROR', 'Uncertainty exceeds quoted error'),
                 (9,  'NON_SCIENCE',
                                 'Pixel not on science portion of detector'),
                 (10, 'DEAD',             'Dead pixel'),
                 (11, 'HOT',              'Hot pixel'),
                 (12, 'WARM',             'Warm pixel'),
                 (13, 'LOW_QE',           'Low quantum efficiency'),
                 (14, 'RC',               'RC pixel'),
                 (15, 'TELEGRAPH',        'Telegraph pixel'),
                 (16, 'NONLINEAR',        'Pixel highly non-linear'),
                 (17, 'BAD_REF_PIXEL',    'Reference pixel cannot be used'),
                 (18, 'NO_FLAT_FIELD',    'Flat field cannot be measured'),
                 (19, 'NO_GAIN_VALUE',    'Gain cannot be measured'),
                 (20, 'NO_LIN_CORR',      'Linearity correction not available'),
                 (21, 'NO_SAT_CHECK',     'Saturation check not available'),
                 (22, 'UNRELIABLE_BIAS',  'Bias variance large'),
                 (23, 'UNRELIABLE_DARK',  'Dark variance large'),
                 (24, 'UNRELIABLE_SLOPE', 'Slope variance large'),
                 (25, 'UNRELIABLE_FLAT',  'Flat variance large'),
                 (26, 'OPEN',             
                                'Open pixel (counts move to adj. pixels)'),
                 (27, 'ADJ_OPEN',         'Adjacent to open pixel'),
                 (28, 'UNRELIABLE_RESET', 'Sensitive to reset anomaly'),
                 (29, 'MSA_FAILED_OPEN', 'Pixel sees light from failed open shutter'),
                 (30, 'OTHER_BAD_PIXEL', 'A catch-all flag'),
                 (31, 'REFERENCE PIXEL', 'Pixel is a reference pixel')]
pixeldq_flags = insert_value_column( pixeldq_setup )

# The master table is a combination of both of the above tables.
# The pixeldq table is now a combination of both.
# master_flags = groupdq_flags + pixeldq_flags
master_flags = pixeldq_flags


# A conversion table mapping the relationship between the old MIRI
# flag names and the new STScI flag names.
# Accept uppercase and lowercase variants.
#                        Old name --> New name
flag_conversion_table = {'unusable' : 'DO_NOT_USE',
                         'non_science' : 'NON_SCIENCE',
                         'low_qual' : 'CDP_LOW_QUAL',
                         'partial_data' : 'CDP_PARTIAL_DATA',
                         'unrel_err' : 'CDP_UNRELIABLE_ERROR',
                         'dead' : 'DEAD',
                         'hot' : 'HOT',
                         'noisy' : 'UNRELIABLE_SLOPE',
                         'UNUSABLE' : 'DO_NOT_USE',
                         'NON_SCIENCE' : 'NON_SCIENCE',
                         'LOW_QUAL' : 'CDP_LOW_QUAL',
                         'PARTIAL_DATA' : 'CDP_PARTIAL_DATA',
                         'UNREL_ERR' : 'CDP_UNRELIABLE_ERROR',
                         'DEAD' : 'DEAD',
                         'HOT' : 'HOT',
                         'NOISY' : 'UNRELIABLE_SLOPE'
                         }

#
# (2) Global functions
#
# The following series of global functions can be imported by any
# software wishing to set, test and combine binary flags. Functions
# containing the name "flag" will set and test one flag at a time.
# Functions containing the name "mask" work with several flags at once.
#
# (a) Functions provided by Ruyman

def convert_dq( dq, old_table, new_table, conversion_map=flag_conversion_table):
    """
    
    Convert a data quality array based on an old table into a new
    data quality array based on a new table, given a conversion map
    between old and new flag names.
    
    The function checks for flags defined in the old table which are
    raised in the old DQ array and raises the corresponding flags
    defined in the new table in the new DQ array.
    
    NOTE: If a flag contained in the old table is not found in the
    conversion map, or a flag defined in the conversion map is not
    defined in the new table, the data quality information contained
    in that old flag is lost.
    
    :Parameters:
    
    dq: numpy array
        An integer array containing the old data quality information,
        encoded according to the old flag table.
    old_table: FlagsTable object
        The flags table from which the old DQ array has been encoded.
    new_table: FlagsTable object
        The flags table which will be used to encode the new array.
    conversion_map: dictionary (optional)
        A dictionary giving the mapping between flag names in the old
        table to flag names in the new table.
        Defaults to the default flag conversion table for the CDP-3
        release.
        
    :Returns:
    
    new_dq: numpy array
        The new data quality array.
    
    """
    assert isinstance( old_table, FlagsTable)
    assert isinstance( new_table, FlagsTable)
    # Begin with a new data quality array full of zeros.
    new_dq = np.zeros_like( dq )
    
    print("Old data quality ranges from %d to %d" % (dq.min(), dq.max()))
    # For each flag in the old table, find where that flag is set in the
    # old DQ array, map it to a flag in the new table and set the new
    # flag in the new DQ array.
    for old_flag in old_table.flagvalues:
        old_flag = old_flag.strip() # Strip off padding at beginning and end
#         print("Checking old flag: %s" % old_flag)
        # Make sure there is a valid conversion
        if old_flag in conversion_map:
            new_flag = conversion_map[old_flag]
            new_flag = new_flag.strip() # Strip off padding at beginning and end
#             print("Translates to new flag: %s" % str(new_flag))
            if new_flag and new_flag in new_table.flagvalues:
                bitmask = bit_to_value[old_table[old_flag]]
                testmask = dq & bitmask
                where = np.where(testmask != 0)
                
                newbitmask = bit_to_value[new_table[new_flag]]
                new_dq[where] = new_dq[where] | newbitmask
#             else:
#                 print("New flag %s not found in new table." % str(new_flag))
#         else:
#             print("Old flag %s not found in conversion map" % old_flag)
        
    print("New data quality ranges from %d to %d" % (new_dq.min(), new_dq.max()))
    return new_dq
        
def flags_table_to_metadata( flags_table, meta_dq ):
    """
    
    Convert bit field definitions from table format to metadata
    keywords.
    
    :Parameters:
    
    flags_table: FlagsTable object
        The flags table to be converted
    meta_dq: schema object
        A reference to the metadata object in the schema
    
    """
    if flags_table is not None and meta_dq is not None:
        assert isinstance( flags_table, FlagsTable)
        for key in flags_table.keys:
            bitnum = flags_table[key]
            metaitem = 'def%d' % bitnum
            #print ("Setting", metaitem, "to", key)
            setattr(meta_dq, metaitem, key)

def make_mask(bitlist):
    """
    
    Function that makes a bitmask out of a list of bits

    :Parameters:
            
    bitlist: list of integer or integer
    
    :Returns:
    
    bitmask: integer
        A bitmask where all the specified bits are set.
    
    """
    # Convert a scalar value into a list
    if not isinstance(bitlist, (tuple,list)):
        bitlist = [bitlist]
    # Set each bit specified in the bitlist
    bitmask = 0
    for bit in bitlist:
        assert isinstance(bit, int) or isinstance(bit, np.integer), "Bit list contains non-integers: %s" % type(bit)
        bitmask |= bit_to_value[bit]
    return bitmask

def format_mask(bitmask):
    """
    
    Formats a numerical bitmask into a binary string.
    
    :Parameters:
    
    bitmask: integer
        A bit mask (e.g. 145)
        
    :Returns:
    
    bitstring: str
        A string displaying the bitmask in binary
        (e.g. '10010001')
    
    """
    assert isinstance(bitmask, int)
    return "{0:b}".format(bitmask)

def raise_mask(dqarr, bitmask):
    """
        
    Function that raises (sets) all the bits in 'dqarr' contained
    in the bitmask.
    
    :Parameters:
            
    dqarr: numpy array or integer
        numpy array which represents a dq plane (or part of it).
        The function also works when dqarr is a scalar integer.
    bitmask: integer
        A bit mask specifying all the bits to be logically "raised"
        in dqarr. For example,
        
        * bitmask=1 = 2**0 will raise bit 0.
        * bitmask=5 = 2**0 + 2**2 will raise bits 0 and 2.
        
    :Returns:
    
    newdqarr: numpy array or integer
        Returns array 'dqarr' with the specified bits raised in all
        elements (pixels).
    
    """
    assert isinstance(bitmask, int)
    # The bits are raised with a binary OR operation.
    return dqarr | bitmask

def lower_mask(dqarr, bitmask):
    """
        
    Function that lowers (unsets) all the bits in 'dqarr' contained
    in the bitmask.
            
    :Parameters:
            
    dqarr: numpy array or integer
        numpy array which represents a dq plane (or part of it).
        The function also works when dqarr is a scalar integer
    bitmask: integer
        A bit mask specifying all the bits to be logically "raised"
        in dqarr. For example,
        
        * bitmask=1 = 2**0 will lower bit 0.
        * bitmask=5 = 2**0 + 2**2 will lower bits 0 and 2.
        
    :Returns:
    
    newdqarr: numpy array or integer
        Returns array 'dqarr' with the specified bits lowered in all
        elements (pixels).

    """
    assert isinstance(bitmask, int)
    # The bits are lowered with a binary AND NOT operation.
    return dqarr & ~bitmask

def test_mask_all(dqarr, bitmask):
    """
        
    Probes 'dqarr' to see whether ALL the bits in the specified mask
    are raised. 
    
    :Parameters:
            
    dqarr:
        numpy array which represents a dq plane (or part of it).
    bitmask: integer
        A bit mask specifying all the bits to be tested
        in dqarr. For example,
        
        * bitmask=1 = 2**0 will test bit 0.
        * bitmask=5 = 2**0 + 2**2 test bits 0 and 2.
        
    :Returns:

    dqmask: numpy bool array
        It returns a numpy array with the same shape as 'dqarr', of bool 
        dtype, with elements set to 'True' where all the bits are set, and 
        'False' elsewhere.
            
    """
    assert isinstance(bitmask, int)
    # The bits are all raised if a binary AND operation makes no change.
    return (dqarr & bitmask) == bitmask

def test_mask_any(dqarr, bitmask):
    """
        
    Probes 'dqarr' to see whether ANY OF the bits in the specified mask
    are raised. 

    :Parameters:
            
    dqarr:
        numpy array which represents a dq plane (or part of it).
    bitmask: integer
        A bit mask specifying any of the bits to be tested
        in dqarr. For example,
        
        * bitmask=1 = 2**0 will test bit 0.
        * bitmask=5 = 2**0 + 2**2 test bits 0 or 2.
        
    :Returns:

    dqmask: numpy bool array
        It returns a numpy array with the same shape as 'dqarr', of bool 
        dtype, with elements set to 'True' where any of the bits are
        set, and 'False' elsewhere.
            
    """
    assert isinstance(bitmask, int)
    # At least one bit is all raised if a binary AND operation leaves at least
    # one bit set.
    return (dqarr & bitmask) != 0

def test_mask_none(dqarr, bitmask):
    """
        
    Probes 'dqarr' to see whether NONE OF the bits in the specified mask
    are raised. 
    
    :Parameters:
            
    dqarr:
        numpy array which represents a dq plane (or part of it).
    bitmask: integer
        A bit mask specifying any of the bits to be tested
        in dqarr. For example,
        
        * bitmask=1 = 2**0 will test bit 0.
        * bitmask=5 = 2**0 + 2**2 test bits 0 or 2.
        
    :Returns:

    dqmask: numpy bool array
        It returns a numpy array with the same shape as 'dqarr', of bool 
        dtype, with elements set to 'True' only where none of the bits
        are set 'False' ir returned elsewhere.
            
    """
    assert isinstance(bitmask, int)
    # None of the bits are set if a binary AND operation leaves nothing set.
    return (dqarr & bitmask) == 0

def mask_array_all(scidata, dqdata, bitmask=1, fill_value=None):
    """
    
    Return a new version of the given science array which is masked
    in all the places where the corresponding dq array matches ALL
    the bits in the defined bitmap.
    
    NOTE: Operations on masked numpy arrays are slower than with
    normal numpy arrays, so it is quicker to work directly on the
    unmasked arrays and propagate the data quality separately.
    However, masked arrays can be useful in situations where the
    range of the data is important, such as when scaling a plot or
    deriving statistical parameters from the data.
    
    :Parameters:
    
    scidata: numpy array
        The science data array to be masked
    dqdata: numpy array
        The data quality array to be used to generate the mask.
        This must be broadcastable onto the scidata array.
    bitmask: integer (optional)
        A bit mask with which to check the data quality array.
        If not specified, bit 0 will be checked by default.
    fill_value: number (optional)
        If specified, the value used by numpy to fill missing entries
        in the data array. If not specified, a numpy default value
        will be used.
        
    :Returns:
    
    masked_data: numpy masked array
        A masked version of the original data array.
    
    """
    maskdq = test_mask_all(dqdata, bitmask)
    return ma.array(scidata, mask=maskdq, fill_value=fill_value)

def mask_array_any(scidata, dqdata, bitmask=1, fill_value=None):
    """
    
    Return a new version of the given science array which is masked
    in all the places where the corresponding dq array matches ANY
    OF the bits in the defined bitmap.
    
    NOTE: Operations on masked numpy arrays are slower than with
    normal numpy arrays, so it is quicker to work directly on the
    unmasked arrays and propagate the data quality separately.
    However, masked arrays can be useful in situations where the
    range of the data is important, such as when scaling a plot or
    deriving statistical parameters from the data.
    
    :Parameters:
    
    scidata: numpy array
        The science data array to be masked
    dqdata: numpy array
        The data quality array to be used to generate the mask.
        This must be broadcastable onto the scidata array.
    bitmask: integer (optional)
        A bit mask with which to check the data quality array.
        If not specified, bit 0 will be checked by default.
    fill_value: number (optional)
        If specified, the value used by numpy to fill missing entries
        in the data array. If not specified, a numpy default value
        will be used.
        
    :Returns:
    
    masked_data: numpy masked array
        A masked version of the original data array.
    
    """
    maskdq = test_mask_any(dqdata, bitmask)
    return ma.array(scidata, mask=maskdq, fill_value=fill_value)

def mask_array_none(scidata, dqdata, bitmask=1, fill_value=None):
    """
    
    Return a new version of the given science array which is masked
    in all the places where the corresponding dq array matches NONE
    OF the bits in the defined bitmap.
    
    NOTE: Operations on masked numpy arrays are slower than with
    normal numpy arrays, so it is quicker to work directly on the
    unmasked arrays and propagate the data quality separately.
    However, masked arrays can be useful in situations where the
    range of the data is important, such as when scaling a plot or
    deriving statistical parameters from the data.
    
    :Parameters:
    
    scidata: numpy array
        The science data array to be masked
    dqdata: numpy array
        The data quality array to be used to generate the mask.
        This must be broadcastable onto the scidata array.
    bitmask: integer (optional)
        A bit mask with which to check the data quality array.
        If not specified, bit 0 will be checked by default.
    fill_value: number (optional)
        If specified, the value used by numpy to fill missing entries
        in the data array. If not specified, a numpy default value
        will be used.
        
    :Returns:
    
    masked_data: numpy masked array
        A masked version of the original data array.
    
    """
    maskdq = test_mask_none(dqdata, bitmask)
    return ma.array(scidata, mask=maskdq, fill_value=fill_value)

def combine_quality( dqarr1, dqarr2 ):
    """

    Combines two data quality arrays to make a third.
    
    The bitwise nature of the data quality flags means that two
    arrays can be combined without needing to know the meaning
    of the underlying flags.
    
    :Parameters:
            
    dqarr1: numpy array or None
        numpy array which represents a dq plane (or part of it).
        Can be None if there is no data quality information.
    dqarr2: numpy array or None
        numpy array which represents another dq plane (or part of it).
        Can be None if there is no data quality information.
    
    :Returns:
    
    newdq: numpy array
        numpy array containing a combination of both DQ arrays.
        Can be None if both of the input arrays are None.
    
    """
    # Check which of the arrays are defined.
    if dqarr1 is not None and dqarr2 is not None:
        # There are two quality arrays - merge them.
        # The bitwise OR operation combines the effect of the flags
        # without the need to know what they mean.
        newdq = dqarr1 | dqarr2
    elif dqarr1 is not None:
        # Only array 1 is defined - return it.
        newdq = dqarr1
    elif dqarr2 is not None:
        # Only array 2 is defined - return it.
        newdq = dqarr2
    else:
        # Neither array is defined - return None.
        newdq = None
    return newdq

#
# (b) Functions provided by Christophe
#
def change_mask( mask, include_in_mask ):
    """

    Change mask intended for a MaskedArray before actually creating the array.
    A real mask is boolean, but in MIRI data you have integers that have a
    complex structure. Each value is a combination of several status that you
    might want You can subtract a combined status, like 928 for instance
    (32, 128, 256, 512) and that will only affect pixels that have this
    specific value. Or you can give a specific status (a power of 2).
    If given that status, we will affect all concerned pixels, whether they
    possess only that status or not.
    
    See also raise_mask, lower_mask (above)
    
    :Parameters:
    
    mask: numpy array
        Array containing values. By default, 0 is correct and all others are not
        
    include_in_mask: list of int
        A list of pixel status values you want to include, for example

        - 1    : 2**0
        - 2    : 2**1
        - 4    : 2**2
        - 8    : 2**3
        - 16   : 2**4
        - 32   : 2**5
        - 64   : 2**6
        - 128  : 2**7
        - 256  : 2**8
        - 512  : 2**9
        - 1024 : 2**10
        - 2048 : 2**11
        
        See pixeldq_setup for the definition of these values.
    
    :Returns:
    
        modified mask: numpy array
    
    """
    out_mask = mask.astype(int)
    for val in include_in_mask:
        # Status can be combined and added.
        # We just want to remove one specific status for all pixels

        # If value is a power of 2, it's a unique status we want to
        # subtract to all concerned pixels (that might contain other status)
        if ((val & (val - 1)) == 0) and val > 0:
            # Identify pixels that contain this specific status
            status_pixel = np.bitwise_and(out_mask, val)
            # Subtract this status to all concerned pixels (all other have 0)
            out_mask -= status_pixel  
        else:
            out_mask[out_mask == val] = 0

    return out_mask

def combine_masks( masks ):
    """
    
    Combine a list of masks into one array. We use what np.MaskedArray uses.
    If 0, we keep the data, if not we mask it.
    
    Similar to combine_quality, except with a list of masks.

    :Parameters:
    
    masks: list of numpy arrays
        List of masks, all having the same shape. Values can be integers
    
    :Returns:
    
    combined_mask: numpy array
        One combined mask where we keep data only if visible throughout
        all the masks. Output mask will be boolean
        
    """
    combined_mask = masks[0].copy()  # By default we keep everything

    # If a pixel was previously masked, or is masked by the current image, we mask
    for mask_tmp in masks[1:]:
        combined_mask = np.logical_or(combined_mask, mask_tmp)

    return combined_mask


def decompose_mask_status( x ):
    """
    
    Given a mask value, decompose what sub-status compose it
    Example:
    One pixel mask value is 928:
    928 decompose into 32, 128, 256, 512

    :Parameters:
    
    x: int
        Mask value

    :Returns:
    
    powers: list of int
        List of powers of 2 making up the mask.

    """
    powers = []
    i = 1
    while i <= x:
        if i & x:
            powers.append(i)
        i <<= 1
    return powers

#
# (3) FlagsTable class
#
class FlagsTable(object):
    """ 
    
    Class to contain a data quality table.
    
    :Parameters:
    
    flagtable: record array or numpy array
        A data structure containing the data quality flag table.
        The object is in the same format as the DQ_DEF array
        read from a MIRI CDP.
        The table may have 3 or 4 columns.
        
        * A 3 column table is assumed to contain (BIT,NAME,DESCRIPTION).
        * A 4 column table is assumed to contain (BIT,VALUE,NAME,DESCRIPTION).
        
    """    
    def __init__(self, flagtable):    
        # The flagtable must be a list, a numpy array, a record array or a FITS record.
        if not isinstance(flagtable, (tuple,list,np.ndarray,np.recarray,FITS_rec)):
            strg = "Flag table must be a tuple, list, numpy array or FITS record"
            strg += " (%s given)" % flagtable.__class__.__name__
            raise TypeError(strg)
            
        # Convert the flag table into 2 dictionaries: one to look
        # up the flag values and a second to look up a flag description
        self.keys = []
        self.flagvalues = {}
        self.flagdescr = {}
        if flagtable is None or len(flagtable) < 1:
            # The flag table may be defined empty and populated later.
            return
        
        if len(flagtable[0]) == 3:
            # Flag table containing BIT, NAME, DESCRIPTION.
            # Calculate the VALUE field.
            for record in flagtable:
                self._check_value_range( record[0] )
                keyw = record[1].strip() # Strip off padding at beginning and end
                self.keys.append(keyw)
                self.flagvalues[keyw] = record[0]
                self.flagdescr[keyw] = record[2]
        elif len(flagtable[0]) == 4:
            # Flag table containing BIT, VALUE, NAME, DESCRIPTION
            # The VALUE field is ignored. It provides no additional information.
            for record in flagtable:
                self._check_value_range( record[0] )
                keyw = record[2].strip() # Strip off padding at beginning and end
                self.keys.append(keyw)
                self.flagvalues[keyw] = record[0]
                self.flagdescr[keyw] = record[3]
        else:
            strg = "Flag table should contain 3 or 4 columns"
            strg += " (it actually contains %d)" % len(flagtable[0])
            raise TypeError(strg)
            
        if not self._check_flags_unique():
            strg = "Flag values given to create the table are not unique"
            raise ValueError(strg)
    
    def _check_value_range(self, value):
        """
        
        Helper function which raises a ValueError if a flag value
        is out of range.
        
        """
        if value < 0 or value >= MAXBITS:
            strg = "Bit number, %d, is outside the allowed range of 0-%d" % \
                (value, MAXBITS-1)
            raise ValueError(strg)
    
    def _check_flags_unique(self):
        """
        
        Helper function which checks that all the flag values are unique.
        
        Returns True if the flag values are unique, False if otherwise.
        
        """
        # Convert the flag values dictionary into a numpy array.
        valuesarray = []
        for key in list(self.flagvalues.keys()):
            valuesarray.append(self.flagvalues[key])
        valuesarray = np.array(valuesarray)
        # Make a unique version of the array. If all the entries are
        # the same length, the two arrays must have exactly the same size.
        uniquearray = np.unique(valuesarray)
        return valuesarray.size == uniquearray.size    
            
    def __str__(self):
        """
        
        Return a string description of the flags table.
        
        """
        strg =  'NAME                 BIT     VALUE      DESCRIPTION\n'
        strg += '------------------   -----   ---------- -----------\n'
        for key in self.keys:
            bit = self.flagvalues[key]
            strg += "\'%-16s\' = 2**%2d = %10d (\'%s\')\n" % (key, bit,
                                                        bit_to_value[bit],
                                                        self.flagdescr[key])
        return strg
    
    # The following methods allow a flag table to be used like
    # a Python dictionary, to look up the value associated with
    # any named flag.
    def __setitem__(self, keyword, value):
        """
        
        Set the named flag value. Called in response to the operation::
        
            flagdata[keyword] = value
            
        Only new entries are allowed, not changes to existing entries.
        
        """
        if keyword not in self.flagvalues:
            self._check_value_range( value )
            self.keys.append(keyword)
            self.flagvalues[keyword] = value
            # Items set individually like this have an unknown description.
            self.flagdescr[keyword] = ''
        else:
            strg = "Changing an existing entry in a flags table "
            strg += "is not allowed. Only new entries may be added."
            raise KeyError(strg)

        if not self._check_flags_unique():
            strg = "Flag value given (%d) is not unique" % value
            raise ValueError(strg)
            
    def __contains__(self, keyword):
        """
        
        Returns True if the given keyword is contained in the flag table.
        Called in response to the operation::
        
            keyword in flagdata
        
        """
        return (keyword in self.flagvalues)
    
    def __delitem__(self, keyword):
        """
        
        Remove the given keyword entry from the flag table with the operation::
        
            del flagdata[keyword]

        This is not allowed.
        
        """
        strg = "Removing an entry from a flags table is not allowed."
        raise KeyError(strg)

    def __getitem__(self, keyword):
        """
        
        Return the named flag value. Called in response to the operation::
        
            value = flagdata[keyword]
            
        """
        return self.flagvalues[keyword]

    def __add__(self, other):
        """
        
        Extend the table contained in this data product. Called in response
        to the operation below. "other" must be a record array, numpy array
        or another FlagsTable object.
        
        result = this + other
        
        """
        if isinstance(other, (tuple,list,np.ndarray,np.recarray,FITS_rec)):
            # A new object is being created with an extended table.
            # Make a duplicate of the present object and extend it
            # using the new table.
            newobject = copy.deepcopy(self)
            newobject.extend( other )
            return newobject
        elif isinstance(other, FlagsTable):
            # Two flag tables are being combined together.
            newobject = copy.deepcopy(self)
            for key in other.keys:
                if key not in newobject.keys:
                    self._check_value_range( other.flagvalues[key] )
                    newobject.keys.append(key)
                    newobject.flagvalues[key] = other.flagvalues[key]
                    newobject.flagdescr[key] = other.flagdescr[key]
                else:
                    strg = "\n***Duplicate key \'%s\' ignored" % key
                    warnings.warn(strg)
            
            if not newobject._check_flags_unique():
                strg = "Flag values given to extend the table are not unique"
                raise ValueError(strg)
        
        else:
            strg = "Cannot add a FlagsTable object to a %s" % \
                other.__class__.__name__
            raise TypeError(strg)
        return newobject
    
    def __sub__(self, other):
        """
        
        Removing entries from the table by subtraction is not allowed.
        Called in response to the operation
        
        result = this - other
        
        """
        strg = "Removing entries from a flags table by subtraction "
        strg += " is not allowed."
        raise ValueError(strg)
            
    def extend(self, newtable):
        """
        
        Extend the flags table in situ by adding a new table to it.
        The new table must not duplicate existing keys.

        :Parameters:
    
        newtable: record array or numpy array
            A data structure containing the data quality flag table.
            The object is in the same format as the DQ_DEF array
            read from a MIRI CDP.
            The table may have 3 or 4 columns.
        
            * A 3 column table is assumed to contain (BIT,NAME,DESCRIPTION).
            * A 4 column table is assumed to contain (BIT,VALUE,NAME,DESCRIPTION).
        
        """
        # The new table must be a record array or a numpy array.
        if not isinstance(newtable, (tuple,list,np.array,np.recarray,FITS_rec)):
            strg = "New table must be a tuple, list or numpy array"
            raise TypeError(strg)
        
        if len(newtable[0]) == 3:
            # Flag table containing BIT, NAME, DESCRIPTION.
            # Calculate the VALUE field.
            for record in newtable:
                if record[1] not in self.keys:
                    self._check_value_range( record[0] )
                    self.keys.append(record[1])
                    self.flagvalues[record[1]] = record[0]  
                    self.flagdescr[record[1]] = record[2]
                else:
                    strg = "\n***Duplicate key \'%s\' ignored" % record[1]
                    warnings.warn(strg)  
        elif len(newtable[0]) == 4:
            # Flag table containing BIT, VALUE, NAME, DESCRIPTION
            # The VALUE field is ignored, since it provides no extra information.
            for record in newtable:
                if record[2] not in self.keys:
                    self._check_value_range( record[0] )
                    self.keys.append(record[2])
                    self.flagvalues[record[2]] = record[0]
                    self.flagdescr[record[2]] = record[3]
                else:
                    strg = "\n***Duplicate key \'%s\' ignored" % record[2]
                    warnings.warn(strg)  
        else:
            strg = "Flag table should contain 3 or 4 columns"
            strg += " (it actually contains %d)" % len(newtable[0])
            raise TypeError(strg)

        if not self._check_flags_unique():
            strg = "Flag values given to extend the table are not unique"
            raise ValueError(strg)

    # The following methods allow the flag values to be used with
    # to a data quality array.
    
    def flags_to_bitmask(self, flags):
        """
        
        Convert a given list of flags and/or flag values into
        a bitmask.
        
        :Examples:
        
        bitmask = flags_to_bitmask( 'unusable' )
        bitmask = flags_to_bitmask( ['unusable', 'partial_data'] )
        bitmask = flags_to_bitmask( 1 )
        bitmask = flags_to_bitmask( [1, 4, 8] )
        
        :Parameters:
        
        flags: string, list of strings, integer or list of integers
            The flag(s) to be converted to a bit mask.
        
            * If a string, the name of the flag to be converted to a bitmask.
            * If a list of strings, the names of all the flags to be combined
              into a bitmask.
            * If a number, the flag value (2**flag-bit) to be used. (A single
              number is assumed already to be a bitmask and is returned.)
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be combined into a bitmask.
        
        :Returns:
        
        bitmask: integer
            A bitmask with all the bits corresponding to the given
            flags set.
        
        """
        if not isinstance(flags, (tuple,list)):
            flags = [flags]
        
        bitmask = 0
        for flag in flags:
            # If the flag is a string, look up its value
            if isinstance(flag, str):
                flagbit = self[flag]
                bitmask |= bit_to_value[flagbit]
            else:
                # Otherwise assume the flag already contains a value.
                bitmask |= flag
        return bitmask
    
    def raiseflags(self, dqarr, flags):
        """
        
        Function that raises one or more flags in 'dqarr'.
            
        :Parameters:
            
        dqarr: numpy array
            numpy array which represents a dq plane (or part of it).
        flags: string, list of strings, integer or list of integers
            The flag(s) to be logically "raised" in dqarr.
            
            * If a string, the name of the flag to be raised.
            * If a list of strings, the names of all the flags to be raised.
            * If a number, the flag value (2**flag-bit) to be raised
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be raised
              
        :Returns:
        
        newdqarr: numpy array
            Returns array 'dqarr' with all elements (pixels) having the
            flag in question raised.
    
        """
        # Translate the flag name to a value
        bitmask = self.flags_to_bitmask(flags)
        return raise_mask( dqarr, bitmask)

    def lowerflags(self, dqarr, flags):
        """ 
    
        Function that lowers one or more flags in 'dqarr'.
               
        :Parameters:
            
        dqarr:
            numpy array which represents a dq plane (or part of it).
        flags: str or number
            The flag(s) to be logically "lowered" in dqarr.
            
            * If a string, the name of the flag to be lowered.
            * If a list of strings, the names of all the flags to be lowered.
            * If a number, the flag value (2**flag-bit) to be lowered
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be lowered
              
        :Returns:
        
        newdqarr: numpy array
            Returns array 'dqarr' with all elements (pixels) having the
            flag in  question lowered. If the flag wasn't raised, does
            nothing, as it should.
         
        """
        # Translate the flag name to a value
        bitmask = self.flags_to_bitmask(flags)
        return lower_mask( dqarr, bitmask)

    def test_flags_all(self, dqarr, flags):
        """
        
        Probes 'dqarr' to see whether ALL the flags in question are raised. 
        
        NOTE: If a list of flags is given this function tests if ALL
        the flags are raised. To test if any one flag is raised use
        the test_flags_any function. If only one value is given, the
        two functions give the same result.
    
        :Parameters:
            
        dqarr:
            numpy array which represents a dq plane (or part of it).
        flags: str or number
            The flag(s) to be to be "tested" against dqarr.

            * If a string, the name of the flag to be tested.
            * If a list of strings, the names of all the flags to be tested.
            * If a number, the flag value (2**flag-bit) to be tested
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be tested
              
        :Returns:
        
        dqmask: numpy array
            It returns a numpy array with the same shape as 'dqarr', of bool 
            dtype, with elements where all the flag are 'on' set as 'True',
            and 'False' elsewhere.
            
        """
        # Translate the flag name to a value
        bitmask = self.flags_to_bitmask(flags)
        return test_mask_all(dqarr, bitmask)

    def test_flags_any(self, dqarr, flags):
        """
        
        Probes 'dqarr' to see whether ANY OF the flags in question are raised. 
        
        NOTE: If a list of flags is given this function tests if ANY
        the flags are raised. To test if all of the flags are raised use
        the test_flags_all function. If only one value is given, the
        two functions give the same result.
    
        :Parameters:
            
        dqarr:
            numpy array which represents a dq plane (or part of it).
        flags: str or number
            The flag(s) to be to be "tested" against dqarr.
            
            * If a string, the name of the flag to be tested.
            * If a list of strings, the names of all the flags to be tested.
            * If a number, the flag value (2**flag-bit) to be tested
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be tested

        :Returns:
        
        dqmask: numpy array
            It returns a numpy array with the same shape as 'dqarr', of
            bool dtype, with elements where any of the flags are 'on' set
            as 'True', and  'False' elsewhere.
            
        """
        # Translate the flag name to a value
        bitmask = self.flags_to_bitmask(flags)
        return test_mask_any(dqarr, bitmask)

    def test_flags_none(self, dqarr, flags):
        """
        
        Probes 'dqarr' to see whether NONE OF the flags in question are raised. 
            
        :Parameters:
            
        dqarr:
            numpy array which represents a dq plane (or part of it).
        flags: str or number
            The flag(s) to be to be "tested" against dqarr.
            
            * If a string, the name of the flag to be tested.
            * If a list of strings, the names of all the flags to be tested.
            * If a number, the flag value (2**flag-bit) to be tested
            * If a list of numbers, the list of flag values (2**flag-bit)
              to be tested

        :Returns:
        
        dqmask: numpy array
            It returns a numpy array with the same shape as 'dqarr', of
            bool  dtype, with elements where none of the flags are 'on'
            set as 'True', and 'False' elsewhere.
            
        """
        # Translate the flag name to a value
        bitmask = self.flags_to_bitmask(flags)
        return test_mask_none(dqarr, bitmask)


#
# A minimal test and some examples of how to use the above utilities
# are run when this file is executed as a main program.
#
if __name__ == '__main__':
    print("Testing the dqflags module.")
    
    # Set up some test science data.
    scidata = np.array([[1.0,  2.0,  3.0,  4.0],
                        [8.0,  7.0,  6.0,  5.0],
                        [9.0, 10.0, 11.0, 12.0]])

    # Set up some test data quality arrays
    dqgood = np.array([[0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]])
    dqdata = np.array([[0,   0,    1,   1],
                       [1, 513, 8193, 513],
                       [1,   0, 8193,   0]])

# ----------------------------------------------------------------------------    

    # Create a collection of bit masks
    zero_one_two_mask = make_mask( [0,1,2] )
    zero_three_four_mask = make_mask( [0,3,4] )
    five_mask = make_mask( 5 )
    
    # First raise all the flags contained in the above bit masks
    print("\nBefore raising 012 mask:\n", dqgood)    
    result = raise_mask( dqgood, zero_one_two_mask)
    print("After raising 012 mask:\n", result)
    result = raise_mask( result, zero_three_four_mask)
    print("After raising 034 mask as well:\n", result)
    result = lower_mask( result, zero_one_two_mask)
    print("After lowering 012 mask:\n", result)
    
    result = raise_mask( dqgood, zero_three_four_mask)
    print("\nArray created by raising 034 mask:\n", result)
    print("Testing all bits in 034 mask:\n", test_mask_all(result, zero_three_four_mask))
    print("Testing all bits in 012 mask:\n", test_mask_all(result, zero_one_two_mask))
    print("Testing any bits in 012 mask:\n", test_mask_any(result, zero_one_two_mask))
    print("Testing any bits in 5 mask:\n", test_mask_any(result, five_mask))

# ----------------------------------------------------------------------------    

    # Create an object containing the MIRI standard flags table.
    # Note: After a MIRI CDP has been read, the same object can be created
    # by calling
    #    flagtable = FlagsTable( dataproduct.dq_def )
    reserved_table = FlagsTable(master_flags)
    print("\nThe JWST master flags table contains...")
    print(reserved_table)
    
    # Raise the 'low_qual' flag in this table
    print("\nBefore raising flag:\n", dqgood)    
    dqgood = reserved_table.raiseflags(dqgood, 'LOW_QE')
    print("After raising flag \'LOW_QE\':\n", dqgood)
    dqgood = reserved_table.raiseflags(dqgood, 'DO_NOT_USE')
    print("After raising flag \'DO_NOT_USE\':\n", dqgood)
    
    print("\nQuality array to be tested:\n", dqdata)    
    print("Where is 'DO_NOT_USE' =", bit_to_value[reserved_table['DO_NOT_USE']], \
        "raised?")
    print(reserved_table.test_flags_all(dqdata, 'DO_NOT_USE'))
    print("Where is 'NON_SCIENCE' =", bit_to_value[reserved_table['NON_SCIENCE']], \
        "raised?")
    print(reserved_table.test_flags_all(dqdata, 'NON_SCIENCE'))
    
    print("Where are 'DO_NOT_USE' and 'NON_SCIENCE' both raised?")
    print(reserved_table.test_flags_all(dqdata, ['DO_NOT_USE', 'NON_SCIENCE']))
    print("Where are either 'DO_NOT_USE' or 'partial_data' raised?")
    print(reserved_table.test_flags_any(dqdata, ['DO_NOT_USE', 'NON_SCIENCE']))
    print("Where are neither 'DO_NOT_USE' nor 'partial_data' raised?")
    print(reserved_table.test_flags_none(dqdata, ['DO_NOT_USE', 'NON_SCIENCE']))

    # Mask the science data in all the places where the 'unusable' or
    # 'partial_data' flags are raised.
    print("\nOriginal science data:\n", scidata)
    print("Data quality array:\n", dqdata)
    bitmask = reserved_table.flags_to_bitmask(['DO_NOT_USE', 'NON_SCIENCE'])
    print("Science data masked where both 'DO_NOT_USE' and 'NON_SCIENCE' are raised")
    print(mask_array_all( scidata, dqdata, bitmask ))
    print("Science data masked where either 'DO_NOT_USE' or 'NON_SCIENCE' raised")
    print(mask_array_any( scidata, dqdata, bitmask ))

# ----------------------------------------------------------------------------    

    # Try extending the flags table by adding new flags
    # NOTE: The master_flags table is now full and cannot be extended
#     bad_pixel_flags = [(30, 'choppy', 'Pixel is pretty choppy'),
#                        (31, 'noisy',  'Pixel is noisy')] #,
#     # There are three ways of making a bad pixel mask table:
#     # (1) Extend an existing table with new flags.
#     bad_pixel_table = FlagsTable( master_flags )
#     bad_pixel_table.extend( bad_pixel_flags )
#     print("\nThe first variation of the bad pixel flags table contains...")
#     print(bad_pixel_table)
# 
#     # (2) Add the new flags to an existing table.
#     combined = reserved_table + bad_pixel_flags
#     print("The second variation of the bad pixel flags table contains...")
#     print(combined)
#     del combined
#     
#     # (3) Create two separate flag table objects then add them together.
#     bad_pixel_extras = FlagsTable( bad_pixel_flags )
#     combined = reserved_table + bad_pixel_extras
#     print("The third variation of the bad pixel flags table contains...")
#     print(combined)
# 
#     # You can add a new flag as long as it is unique
#     #reserved_table['newflag'] = 2
#     reserved_table['newflag'] = 31
#     
#     # But you can't modify an existing one
#     #reserved_table['unusable'] = 23
#     
#     # And you can't create a new flag
#     # with a value that's already used    
#     # reserved_table['freakout'] = 5

# ----------------------------------------------------------------------------    

    # Insert a value field into the JWST master flags table
    table_with_values = insert_value_column(master_flags)
    print("\nInsert value column into JWST master flags table. Before...")
    print(master_flags)
    print("After...")
    print(table_with_values)

    # Create an object containing the JWST master flags table.
    master_table = FlagsTable(master_flags)
    print("\nThe JWST master flags table contains...")
    print(master_table)
    
    print("Value of DEAD is", master_table['DEAD'])
    print("Value of UNRELIABLE_ERROR is", master_table['UNRELIABLE_ERROR'])

    print("Test finished.")
