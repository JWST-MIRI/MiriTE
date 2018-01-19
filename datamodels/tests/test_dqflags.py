#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""

Simple tests of classes in "dqflags" module.

:History:

25 Sep 2013: created
03 Oct 2013: Converted into unit tests.
06 Aug 2015: Updated following changes to the JWST flag tables.

@author: Ruyman Azzollini (DIAS), Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import unittest
import numpy as np

import miri.datamodels.dqflags as dqflags


class TestDqFlagsFunctions(unittest.TestCase):
    
    # Test the MiriDroopModel class.
    
    def setUp(self):
        # Set up some test science data.
        self.scidata = np.array([[1.0,  2.0,  3.0,  4.0],
                                 [8.0,  7.0,  6.0,  5.0],
                                 [9.0, 10.0, 11.0, 12.0]])

        # Set up some test data quality arrays
        self.dqgood = np.array([[0, 0, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 0, 0]])
        self.dqdata = np.array([[0, 0, 1, 1],
                                [1, 9, 4, 2],
                                [1, 0, 8, 0]])
       
    def tearDown(self):
        # Tidy up
        del self.scidata, self.dqgood, self.dqdata

    def test_raising_lowering(self):
        # Check that the bit raising and lowering functions
        # work as expected.
        # Begin by making some bit masks.
        zero_one_two_mask = dqflags.make_mask( [0,1,2] )
        expected = 2**0 + 2**1 + 2**2
        self.assertEqual( zero_one_two_mask, expected )
        zero_three_four_mask = dqflags.make_mask( [0,3,4] )
        expected = 2**0 + 2**3 + 2**4
        self.assertEqual( zero_three_four_mask, expected )
        five_mask = dqflags.make_mask( 5 )
        expected = 2**5
        self.assertEqual( five_mask, expected )
        
        # Now mask the test data using these bit masks
        result = dqflags.raise_mask( self.dqgood, zero_one_two_mask)
        self.assertTrue( np.all( result == zero_one_two_mask) )
        result = dqflags.raise_mask( result, zero_three_four_mask)
        self.assertTrue( np.all( result == \
                                 zero_one_two_mask | zero_three_four_mask) )
        result = dqflags.lower_mask( result, zero_one_two_mask)
        self.assertTrue( np.all( result == zero_three_four_mask - 1 ) )

    def test_testing(self):
        # Check that the bit testing functions work
        zero_three_mask = dqflags.make_mask( [0,3] )
        
        data_with_both_flags = dqflags.test_mask_all( self.dqdata,
                                                      zero_three_mask )
        expected = np.array([[0, 0, 0, 0],
                             [0, 1, 0, 0],
                             [0, 0, 0, 0]], dtype=np.bool)
        self.assertTrue( np.all( data_with_both_flags == expected ) )

        data_with_either_flag = dqflags.test_mask_any( self.dqdata,
                                                      zero_three_mask )
        expected = np.array([[0, 0, 1, 1],
                             [1, 1, 0, 0],
                             [1, 0, 1, 0]], dtype=np.bool)
        self.assertTrue( np.all( data_with_either_flag == expected ) )

        data_with_neither_flag = dqflags.test_mask_none( self.dqdata,
                                                      zero_three_mask )
        expected = np.array([[1, 1, 0, 0],
                             [0, 0, 1, 1],
                             [0, 1, 0, 1]], dtype=np.bool)
        self.assertTrue( np.all( data_with_neither_flag == expected ) )
        
    def test_masking(self):
        # Check that the bit masking functions work.
        zero_three_mask = dqflags.make_mask( [0,3] )

        masked_array = dqflags.mask_array_all(self.scidata, self.dqdata,
                                             zero_three_mask)
        original_sum = np.sum( self.scidata )
        new_sum = np.sum( masked_array )
        # The new array should be the same as the old array, but with
        # the 7.0 element masked off.
        self.assertAlmostEqual( original_sum - new_sum, 7.0 )

        masked_array = dqflags.mask_array_any(self.scidata, self.dqdata,
                                             zero_three_mask)
        original_sum = np.sum( self.scidata )
        new_sum = np.sum( masked_array )
        # The new array should be the same as the old array, but with
        # the 3, 4, 8, 7, 9 and 11 elements masked off.
        self.assertAlmostEqual( original_sum - new_sum, 3+4+8+7+9+11 )

        masked_array = dqflags.mask_array_none(self.scidata, self.dqdata,
                                             zero_three_mask)
        original_sum = np.sum( self.scidata )
        new_sum = np.sum( masked_array )
        # The new array should be the same as the old array, but with
        # the 1, 2, 6, 5, 10 and 12 elements masked off.
        self.assertAlmostEqual( original_sum - new_sum, 1+2+6+5+10+12 )

class TestMiriFlagsTable(unittest.TestCase):
    
    # Test the MiriFlagsTable class.
    
    def setUp(self):
        # Create an example flag table using the MIRI reserved flags
        self.flagtable = dqflags.FlagsTable( dqflags.master_flags )
       
    def tearDown(self):
        # Tidy up
        del self.flagtable

    def test_creation(self):
        # A FlagsTable object must be created from a list or numpy array
        # containing a table.
        self.assertRaises(TypeError, dqflags.FlagsTable, None )
        self.assertRaises(TypeError, dqflags.FlagsTable, 42 )
        self.assertRaises(TypeError, dqflags.FlagsTable, 'Silly' )

        # The table must contain records, not scalars
        self.assertRaises(TypeError, dqflags.FlagsTable, [1,2,3] )
        
        # The table must contain unique values
        test_flags = [(0,    'unusable',     'Unusable data'),
                      (1,    'non_science',  'Non science data'),
                      (0,    'duplicate',    'Duplicated value')]
        self.assertRaises(ValueError, dqflags.FlagsTable, test_flags )
        
        # The table can be empty and appended to later
        emptytable = dqflags.FlagsTable( [] )
        descr = str(emptytable)
        self.assertIsNotNone(descr)
        del descr
        emptytable.extend( dqflags.master_flags )
        descr = str(emptytable)
        self.assertIsNotNone(descr)
        del descr
        del emptytable, test_flags

    def test_description(self):
        # Test that the querying and description functions work.
        # For this test to pass these only need to run without error.
        descr = str(self.flagtable)
        self.assertIsNotNone(descr)
        del descr
        descr = repr(self.flagtable)
        self.assertIsNotNone(descr)
        del descr
        
    def test_getting_setting(self):
        # The table can be added to, as long as new flags are added
        self.flagtable['newflag'] = 31
        
        # Attempt to define an out of range value
        #self.flagtable['bigflag'] = 32
        #self.flagtable['negflag'] = -1
        self.assertRaises(ValueError, self.flagtable.__setitem__, 'bigflag', 32)
        self.assertRaises(ValueError, self.flagtable.__setitem__, 'negflag', -1)   
        
        # Attempt to overwrite an existing flag
        #self.flagtable['unusable'] = 29
        self.assertRaises(KeyError, self.flagtable.__setitem__, 'DO_NOT_USE', 29)
        # Attempt to duplicate an existing value
        #self.flagtable['yetanotherflag'] = 11
        self.assertRaises(ValueError, self.flagtable.__setitem__,
                          'yetanotherflag', 11)
        # Deletion is not allowed
        #del self.flagtable['unusable']
        self.assertRaises(KeyError, self.flagtable.__delitem__, 'DO_NOT_USE')
        
        # Get the value of some existing items
        value = self.flagtable['DO_NOT_USE']
        self.assertEqual( value, 0 )
        value = self.flagtable['newflag']
        self.assertEqual( value, 31 )
        
        # Get the value of a non-existent flag
        #value = self.flagtable['nosuchflag']
        self.assertRaises(KeyError, self.flagtable.__getitem__, 'nosuchflag')
        
# Test no longer possible now self.flagstable contains the maximum
# number of flags.
#     def test_operations(self):
#         # Test the arithmetic operations
#         # It is possible to extend an existing table by adding a new one
#         bad_pixel_flags = [(29, 'choppy', 'Pixel is pretty choppy'),
#                            (30, 'noisy',  'Pixel is noisy')] #,
#         newtable = self.flagtable + bad_pixel_flags
#         descr = str(newtable)
#         self.assertIsNotNone(descr)
#         
#         value = newtable['choppy']
#         self.assertEqual(value, 29)
#         value = newtable['noisy']
#         self.assertEqual(value, 30)
#         del newtable
#         
#         # The same operation is possible by adding two FlagsTable objects
#         bad_pixel_table = dqflags.FlagsTable( bad_pixel_flags )
#         newtable = self.flagtable + bad_pixel_table
#         descr = str(newtable)
#         self.assertIsNotNone(descr)
#         
#         value = newtable['choppy']
#         self.assertEqual(value, 29)
#         value = newtable['noisy']
#         self.assertEqual(value, 30)
#         del newtable
        
    def test_raising_lowering(self):
        # Test the use of the flags table to raise and lower flags
        # in a data array.
        datashape = (4,6)  # shape of SCI and DQ arrays
        ndreg = [2,3,0,5]   # [x0,x1,y0,y1] region where the CDP is not defined. 
        badsolpix = (1,4)    # place where there is a bad solution
    
        # Start with a SCI array full of ones and a DQ array full of zeros
        sci = np.ones(datashape,dtype='float32') # SCI array 
        dq = np.zeros(datashape,dtype='int32') # DQ array

        # Create a zero-filled copy of the SCI array and fill it with
        # bad values over the undefined region. Then get the index
        # numbers for that region
        fksci = sci * 0.
        fksci[ndreg[0]:ndreg[1]+1,ndreg[2]:ndreg[3]+1] = np.nan
        ixndreg = np.where(np.isnan(fksci))

        # Zero the SCI array where the CDP is not defined.
        sci[ixndreg] = 0. 

        # A bad solution for a single pixel
        sci[badsolpix] = np.nan 
        
        # The dq array starts off full of zeros.
        self.assertEqual( dq.sum(), 0 )
        # Raise the 'non_science' and 'unusable' flags in all the
        # places where the CDP is not allowed.
        dq[ixndreg] = self.flagtable.raiseflags(dq[ixndreg],
                                                ['SATURATED', 'DO_NOT_USE'])
        # The operation should have set the bottom two rows
        # to 2**0 + 2**1 = 3
        expected = [[0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [3, 3, 3, 3, 3, 3],
                    [3, 3, 3, 3, 3, 3]]
        self.assertTrue( np.all( dq == expected ))

        # Raise the "bad solution" and "unusable" flags where the
        # SCI array is defined as NaN        
        dq[np.where(np.isnan(sci))] = \
            self.flagtable.raiseflags(dq[np.where(np.isnan(sci))],
                                      ['NON_SCIENCE', 'DO_NOT_USE']) 
        # The operation should have added a 513 at element [1,4]
        expected = [[0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0,513, 0],
                    [3, 3, 3, 3, 3, 3],
                    [3, 3, 3, 3, 3, 3]]
        self.assertTrue( np.all( dq == expected ))

# If being run as a main program, run the tests.
if __name__ == '__main__':
    unittest.main()
