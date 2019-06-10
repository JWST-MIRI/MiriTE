#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

A package of compatibility functions helping to ease the transition between
old MIRI data model and the new one. The module serves two purposes:

1) To provide replacements for methods which are not implemented in the
   new data model. If any of these methods are still popular (e.g. get_title
   might be generally useful since it generates a nicely formatted string)
   they could be moved to the new data model.

2) To provide documentation on what sort of functions in the new data model
   replace the old methods and attributes. For example the .ndim attribute
   can in most cases just be replaced by .data.ndim. It isn't worth keeping
   these trivial aliases. 

:History:

31 Jan 2013: Created
06 Feb 2013: Most of the contents moved to MiriDataModel.
09 Jul 2014: field_def changed to dq_def.
09 Sep 2015: Brought up to date with changes to the dq_def table.

@author: Steven Beard (UKATC)

"""

import sys, os

import numpy as np
import numpy.ma as ma

# Import the STScI image model and utilities
import jwst.datamodels.util as jmutil
from jwst.datamodels.model_base import DataModel

from miri.datamodels.dqflags import insert_value_column

# List all classes and global functions here.
__all__ = ['WithCompatibility']

#
# Private helper functions.
#
def _truncate_string_left(strg, maxlen):
    """
        
    Helper function which truncates the left hand side of a string
    to the given length and adds a continuation characters, "...".
        
    """
    if len(strg) > maxlen:
        lhs = maxlen - 4
        strg = "... %s" % strg[-lhs:]
    else:
        return strg

def _truncate_string_right(strg, maxlen):
    """
        
    Helper function which truncates the right hand side of a string
    to the given length and adds a continuation characters, "...".
        
    """
    if len(strg) > maxlen:
        rhs = maxlen - 4
        return "%s ..." % strg[:rhs]
    else:
        return strg

def _truncate_filename(filename, maxlen):
    """
    
    Helper function which truncates a filename if it is longer than
    the given maximum length.
    
    """
    if len(filename) > maxlen:
        filebits = os.path.split(filename)
        return filebits[-1]
    else:
        return filename


class WithCompatibility(object):
    """
    
    An abstract "bolt on" class which adds various compatibility functions
    to a MIRI data model.
    
    This class has no initializer.
    
    """
    def __init__(self):
        # A dummy initializer. Having no initializer means the class
        # must not assume that any attribute exists - each one must
        # be tested.
        pass


    # The size, shape and ndim attributes map onto the size, shape and
    # ndim attributes of 'data' array, if it exists (None if it does not).
    # Modifying these attributes is not allowed.
    @property
    def size(self):
        if hasattr(self, 'data'):
            return self.data.size
        else:
            return None

    @size.setter
    def size(self, data):
        raise AttributeError("MIRI data product size attribute is read-only")

#
# The shape attribute no longer works. It must clash with something
# already implemented in the STScI model. Just use .data.shape instead.
#
#    @property
#    def shape(self):
#        print("SHAPE!")
#        if hasattr(self, 'data'):
#            return self.data.shape
#        else:
#            return None
#
#    @shape.setter
#    def shape(self, data):
#        raise AttributeError("MIRI data product shape attribute is read-only")

    @property
    def ndim(self):
        if hasattr(self, 'data'):
            return self.data.ndim
        else:
            return None

    @ndim.setter
    def ndim(self, data):
        raise AttributeError("MIRI data product ndim attribute is read-only")


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the compatibility  module.")
    
    # Get the MiriMeasuredModel class and add the compatibility bolt-on.
    from miri.datamodels.miri_measured_model import \
        MiriRampModel
    from miri.datamodels.miri_badpixel_model import \
        MiriBadPixelMaskModel
    # Add the compatibility bolt-ons
    class NewRampModel(MiriRampModel, WithCompatibility):
        # No new implementation is needed - it's all in the bolt-on.
        pass
    class NewMaskModel(MiriBadPixelMaskModel, WithCompatibility):
        # No new implementation is needed - it's all in the bolt-on.
        pass

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    data3d = [data3x3,data3x3]
    data4d = [data3d,data3d]
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
    
    print("\nMeasured ramp data with data + err + dq:")
    testdata = NewRampModel(data=data4d, err=err3x3, pixeldq=dq3x3)
#    print(testdata)
    
    print("Ramp data structure title is:\n" + testdata.get_title(underline=True))
    
    hasdata = testdata.has_dataarray('data')
    if hasdata:
        print("Ramp data has data array with title \'" + \
              testdata.get_data_title('data') + \
              "\', ndim=" + str(testdata.ndim) + \
              ", size=" + str(testdata.size) + \
              " and shape=" + str(testdata.shape))
        print("Data labels: " + testdata.get_data_labels('data'))
    else:
        print("Ramp data has no data array")

    hassci = testdata.has_dataarray('SCI')
    if hassci:
        print("Ramp data has SCI array with title \'" + \
              testdata.get_data_title('SCI') + "\'")
    else:
        print("Ramp data has no SCI array")

    bananas = testdata.has_dataarray('bananas')
    if bananas:
        print("Ramp data has bananas")
    else:
        print("Ramp data has no bananas")

    hastable = testdata.has_datatable('data')
    if hastable:
        print("Ramp data has data table")
    else:
        print("Ramp data has no data table")
    
    mask = np.array([[0,1,0],[1,0,1],[0,1,0]])
    maskdef = [(0, 'good',    'Good data'),
               (1, 'dead',    'Dead Pixel'),
               (2, 'hot',     'Hot Pixel'),
               (3, 'noisy',   'Noisy pixel'),
               (4, 'globsat', 'Above global saturation'),
               (5, 'pixsat',  'Above pixel saturation'),
               (6, 'cosmic',  'Cosmic ray detected'),
               (7, 'anomaly', 'Noise anomaly or cosmic ray'),
               (8, 'spike',   'Noise spike detected'),
               (9, 'nolin',   'No linearity correction possible'),
               (10,'outlin',  'Outside linearity correction range')
               ]
    maskdef_4 = insert_value_column( maskdef )

    print("\nMask with field values derived from list of tuples:")
    testmask = NewMaskModel( dq=mask, dq_def=maskdef_4 )
#    print(testmask1)

    print("Mask data structure title is:\n" + testmask.get_title(underline=True))

    hasdata = testmask.has_dataarray('data')
    if hasdata:
        print("Mask data has data array with title \'" + \
              testmask.get_data_title('data') + \
              "\', ndim=" + str(testmask.ndim) + \
              ", size=" + str(testmask.size) + \
              " and shape=" + str(testmask.shape))
        print("Data labels: " + testmask.get_data_labels('data'))
    else:
        print("Mask data has no data array")

    hassci = testmask.has_dataarray('MASK')
    if hassci:
        print("Mask data has MASK array with title \'" + \
              testmask.get_data_title('MASK') + "\'")
    else:
        print("Mask data has no MASK array")

    bananas = testmask.has_dataarray('bananas')
    if bananas:
        print("Mask data has bananas")
    else:
        print("Mask data has no bananas")

    hastable = testmask.has_datatable('dq_def')
    if hastable:
        print("Mask data has dq_def table with title \'" + \
              testmask.get_data_title('dq_def') + "\'")
    else:
        print("Mask data has no dq_def table")
