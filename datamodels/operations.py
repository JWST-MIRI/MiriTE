#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Arithmetic and binary operator functions for the MIRI data model.

:Reference:

The STScI jwst.datamodels documentation.

:History:

24 Jan 2011: Created
22 Feb 2013: Zero dimensional arrays cannot be masked.
08 Oct 2013: _shrink_dq function added, to allow masking when a DQ array
             is larger than a data array.
31 Oct 2013: Improved memory management by starting a mathematical
             operation with an empty object rather than a copy.
             Corrected the formula used by _combine_errors_divisive.
11 Dec 2013: Mask an array only when its data quality value contains
             an odd number. Corrected a typo in _generate_mask().
21 May 2014: Make sure data quality arrays are integer before using
             bitwise operations.
25 Sep 2014: reserved_flags replaced by master_flags.
30 Nov 2015: Tightened up a few data type conversions, to ensure that
             bit masks have the same data type before being combined.
28 Jan 2016: Changed HasMask to use the .dq attribute instead of .mask
             (which is now defined as an alias).
23 Mar 2016: Documentation correction.
06 Apr 2016: Replaced throughout the use of _real_cls() by __class__(),
             following changes to jwst.datamodels.model_base.DataModel.
04 May 2016: noerr option added to HasDataErrAndDq.
12 Jul 2017: set_data_fill and set_err_fill options added to HasDataErrAndDq.
27 Jun 2018: Added HasDataErrAndGroups class to be used with ramp data.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""

import sys
from astropy.extern import six

import numpy as np
import numpy.ma as ma

#import miri.datamodels.dqflags as dqflags
from miri.datamodels.dqflags import master_flags, combine_quality

# Import the STScI image model and utilities
import jwst.datamodels.util as jmutil
from jwst.datamodels.model_base import DataModel

# List all classes and global functions here.
__all__ = ['HasMask', 'HasData', 'HasDataErrAndDq']


class HasMask(object):
    """
    
    An abstract class which provides the binary operations relevant for 
    data models containing a primary mask array.
    
    The primary mask array is assumed to be stored in an attribute 
    called dq.
    
    """
    def __init__(self, dq):
        if dq is not None:
            self.dq = dq
 
    # "mask" is an alias for the "dq" attribute.
    @property
    def mask(self):
        if hasattr(self, 'dq'):
            return self.dq
        else:
            return None

    @mask.setter
    def mask(self, dq):
        self.dq = dq

    def _check_for_mask(self):
        """
        
        Helper function which raises an exception if the object
        does not contain a valid data array.
        
        """
        if not self._isvalid(self.dq):
            strg = "%s object does not contain a valid mask array" % \
                self.__class__.__name__
            raise AttributeError(strg)

    def _isvalid(self, data):
        """
        
        Helper function to verify that a given array, tuple or list is
        not empty and has valid content.
        
        """
        if data is None:
            return False
        elif isinstance(data, (list,tuple)):
            if len(data) <= 0:
                return False
            else:
                return True
        elif isinstance(data, (np.ndarray)):
            if data.size <= 0:
                return False
            else:
                return True
        elif not data:
            return False
        else:
            return True

    def __or__(self, other):
        """
        
        Bitwise OR operation between this mask and another
        data product or scalar.
        
        """
        # Check this object is capable of binary operation.
        self._check_for_mask()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being operated with.
            # This is only a sensible operation when the scalar
            # quantity is converted to an integer.
            newobject.dq = self.dq | int(other)
            
        elif isinstance(other, (np.ndarray,list,tuple)):
            # A data array is being combined with this product. This should
            # work provided the two arrays are broadcastable.
            newobject.dq = self.dq | np.asarray(other, dtype=self.dq.dtype)

        elif isinstance(other, DataModel) and hasattr(other, 'dq'):
            # Two mask data products are being combined together.
            newobject.dq = self.dq | other.dq

        else:
            strg = "Cannot bitwise combine " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject

    def __xor__(self, other):
        """
        
        Bitwise EXCLUSIVE OR operation between this mask and another
        data product or scalar.
        
        """
        # Check this object is capable of binary operation.
        self._check_for_mask()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being operated with.
            # This is only a sensible operation when the scalar
            # quantity is converted to an integer.
            newobject.dq = self.dq ^ int(other)
            
        elif isinstance(other, (np.ndarray,list,tuple)):
            # A data array is being combined with this product. This should
            # work provided the two arrays are broadcastable.
            newobject.dq = self.dq ^ np.asarray(other, dtype=self.dq.dtype)

        elif isinstance(other, DataModel) and \
             hasattr(other, 'dq') and self._isvalid(other.dq):
            # Two mask data products are being combined together.
            newobject.dq = self.dq ^ other.dq

        else:
            strg = "Cannot bitwise combine " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject

    def __and__(self, other):
        """
        
        Bitwise AND operation between this mask and another
        data product or scalar.
        
        """
        # Check this object is capable of binary operation.
        self._check_for_mask()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being operated with.
            # This is only a sensible operation when the scalar
            # quantity is converted to an integer.
            newobject.dq = self.dq & int(other)
            
        elif isinstance(other, (np.ndarray,list,tuple)):
            # A data array is being combined with this product. This should
            # work provided the two arrays are broadcastable.
            newobject.dq = self.dq & np.asarray(other, dtype=self.dq.dtype)

        elif isinstance(other, DataModel) and \
             hasattr(other, 'dq') and self._isvalid(other.dq):
            # Two mask data products are being combined together.
            newobject.dq = self.dq & other.dq

        else:
            strg = "Cannot bitwise combine " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject


class HasData(object):
    """
    
    An abstract class which provides the arithmetic operations
    relevant for data models containing a primary data array.
    
    The primary data array is assumed to be stored in an attribute
    called data.
    
    """
    def __init__(self, data):
        if data is not None:
            self.data = data

    def _check_for_data(self):
        """
        
        Helper function which raises an exception if the object
        does not contain a valid data array.
        
        """
        if not self._isvalid(self.data):
            strg = "%s object does not contain a valid data array" % \
                self.__class__.__name__
            raise AttributeError(strg)

    def _isvalid(self, data):
        """
        
        Helper function to verify that a given array, tuple or list is
        not empty and has valid content.
        
        """
        if data is None:
            return False
        elif isinstance(data, (list,tuple)):
            if len(data) <= 0:
                return False
            else:
                return True
        elif isinstance(data, (ma.masked_array,np.ndarray)):
            if data.size <= 0:
                return False
            else:
                return True
        elif not data:
            return False
        else:
            return True

    def __add__(self, other):
        """
        
        Add a scalar, an array or another MiriMeasuredModel object to
        this MiriMeasuredModel object.
        
        """
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
        
        if isinstance(other,(float,int)):
            # A scalar quantity is being added.
            newobject.data = self.data + other
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being added to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data + np.asarray(other)
            
        elif isinstance(other, DataModel):
            # Two data products are being added together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data + other.data
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot add " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject

    def __sub__(self, other):
        """
        
        Subtract a scalar, an array or another MiriMeasuredModel object 
        from this MiriMeasuredModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being subtracted.
            newobject.data = self.data - other
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being subtracted to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data - np.asarray(other)
            
        elif isinstance(other, DataModel):
            # Two data products are being subtracted together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data - other.data
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot subtract " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    def __mul__(self, other):
        """
        
        Multiply this MiriMeasuredModel object by a scalar, an array or
        another MiriMeasuredModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being multiplied.
            newobject.data = self.data * other
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data * np.asarray(other)
            
        elif isinstance(other, DataModel):
            # Two data products are being multiplied together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data * other.data
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot multiply " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject
       
    def __truediv__(self, other):
        """
        
        Divide this MiriMeasuredModel object by a scalar, an array or
        another MiriMeasuredModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
 
        if isinstance(other,(float,int)):
            # A scalar quantity is being divided.
            if np.abs(other) <= sys.float_info.epsilon:
                strg = "%s: Divide by scalar zero!" % self.__class__.__name__
                del newobject
                raise ValueError(strg)
            newobject.data = self.data / other
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            newobject.data = self.data / np.asarray(other)
            
        elif isinstance(other, DataModel):
            # The data product is being divided by another. Ensure they
            # both have a valid primary data array.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data / other.data
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot divide " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    # In Python 3, division is the same as true division.
    def __div__(self, other):
        return self.__truediv__(other)


class HasDataErrAndDq(HasData):
    """
    
    An abstract class which provides the arithmetic operations and 
    masking functions relevant for data models containing a data array, 
    error array and data quality array.

    The primary, error and quality arrays are assumed to be stored in
    attributes called data, err and dq.

    """
    def __init__(self, data, err, dq, noerr=False):
        super(HasDataErrAndDq, self).__init__(data=data)
        self._data_mask = None
        self._data_fill = 0.0
        self._data_fill_value = None

        self.noerr = noerr
        if not self.noerr:
            if err is not None:
                self.err = err
            self._err_mask = None
            self._err_fill = 'max'
            self._err_fill_value = None

        if dq is not None:
            self.dq = dq

    def set_data_fill(self, data_fill):
        """
        
        Set the data fill instruction to something other than the default
        of 0.0.
    
        :Parameters:
        
        data_fill: str or number
            An instruction for how to fill the missing values within
            a masked array:
        
            * 'min': Fill with the minimum value.
            * 'max': Fill with the maximum value.
            * 'mean': Fill with the mean value
            * 'median': Fill with the median value
            * '': Fill with the default numpy value.
            * Any other value is assumed to be the fill value.
        
        """
        self._data_fill = data_fill

    def set_err_fill(self, err_fill):
        """
        
        Set the error fill instruction to something other than the default
        of 'max'.
    
        :Parameters:
        
        err_fill: str or number
            An instruction for how to fill the missing values within
            a masked array:
        
            * 'min': Fill with the minimum value.
            * 'max': Fill with the maximum value.
            * 'mean': Fill with the mean value
            * 'median': Fill with the median value
            * '': Fill with the default numpy value.
            * Any other value is assumed to be the fill value.
        
        """
        self._err_fill = err_fill

    def _shrink_dq(self, dqarray):
        """
        
        Helper function which shrinks a data quality array along
        its highest axis to generate a new array of smaller size.
        
        For example, a 3-D array of shape (3 x 3 x 2) is shrunk to
        a 2-D array of shape (3 x 3). Quality flags are combined
        in a bitwise manner.
        
        """
        # Ensure the input array is of unsigned integer type
        dqarray = np.asarray(dqarray, dtype=np.uint)
        # The new shape has the highest dimension removed
        newshape = dqarray.shape[1:]
        # Start with a DQ array full of zeros
        newdq = np.zeros( newshape, dtype=np.uint)
        # Split the data quality array along the highest
        # axis into a list of pieces.
        npieces = dqarray.shape[0]
        for piece in np.split(dqarray, npieces, 0):
            # Convert each piece into an N-1 dimensional array of integers.
            # Each should be the same size and shape as the new DQ array.
            npiece = np.asarray( np.squeeze( piece ), dtype=np.uint)
            # Merge each new piece into the new DQ array with a bitwise OR
            newdq |= npiece
        # The result should be a new mask with reduced dimensionality
        return newdq        

    def _generate_mask(self, data, dq, bitmask=1):
        """
    
        Use the contents of the dq array to generate a numpy mask of the
        same shape as the data array.
    
        :Parameters:
    
        data: numpy array
            The data array to be masked
        dq: numpy array
            The data quality array to be used to generate the mask
        bitmask: unsigned int
            If specified, a mask for selecting particular bits
            from the data quality values.
            The default of 1 will match only bit zero.
            None will match any non-zero data quality value.
        
        :Returns:
    
        mask: numpy mask
            A mask which can be used with the data array.
    
        """
#         print("+++ Generating mask from", data, "\nand", dq,
#               "\nwith bitmask", bitmask)
        # A mask can only be generated when both arrays exist and
        # are not empty. The DATA array and DQ array must also be
        # broadcastable.
        if self._isvalid(data) and dq is not None:
            # Ensure the data quality array is of unsigned integer type
            # so bitwise operations are possible.
            dq = np.asarray(dq, dtype=np.uint)
            if data.ndim < dq.ndim and jmutil.can_broadcast(dq.shape, data.shape):
                # The DQ array is larger than the array being masked.
                # This is a special case.
                # Shrink down the DQ array until the dimensions match.
                shrunk_dq = self._shrink_dq(dq)
                while (shrunk_dq.ndim > data.ndim):
                    shrunk_dq = self._shrink_dq(shrunk_dq)

                # Start with a zero (False) mask and mask off (set to True)
                # all the pixels indicated by the DQ array.
                maskdq = np.zeros(data.shape, dtype=np.bool)
                if bitmask is None:
                    # None means all bits set.
                    bad = np.where(shrunk_dq != 0)
                else:
                    bad = np.where((shrunk_dq & bitmask) != 0)
                maskdq[bad] = True
                return maskdq
      
            elif data.size >= dq.size and jmutil.can_broadcast(data.shape, dq.shape):
                # Broadcast the DQ array onto something the same shape
                # as the data array.
                datadq = np.zeros(data.shape, dtype=np.uint) + dq
                # Start with a zero (False) mask and mask off (set to True)
                # all the pixels indicated by the DQ array.
                maskdq = np.zeros(data.shape, dtype=np.bool)
                if bitmask is None:
                    # None means all bits set.
                    bad = np.where(datadq != 0)
                else:
                    bad = np.where((datadq & bitmask) != 0)
                maskdq[bad] = True
                return maskdq
            else:
                return ma.nomask # or None
        else:
            return ma.nomask # or None

    def _generate_fill(self, data, fill_descr):
        """
    
        Generate a fill value for a data array based on the masked array
        plus a fill description.
    
        :Parameters:
    
        data: numpy array
            The data array to be examined.
        fill_descr: str or number
            An instruction for how to fill the missing values within
            a masked array:
        
            * 'min': Fill with the minimum value.
            * 'max': Fill with the maximum value.
            * 'mean': Fill with the mean value
            * 'median': Fill with the median value
            * '': Fill with the default numpy value.
            * Any other value is assumed to be the fill value.
        
        :Returns:
    
        fill_value: number
            The fill value
    
        """
        # The data array must exist and must not be empty.
        if self._isvalid(data):
            if isinstance(fill_descr, six.string_types):
                if fill_descr == 'min':
                    # Use the minimum unmasked value as the fill value
                    fill_value = data.min()
                elif fill_descr == 'max':
                    # Use the maximum unmasked value as the fill value
                    fill_value = data.max()
                elif fill_descr == 'mean':
                    # Use the mean unmasked value as the fill value
                    fill_value = data.mean()
                elif fill_descr == 'median':
                    # Use the median unmasked value as the fill value
                    fill_value = data.median()
                else:
                    # Use the default numpy fill value
                    fill_value = None
            else:
                # Assume the fill description is a number or None
                fill_value = fill_descr
        else:
            fill_value = None
        return fill_value

    def _mask_array(self, data, dq, fill_value=None):
        """
    
        Return a masked version of the given array.
        
        NOTE: This function might introduce small rounding errors into
        floating point data, so a value displayed as 3.00000005 before
        masking might display as 3.000000048 afterwards. The difference
        is insignificant, but it looks worse when displayed.
    
        :Parameters:
    
        data: numpy array
            The data array to be masked
        dq: numpy array
            The data quality array to be used to generate the mask
        fill_value: number
            If specified, the value used to fill missing entries in the
            data array. If not specified, a numpy default value will be
            used.
        
        :Returns:
    
        masked_data: numpy masked array
            A masked version of the original data array.
    
        """
        maskdq = self._generate_mask(data, dq)
        return ma.array(data, mask=maskdq, fill_value=fill_value)

    def _combine_errors_maximum(self, error1, error2):
        """
        
        Helper function to combine two error arrays and return the maximum.
        Can be used when two data arrays are combined with a min or max
        function, or are combined by resampling.
        
        NOTE: This function is valid only when both error arrays are sampling
        the same error source and you prefer to believe the most pessimistic
        estimate. Use with care.
        
        """
        # The end product will have an ERR unit only if both products
        # started with an ERR unit.
        if error1 is not None and error2 is not None:
            newerr = np.maximum(error1, error2)
        else:
            newerr = None
        return newerr

    def _combine_errors_quadrature(self, error1, error2):
        """
        
        Helper function to combine two error arrays in quadrature.
        Can be used when two data arrays are added or subtracted.
        
        NOTE: This function is valid only when combining two sets
        of data with independent errors. This assumption might not
        be valid in all circumstances, so use with care.
        
        """
        # The end product will have an ERR unit only if both products
        # started with an ERR unit.
        if error1 is not None and error2 is not None:
            # NOTE: These operations might cause an overflow
            # for some data types.
            err1sq = np.square(error1)
            err2sq = np.square(error2)
            sumsq = err1sq + err2sq
            newerr = np.sqrt(sumsq)
        else:
            newerr = None
        return newerr

    def _combine_errors_multiplicative(self, error1, error2, data1, data2):
        """
        
        Helper function to combine two error arrays in quadrature,
        where each error array is weighted by a sensitivity
        coefficient.
        This functions can be used when two data arrays are multiplied, 
        so the sensitivity coefficient is proportional to the other 
        array's measurement data.

        NOTE: This function is valid only when combining two sets
        of data with independent errors. This assumption might not
        be valid in all circumstances, so use with care.
        
        """
        # The end product will have an ERR unit only if both products
        # started with an ERR unit.
        if error1 is not None and error2 is not None:
            if data1 is not None and data2 is not None:
                # NOTE: These operations might cause an overflow
                # for some data types.
                data1sq = np.square(data1)
                data2sq = np.square(data2)
                err1sq = np.square(error1)
                err2sq = np.square(error2)
                sumsq = (data2sq * err1sq) + (data1sq * err2sq)
                #newerr = np.sqrt(sumsq) / (data1sq+data2sq) ???
                newerr = np.sqrt(sumsq)
            else:
                # Without the data arrays the weighting is unknown.
                return self._combine_errors_quadrature(error1, error2)
        else:
            newerr = None
        return newerr

    def _combine_errors_divisive(self, error1, error2, data1, data2):
        """
        
        Helper function to combine two error arrays in quadrature,
        where each error array is weighted by a sensitivity
        coefficient.
        This functions is used when one data array is divided by
        another, so the sensitivity coefficient for the first array
        is proportional to the inverse of the second but the
        sensitivity coefficient for the second array is proportional
        to the first.
        
        CHECK THE MATHS

        NOTE: This function is valid only when combining two sets
        of data with independent errors. This assumption might not
        be valid in all circumstances, so use with care.
        
        """
        # The end product will have an ERR unit only if both products
        # started with an ERR unit.
        if error1 is not None and error2 is not None:
            if data1 is not None and data2 is not None:
                # NOTE: These operations might cause an overflow
                # for some data types.
                data1sq = np.square(data1)
                data2sq = np.square(data2)
                # NOTE: The errors will blow up if any of the data2sq values
                # are close to zero. There might be a divide by zero.
                err1sq = np.square(error1)
                err2sq = np.square(error2)
                sumsq = (err1sq / data2sq) + \
                        ((err2sq * data1sq) / (data2sq * data2sq))
#                 sumsq = (data2weight * err1sq) + (data1sq * err2sq)
                
                # Comment by Juergen Schreiber:
                # Shouldn't the error propagation according to Gauss be
                # sqrt(err1sq*sci2weight + err2sq*sci1sq/(sci2sq*sci2sq))
                # since the partial derivation of a/b on b is -a/(b*b)
                newerr = np.sqrt(sumsq)
            else:
                # Without the data arrays the weighting is unknown.
                return self._combine_errors_quadrature(error1, error2)
        else:
            newerr = None
        return newerr

    def _combine_quality(self, dq1, dq2):
        """
        
        Helper function to combine the quality arrays of two
        MiriMeasuredModel objects. Any point flagged as bad in
        either of the two products is flagged as bad in the
        result.
        
        """
        return combine_quality(dq1, dq2)

    def __add__(self, other):
        """
        
        Add a scalar, an array or another DataModel object to
        this MiriMeasuredModel object.
        
        """
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
        
        if isinstance(other,(float,int)):
            # A scalar quantity is being added. Add to the SCI array but
            # leave the ERR and DQ arrays as they are.
            newobject.data = self.data + other
            if not self.noerr:
                newobject.err = self.err
            newobject.dq = self.dq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being added to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data + np.asarray(other)
            # Adding a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.dq = self.dq
             
        elif isinstance(other, DataModel):
            # Two data products are being added together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data + other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = \
                            self._combine_errors_quadrature(self.err,
                                                            other.err)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)
                    
                if hasattr(other, 'dq') and self._isvalid(other.dq):
                    newobject.dq = self._combine_quality(self.dq, other.dq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot add " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject

    def __sub__(self, other):
        """
        
        Subtract a scalar, an array or another DataModel object from
        this MiriMeasuredModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being subtracted. Subtract from the SCI
            # array but leave the ERR and DQ arrays as they are.
            newobject.data = self.data - other
            if not self.noerr:
                newobject.err = self.err
            newobject.dq = self.dq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being subtracted to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data - np.asarray(other)
            # Adding a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.dq = self.dq
            
        elif isinstance(other, DataModel):
            # Two data products are being subtracted. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data - other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = \
                            self._combine_errors_quadrature(self.err, other.err)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)
                    
                if hasattr(other, 'dq') and self._isvalid(other.dq):
                    newobject.dq = self._combine_quality(self.dq, other.dq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot subtract " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    def __mul__(self, other):
        """
        
        Multiply this MiriMeasuredModel object by a scalar, an array or
        another DataModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being multiplied. Multiply the SCI and ERR
            # arrays but leave the DQ array as it is.
            newobject.data = self.data * other
            if not self.noerr:
                newobject.err = self.err * other
            newobject.dq = self.dq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data * np.asarray(other)
            # Multiplying a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.dq = self.dq
            
        elif isinstance(other, DataModel):
            # Two data products are being multiplied together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data * other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = self._combine_errors_multiplicative( \
                                            self.err, other.err, self.data,
                                            other.data)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)

                if hasattr(other, 'dq') and self._isvalid(other.dq):
                    newobject.dq = self._combine_quality(self.dq, other.dq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot multiply " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject
       
    def __truediv__(self, other):
        """
        
        Divide this MiriMeasuredModel object by a scalar, an array or
        another DataModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being divided. Divide the SCI and ERR
            # arrays but leave the DQ array as it is.
            # Trap a divide by zero..
            if np.abs(other) <= sys.float_info.epsilon:
                strg = "%s: Divide by scalar zero!" % self.__class__.__name__
                del newobject
                raise ValueError(strg)
            newobject.data = self.data / other
            if not self.noerr:
                newobject.err = self.err / other
            newobject.dq = self.dq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            newobject.data = self.data / np.asarray(other)
            # Dividing by a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.dq = self.dq
            
        elif isinstance(other, DataModel):
            # The data product is being divided by another. Ensure they
            # both have a valid primary data array.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data / other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = self._combine_errors_divisive( \
                                            self.err, other.err, self.data,
                                            other.data)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)

                if hasattr(other, 'dq') and self._isvalid(other.dq):
                    newobject.dq = self._combine_quality(self.dq, other.dq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot divide " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    # From Python 3, division is the same as true division.
    def __div__(self, other):
        return self.__truediv__(other)

    @property
    def data_masked(self):
        # Generate the masked data on the fly. This ensures the
        # masking is always up to date with the latest dq array.
        # TODO: Can this result be cached and the cache invalidated
        # when either the data or dq arrays change?
        if self.data is not None and self.data.ndim > 0 and self.dq is not None:
            if np.all(self.dq == 0):
                # All data good.
                return self.data
            else:
                self._data_mask = self._generate_mask(self.data, self.dq)
                self._data_fill_value = self._generate_fill(self.data,
                                                            self._data_fill)
                return ma.array(self.data, mask=self._data_mask,
                                fill_value=self._data_fill_value)
        else:
            return self.data

    @property
    def err_masked(self):
        # Generate the masked error array on the fly. This ensures the
        # masking is always up to date with the latest dq array.
        # TODO: Can this result be cached and the cache invalidated
        # when either the err or dq arrays change?
        if self.noerr:
            return None
        if self.err is not None and self.err.ndim > 0 and self.dq is not None:
            if np.all(self.dq == 0):
                # All data good.
                return self.err
            else:
                self._err_mask = self._generate_mask(self.err, self.dq)
                self._err_fill_value = self._generate_fill(self.err,
                                                           self._err_fill)
                return ma.array(self.err, mask=self._err_mask,
                                fill_value=self._err_fill_value)
        else:
            return self.err
        
    @property
    def data_filled(self):
        masked = self.data_masked
        if masked is not None and isinstance(masked, ma.masked_array):
            return masked.filled(self._data_fill_value)
        else:
            return self.data

    @property
    def err_filled(self):
        if self.noerr:
            return None
        masked = self.err_masked
        if masked is not None and isinstance(masked, ma.masked_array):
            return masked.filled(self._err_fill_value)
        else:
            return self.err


class HasDataErrAndGroups(HasDataErrAndDq):
    """
    
    An abstract class which overrides the data quality masking functions
    of HasDataErrAndDq for ramp data which contains PIXELDQ and RAMPDQ
    arrays instead of DQ. The DQ array for ramp data is read-only.

    """
    def __init__(self, data, err, noerr=False):
        super(HasDataErrAndGroups, self).__init__(data=data, err=err, dq=None,
                                                  noerr=noerr )

    def __add__(self, other):
        """
        
        Add a scalar, an array or another DataModel object to
        this MiriMeasuredModel object.
        
        """
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )
        
        if isinstance(other,(float,int)):
            # A scalar quantity is being added. Add to the SCI array but
            # leave the ERR and DQ arrays as they are.
            newobject.data = self.data + other
            if not self.noerr:
                newobject.err = self.err
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being added to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data + np.asarray(other)
            # Adding a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
             
        elif isinstance(other, DataModel):
            # Two data products are being added together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data + other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = \
                            self._combine_errors_quadrature(self.err,
                                                            other.err)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)
                    
                if hasattr(other, 'pixeldq') and self._isvalid(other.pixeldq):
                    newobject.pixeldq = self._combine_quality(self.pixeldq, other.pixeldq)
                if hasattr(other, 'groupdq') and self._isvalid(other.groupdq):
                    newobject.groupdq = self._combine_quality(self.groupdq, other.groupdq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot add " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)
        
        return newobject

    def __sub__(self, other):
        """
        
        Subtract a scalar, an array or another DataModel object from
        this MiriMeasuredModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being subtracted. Subtract from the SCI
            # array but leave the ERR and DQ arrays as they are.
            newobject.data = self.data - other
            if not self.noerr:
                newobject.err = self.err
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being subtracted to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data - np.asarray(other)
            # Adding a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, DataModel):
            # Two data products are being subtracted. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data - other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = \
                            self._combine_errors_quadrature(self.err, other.err)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)
                    
                if hasattr(other, 'pixeldq') and self._isvalid(other.pixeldq):
                    newobject.pixeldq = self._combine_quality(self.pixeldq, other.pixeldq)
                if hasattr(other, 'groupdq') and self._isvalid(other.groupdq):
                    newobject.groupdq = self._combine_quality(self.groupdq, other.groupdq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot subtract " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    def __mul__(self, other):
        """
        
        Multiply this MiriMeasuredModel object by a scalar, an array or
        another DataModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being multiplied. Multiply the SCI and ERR
            # arrays but leave the DQ array as it is.
            newobject.data = self.data * other
            if not self.noerr:
                newobject.err = self.err * other
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            newobject.data = self.data * np.asarray(other)
            # Multiplying a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, DataModel):
            # Two data products are being multiplied together. Ensure they
            # both have a valid primary data array.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data * other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = self._combine_errors_multiplicative( \
                                            self.err, other.err, self.data,
                                            other.data)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)

                if hasattr(other, 'pixeldq') and self._isvalid(other.pixeldq):
                    newobject.pixeldq = self._combine_quality(self.pixeldq, other.pixeldq)
                if hasattr(other, 'groupdq') and self._isvalid(other.groupdq):
                    newobject.groupdq = self._combine_quality(self.groupdq, other.groupdq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot multiply " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject
       
    def __truediv__(self, other):
        """
        
        Divide this MiriMeasuredModel object by a scalar, an array or
        another DataModel object.
        
        """  
        # Check this object is capable of mathematical operation.
        self._check_for_data()
        # Start with an empty version of the current object and clone
        # the metadata.
        newobject = self.__class__()
        newobject.update( self )

        if isinstance(other,(float,int)):
            # A scalar quantity is being divided. Divide the SCI and ERR
            # arrays but leave the DQ array as it is.
            # Trap a divide by zero..
            if np.abs(other) <= sys.float_info.epsilon:
                strg = "%s: Divide by scalar zero!" % self.__class__.__name__
                del newobject
                raise ValueError(strg)
            newobject.data = self.data / other
            if not self.noerr:
                newobject.err = self.err / other
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, (ma.masked_array,np.ndarray,list,tuple)):
            # A data array is being multiplied to this product. This should
            # work provided the two arrays are broadcastable.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            newobject.data = self.data / np.asarray(other)
            # Dividing by a plain data array erases the error information.
            if not self.noerr:
                newobject.err = np.zeros_like(self.err)
            newobject.pixeldq = self.pixeldq
            newobject.groupdq = self.groupdq
            
        elif isinstance(other, DataModel):
            # The data product is being divided by another. Ensure they
            # both have a valid primary data array.
            # NOTE: Any divide by zero operations will be trapped by numpy.
            if hasattr(other, 'data') and self._isvalid(other.data):
                newobject.data = self.data / other.data
                if not self.noerr:
                    if hasattr(other, 'err') and self._isvalid(other.err):
                        newobject.err = self._combine_errors_divisive( \
                                            self.err, other.err, self.data,
                                            other.data)
                    else:
                        # If only one error array is known, the combined error
                        # becomes unknown.
                        newobject.err = np.zeros_like(self.err)

                if hasattr(other, 'pixeldq') and self._isvalid(other.pixeldq):
                    newobject.pixeldq = self._combine_quality(self.pixeldq, other.pixeldq)
                if hasattr(other, 'groupdq') and self._isvalid(other.groupdq):
                    newobject.groupdq = self._combine_quality(self.groupdq, other.groupdq)
            else:
                raise TypeError("Both data products must contain a " + \
                                "primary data array.")            
        else:
            strg = "Cannot divide " + str(self.__class__.__name__)
            strg += " and " + str(other.__class__.__name__) + "objects."
            del newobject
            raise TypeError(strg)

        return newobject

    # From Python 3, division is the same as true division.
    def __div__(self, other):
        return self.__truediv__(other)

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the operations module.")

    import math

    # Check that dqflags has been imported properly
    print("Master data quality flags:")
    for flags in master_flags:
        print(flags)
    
    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dqtest = [[0,1,0],
              [0,1,1],
              [0,0,0]]
    dqtest2 = np.array([dqtest,dqtest,dqtest,dqtest])
    
    testobj = HasDataErrAndDq( data3x3, err3x3, dqtest2)
    newdq1 = testobj._shrink_dq( dqtest2 )
    print("\nData quality array:\n", dqtest2)
    print("has shrunk to:\n", newdq1)
    newdq2 = testobj._shrink_dq( newdq1 )
    print("and has shrunk again to:\n", newdq2)
    newdq3 = testobj._shrink_dq( newdq2 )
    print("and has shrunk finally to:\n", newdq3)
    del newdq1, newdq2, newdq3
    
    print("Testing combination and masking of data quality arrays")
    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dqtest = np.array([[0,1,0], [4,2,1], [0,3,0]])
    testobj = HasDataErrAndDq( data3x3, err3x3, dqtest2)
    mask1 = testobj._generate_mask(data3x3, dqtest, bitmask=None)
    print("\nGenerating mask from:\n", dqtest)
    print("with no bitmask gives:\n", str(mask1))
    mask2 = testobj._generate_mask(data3x3, dqtest, bitmask=1)
    print("\nGenerating mask from:\n", dqtest)
    print("with bitmask 1 gives:\n", str(mask2))
    mask3 = testobj._generate_mask(data3x3, dqtest, bitmask=3)
    print("\nGenerating mask from:\n", dqtest)
    print("with bitmask 3 gives:\n", str(mask3))
    del mask1, mask2, mask3
    
    # Testing error combination functions
    sq0 = 0.0
    sq1 = 1.0
    sq2 = math.sqrt(2.0)
    sq3 = math.sqrt(3.0)
    sq4 = 4.0
    sq5 = math.sqrt(5.0)
    sq6 = math.sqrt(6.0)
    sq7 = math.sqrt(7.0)
    sq8 = math.sqrt(8.0)
    sq9 = 3.0
    error1 = np.array([[sq0,sq1,sq2],[sq3,sq4,sq5],[sq7,sq8,sq9]])
    error2 = np.array([[sq9,sq8,sq7],[sq5,sq4,sq3],[sq2,sq1,sq0]])
    error0 = np.zeros_like(error1)
    
    print("\nCombining error array with itself:\n", error1)
    newerr = testobj._combine_errors_quadrature(error1, error1)
    print("by quadrature:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error0)
    newerr = testobj._combine_errors_quadrature(error1, error0)
    print("by quadrature:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    newerr = testobj._combine_errors_quadrature(error1, error2)
    print("by quadrature:\n", newerr)

    data0 = np.array([[0,0,0],[0,0,0],[0,0,0]])
    data1 = np.array([[1,1,1],[1,1,1],[1,1,1]])
    data2 = np.array([[2,2,2],[2,2,2],[2,2,2]])
    data_bad = np.array([[1,1,1],[1,0,1],[1,1,1]])

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted twice by:\n", data1)
    newerr = testobj._combine_errors_multiplicative(error1, error2, data1, data1)
    print("multiplicative:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted by:\n", data1)
    print("and:\n", data0)
    newerr = testobj._combine_errors_multiplicative(error1, error2, data1, data0)
    print("multiplicative:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted by:\n", data1)
    print("and:\n", data2)
    newerr = testobj._combine_errors_multiplicative(error1, error2, data1, data2)
    print("multiplicative:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted twice by:\n", data1)
    newerr = testobj._combine_errors_divisive(error1, error2, data1, data1)
    print("divisive:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted by:\n", data0)
    print("and:\n", data1)
    newerr = testobj._combine_errors_divisive(error1, error2, data0, data1)
    print("divisive:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted by:\n", data1)
    print("and:\n", data2)
    newerr = testobj._combine_errors_divisive(error1, error2, data1, data2)
    print("divisive:\n", newerr)

    print("\nCombining error array:\n", error1)
    print("with:\n", error2)
    print("weighted by:\n", data1)
    print("and:\n", data_bad)
    newerr = testobj._combine_errors_divisive(error1, error2, data1, data_bad)
    print("divisive:\n", newerr)

    print("Test finished.")
