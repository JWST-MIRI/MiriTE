#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""

Module util - Utility module containing functions shared by all the unit tests.

:History:

02 Sep 2013: Created.
11 Oct 2013: now assert_products_equal is able to recognize two arrays as
             equal even if they have nans (Ruyman Azzollini, DIAS).
29 Jun 2017: Explicitly test for identical NaN distributions using np.all()
             instead of using a subtraction.

@author: Steven Beard (UKATC), Michael Droettboom (STScI),
         Ruyman Azzollini (DIAS)

"""

from astropy.extern import six

import numpy as np
from numpy.testing import assert_array_equal


def assert_recarray_equal(a, b, msg=None):
    """
    
    Make sure that two numpy record arrays are identical.

    Numpy has a bug comparing recarrays, in that it assumes two
    recarrays with different byteorders to be different. This
    function works around that bug.
    
    :Parameters:
    
    a: numpy record array
        First record array
    b: numpy record array
        First record array
    msg: str, optional
        If specified, a message to be displayed if the assertion
        fails.
        
    Raises an exception if the arrays are not identical.
     
    """
    assert_array_equal(a, np.asarray(b, a.dtype), err_msg=msg)

def assert_products_equal(testobj, a, b, arrays='data', tables=''):
    """
    
    Test and ensure that two jwst_lib data products are identical.
    A jwst_lib data product is assumed to contain a set of
    attributes containing either a data array or a data table.
    
    Metadata isn't compared.
    
    :Parameters:
    
    testobj: TestCase object
        The test case being executed
    a: jwst_lib data product
        First data product
    b: jwst_lib data product
        Second data product
    arrays: str or list of str, optional
        A list of array attributes to be compared. By default,
        just the 'data' attribute will be tested.
        The arrays test will be skipped if an empty string or
        empty list is provided.
    tables: str or list of str, optional
        A list of table attributes to be compared. By default,
        no tables will be tested.
        The tables test will be skipped if an empty string or
        empty list is provided.
        
    Raises an exception if any of the arrays or tables are not identical.
    
    """
    # Make sure the attribute names are in list form.
    if arrays and isinstance(arrays, six.string_types):
        arrays = [arrays]
    if tables and isinstance(tables, six.string_types):
        tables = [tables]
        
    if arrays:
        for array_attribute in arrays:
            first = getattr(a, array_attribute)
            second = getattr(b, array_attribute)
            msg = "%s: First .%s array attribute is None" % \
                (a.__class__.__name__, array_attribute)
            testobj.assertIsNotNone( first, msg )
            msg = "%s: Second .%s array attribute is None" % \
                (b.__class__.__name__, array_attribute)
            testobj.assertIsNotNone( second, msg )
            msg = "%s: .%s 'nan' distributions in array attributes differ" % \
                (a.__class__.__name__, array_attribute)
            testobj.assertTrue(np.all(np.isnan(first) == np.isnan(second)), msg)
            
            msg = "%s: .%s array attributes differ" % \
                (a.__class__.__name__, array_attribute)
            testobj.assertTrue( np.allclose(np.nan_to_num(first),\
                np.nan_to_num(second) ), msg )
    if tables:
        for table_attribute in tables:
            first = getattr(a, table_attribute)
            second = getattr(b, table_attribute)
            msg = "%s: First .%s table attribute is None" % \
                (a.__class__.__name__, table_attribute)
            testobj.assertIsNotNone( first, msg )
            msg = "%s: Second .%s table attribute is None" % \
                (b.__class__.__name__, table_attribute)
            testobj.assertIsNotNone( second, msg )
            msg = "%s: .%s table attributes differ" % \
                (a.__class__.__name__, table_attribute)
            assert_recarray_equal(np.asarray(first), np.asarray(second),
                                  msg=msg)
