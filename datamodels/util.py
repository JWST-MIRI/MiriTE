#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Data product utilities.

Contains a common set of functions and utilities for testing and
processing the MIRI data models.

:History:

09 Oct 2013: Created
11 Oct 2013: now assert_products_equals is able to recognize two arrays as
             equal even if they have nans (Ruyman Azzollini, DIAS).
17 Oct 2013: Removed the warning about there being no tables.
11 Dec 2013: Added an open function, which attempts to open a data product
             with the correct data type.
13 Feb 2014: Corrected a lingering mistake where all occurences of "datatype"
             had not been converted to "astype".
17 Feb 2014: Delete a dataproduct before removing the file saved from it.
01 Oct 2014: Added the CDP conversion utilities. Use 'ANY' for CDP files
             where the detector settings do not matter.
01 Oct 2014: Added verify_metadata. Global constants from cdplib.py moved
             here to avoid a circular reference.
10 Oct 2014: Check for REFTYPE as well as TYPE in get_data_class.
17 Oct 2014: Update to CDP 2 to 3 conversion routine, adding functionality to 
             update values in VERSION, ORIGFILE, and FILENAME.
31 Oct 2014: Update in verify_cdp_file function to ignore the presence of an
             empty FIELD_DEF table, if a DQ_DEF is also present.
09 Dec 2014: Added USEAFTER and AUTHOR to CDP metadata check.
27 May 2015: Replaced pyfits with astropy.io.fits
25 Jun 2015: Updated the data model checks for CDP-4. Removed some obsolete code.
06 Jul 2015: Convert all data types to uppercase before looking up in CDP
             dictionary.
07 Oct 2015: Made exception catching Python 3 compatible.
04 Dec 2015: Worked around problem with HISTORY checking in verify_metadata.
             CDP files will now fail the metadata test if any HISTORY record
             is absent.
07 Dec 2015: Added verify_fits_file test.
08 Dec 2015: Corrected problem where get_data_class would fail when kwlist is
             completely empty. Added check for missing REFTYPE keywords to
             verify_fits_file. Added add_subarray_metadata function.
10 Dec 2015: PEDIGREE can be GROUND, DUMMY or SIMULATION.
             Corrected mix-up between TYPE and REFTYPE.
             Changed default subarray mode to 'GENERIC' in add_subarray.
20 May 2016: Corrected a bug which overwrote the kwlist and prevented
             a file of simulation data from being correctly opened.
22 Jun 2016: Corrected bug in add_subarray_metadata.
20 Jan 2017: Corrected depecation warning from astropy. Added a subarray
             definition to data shape integrity check.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
04 Oct 2017: Check for data type in the FILETYPE keyword.
29 Nov 2017: Removed obsolete find_fits_values function from CDP history
             check. REFTYPE trumps FILETYPE when assessing data type.
18 Jan 2018: Corrected a problem where the open() function closed an
             hdulist provided as input.
09 Apr 2018: Added 'N/A' to the lists of valid CDP options. 'N/A' is
             replaced by 'ANY' when looking up a CDP. Added FASTGRPAVG.
27 Apr 2018: Use str(e) to obtain an exception message rather than e.message.
             Replaced deprecated boolean - operator with logical_xor
             in assert_products_equal function. 
17 May 2018: Python 3: Converted dictionary keys return into a list.
29 Jun 2018: Global parameters moved to miri.parameters.
06 Sep 2018: Updated verify_metadata to check that FLAT, LINEARITY and PSF
             CDPs contain both FILTER, CHANNEL and BAND keywords.
04 Oct 2018: Added DISTORTION to the list of CDPs which must contain both
             FILTER, CHANNEL and BAND keywords.
12 Nov 2018: Added compulsory metadata for imager and MRS CDPs. Reworked the
             metadata warning message so it doesn't appear twice in the output.
29 Jan 2019: Check for a legacy DATAMODL keyword in verify_cdp_file.
08 Feb 2019: Added more CDP verification checks. Ensure all subarray coordinates
             match the MIRI definitions. Ensure a SKYFLAT uses PEDIGREE=DUMMY.
             Ensure that all data files use the correct dq_def field names.
11 Mar 2019: Added yet more CDP verification checks. Every CDP file is now
             analysed with fitsverify to ensure it meets FITS standards.
             The USAFTER keyword must exist and contain both a date and a time.
12 Mar 2019: Removed use of astropy.extern.six (since Python 2 no longer used).

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""

import os
import copy
import warnings
import subprocess

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

import astropy.io.fits as pyfits

# Import the JWST data models and the CDP and SIM dictionaries
import jwst.datamodels
from miri.parameters import SUBARRAY, CDP_METADATA, CDP_METADATA_IMAGER, \
    CDP_METADATA_MRS, CDP_SUBARRAY, CDP_HISTORY
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.cdp import CDP_DICT
from miri.datamodels.sim import SIM_DICT

# -----------------------------------------------------------------------------
#
# Common utility functions
#
def read_fits_header( fitsobject, keyword, hduname='' ):
    """
    
    Returns the value of a named keyword contained in a FITS header.
    
    :Parameters:
    
    fitsobject: str or pyfits.HDUList object
        Either: The name of a FITS file (string);
        Or: A file descriptor from an already-opened FITS file;
        Or: A pyfits.HDUList object from an already-opened file.
    keyword: str
        The name of the keyword to be read.
    hduname: str (optional)
        The name of the HDU containing the header.
        Defaults to the primary HDU.
        Not relevant if fitsobject is already an pyfits.HDUList object.
    
    :Returns:
    
    value: object
        The value contained in the file header.
        None if the header keyword could not be found
    
    """
    value = None
    hdulist = None
    file_open = False

    try:
        if isinstance(fitsobject, bytes):
            fitsobject = fitsobject.decode('utf-8')
        if isinstance(fitsobject, str) or hasattr(fitsobject, "read"):
            hdulist = pyfits.open( fitsobject )
            file_open = True
        elif isinstance(fitsobject, pyfits.HDUList):
            hdulist = fitsobject
        if hdulist is not None:
            if hduname:
                header = hdulist[hduname].header
            else:
                header = hdulist[0].header
            if keyword in header:
                value = header[keyword]
    except Exception as e:
        strg = "Failed to open FITS object, \'%s\'\n" % str(fitsobject)
        strg += "  %s: %s" % (e.__class__.__name__, str(e))
        raise IOError(strg)
    finally:
        if file_open and hdulist is not None:
            hdulist.close()
            del hdulist
    return value

def get_data_class( kwlist, dictionary=CDP_DICT ):
    """
    
    Return a data model class which matches the given list of
    keywords. The class is looked up in a hierarchical dictionary.
    
    :Parameters:
    
    kwlist: str or list of str
        A keyword, or a list of keywords, identifying the data type
        within the given dictionary. The first keyword is the data
        type identifier. More keywords (such as detector name and filter)
        can be provided to narrow down the search
    dictionary, dict (optional)
        The dictionary of data models to be used. If not specified, the
        CDP dictionary will be used.
        
    :Returns:
    
    dataclass: Class
        The class with which to open the defined data model.
        None is returned if no class could be matched.
    
    """
#     print("get_data_class: looking up", kwlist)
    # Make sure kwlist is a list.
    if isinstance(kwlist, bytes):
        kwlist = kwlist.decode('utf-8')
    if isinstance(kwlist, str):
        kwlist = [kwlist]
    if kwlist is None or len(kwlist) < 1:
        # No keywords at all.
        return None
    # Convert the data type keyword to uppercase, to save duplicate entries
    # in the dictionary ("Distortion" is the same as "DISTORTION").
    kwlist[0] = kwlist[0].upper()
    # Reverse the list so it may be popped from the back.
    kwlist.reverse()
        
    def lookup( klist, kdict ):
        if len(klist) > 0:
            key = klist.pop()
        else:
            key = 'ANY'
        # 'N/A' is an alias for 'ANY'
        if key == 'N/A':
            key = 'ANY'            
        if key in kdict:
            thing = kdict[key]
            if isinstance(thing, dict):
                return lookup( klist, thing )
            else:
                return thing
        elif 'ANY' in kdict:
            return kdict['ANY']
        else:
            return None
    
    # Look up the class by a recursive search through the dictionary.
    dataclass = lookup( kwlist, dictionary )
    return dataclass

# -----------------------------------------------------------------------------
#
# Universal data model opening function
#
def open( init=None, astype=None):
    """
    
    Creates a MIRI Data Product object from the given initializer using
    the class identified by its metadata (c.f. the jwst.datamodels.open
    function).
    
    Examples:
    
    datamodel = open( "myfile.fits" )
    
    datamodel = open( "myfile.fits", astype="PSF" )
    
    datamodel = open( "myfile.fits", astype=MiriImagerPointSpreadFunction )
    
    datamodel = open( (3,4) )

    :Parameters:

    init: file path, file object, pyfits.HDUList, shape tuple, numpy array, None
        The object from which to open a MIRI data product. Normally,
        this will be a string containing the name of a FITS file, or a
        pyfits.HDUList object, containing metadata. Other kinds of object
        are passed to the jwst.datamodels.open function, unless a specific
        data type is specified. The full list of possibilities is as follows:

        * file path: Initialize from the given file
        * readable file object: Initialize from the given file object
        * pyfits.HDUList: Initialize from the given `pyfits.HDUList`
        * A numpy array: A new model with the data array initialized
          to what was passed in.
        * shape tuple: Initialize with empty data of the given shape.
        * None: A default data model with no shape.
          
    astype: str or class, optional
        If specified, overrides the data type specified in the metadata.
        If there is no metadata, this parameter can be used to specify
        the name of the data model to be created.
        The parameter is normally a string describing the data product to be
        created.
        Examples 'Bad', 'Dark', 'PixelFlat', 'Lin', 'RSRF', 'PSF', etc....
        See the MIRI wiki for a complete list.
        The data type can also be specified by providing a class explicitly
        in this parameter.

    :Returns:

    mirimodel: MIRI data model instance
        A new data model object.
    
    """
    # Initialise defaults.
    hdulist = None
    preserve_hdulist = False
    kwlist = []
    if astype is None:
        datatype = ''
        mirimodel = None
    elif isinstance(astype, str):
        datatype = astype
        mirimodel = None
    elif isinstance(astype, bytes):
        datatype = astype.decode('utf-8')
        mirimodel = None
    else:
        datatype = ''
        mirimodel = astype

    # Check whether a class has been specified explicitly
    if mirimodel is None:
        # No class has been specified.
        # If the initialiser is a FITS-type object, the data type may
        # be distinguished using the file header.
        try:
            if isinstance(init, str) or hasattr(init, "read"):
                hdulist = pyfits.open(init)
            elif isinstance(init, bytes):
                hdulist = pyfits.open(init.decode('utf-8'))
            elif isinstance(init, pyfits.HDUList):
                hdulist = init
                preserve_hdulist = True

            if hdulist is not None:
                if len(hdulist) == 1:
                    # MIRI and JWST data models have more than 1 HDU.
                    # This may be a FITSWriter data product.
                    strg = "\n***The data file contains only a PRIMARY HDU. "
                    strg += "The input may be a FITSWriter data product. "
                    warnings.warn(strg)
                
                header = hdulist[0].header
                header_keys = list(header.keys())
#                 for key in header_keys:
#                     if key:
#                         print("%s = %s" % (key, str(header[key])))
                if datatype:
                    # The data type in the file header is overriden.
                    kwlist.append(datatype)
                elif 'REFTYPE' in header or 'REFTYPE' in header_keys:
                    # There is a reference data type keyword in the header.
                    kwlist.append(header['REFTYPE'])
                elif 'FILETYPE' in header or 'FILETYPE' in header_keys:
                    # There is a new data type keyword in the header.
                    kwlist.append(header['FILETYPE'])
                elif 'DATAMODL' in header or 'DATAMODL' in header_keys:
                    # There is an old data type keyword in the header.
                    kwlist.append(header['DATAMODL'])
                elif 'TYPE' in header or 'TYPE' in header_keys:
                    # There is an old data type keyword in the header.
                    kwlist.append(header['TYPE'])
                if 'DETECTOR' in header or 'DETECTOR' in header_keys:
                    # There is a detector keyword in the header.
                    detector = header['DETECTOR']
                    if detector:
                        kwlist.append(detector)
                if 'FILTER' in header or 'FILTER' in header_keys:
                    # There is a filter keyword in the header.
                    mirifilter = header['FILTER']
                    if mirifilter:
                        kwlist.append(mirifilter)
#                 elif 'CHANNEL' in header or 'CHANNEL' in header_keys:
#                     # There is a channel keyword in the header.
#                     mirichannel = header['CHANNEL']
#                     if mirichannel:
#                         kwlist.append(mirichannel)
            elif datatype:
                # There is no metadata, but the data type is explicitly
                # specified.
                kwlist.append(datatype)
            
        except Exception as e:
            strg = "Failed to open FITS object, \'%s\'\n" % str(init)
            strg += "  %s: %s" % (e.__class__.__name__, str(e))
            raise IOError(strg)
        finally:
            if hdulist is not None and not preserve_hdulist:
                hdulist.close()
                del hdulist

        # Attempt to convert the keyword list to a MIRI data model class
        mirimodel = get_data_class(copy.copy(kwlist), dictionary=CDP_DICT)
        if mirimodel is None:
            # Not a CDP. Try the simulation data products.
            mirimodel = get_data_class(copy.copy(kwlist), dictionary=SIM_DICT)
        if mirimodel is None:    
            # Unknown initialisers or data types are interpreted by jwst.datamodels
            strg = "\n***Data type could not be determined from metadata. "
            strg += "Opening as a plain jwst.datamodels data model."
            warnings.warn(strg)
            mirimodel = jwst.datamodels.open( init )
            return mirimodel

    return mirimodel(init)

# -----------------------------------------------------------------------------
#
# Data model testing support functions
#
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
        
    Raises an exception if the recarrays are not identical.
     
    """
    anames = a.dtype.names
    bnames = b.dtype.names
    if not np.alltrue(anames == bnames):
        strg = "Sets of columns in recarrays are different\n"
        raise TypeError(strg)
    for item in anames:
        if np.issubsctype(a[item],str):
            assert_array_equal(a[item], b[item], err_msg=msg)
        else:
            assert_array_equal(np.nan_to_num(a[item]), \
                np.nan_to_num(b[item]), err_msg=msg)

def assert_products_equal(a, b, arrays='data', tables=''):
    """
    
    Test and ensure that two jwst.datamodels data products are identical.
    A jwst.datamodels data product is assumed to contain a set of
    attributes containing either a data array or a data table.
    
    Raises an exception if any of the arrays or tables are not identical.

    Metadata isn't compared.
    
    :Parameters:
    
    a: jwst.datamodels data product
        First data product
    b: jwst.datamodels data product
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
    
    """
    # Make sure the attribute names are in list form.
    if arrays and isinstance(arrays, str):
        arrays = [arrays]
    if tables and isinstance(tables, str):
        tables = [tables]
        
    if arrays:
        for array_attribute in arrays:
            first = getattr(a, array_attribute)
            second = getattr(b, array_attribute)

            if first is not None and second is not None:
#                 discr_nan = (np.isnan(first) - np.isnan(second)).sum() == 0
                discr_nan = not np.any(np.logical_xor(np.isnan(first), np.isnan(second)))
                discr_all = np.allclose(np.nan_to_num(first),
                                        np.nan_to_num(second))

                if not (discr_nan and discr_all):
                    strg = "Data objects \'%s\' and \'%s\' have failed the test." % \
                        (a.__class__.__name__, b.__class__.__name__)
                    strg += "\n  \'%s\' array attribute has changed." % \
                        array_attribute
                    strg += "\n" + str(first)
                    strg += "\n   c.f.\n" + str(second)
                    raise TypeError(strg)   
            else:
                strg = "Data objects \'%s\' and \'%s\' have failed the test." % \
                    (a.__class__.__name__, b.__class__.__name__)
                if first is None and second is None:
                    strg += "\n  Both array attributes \'%s\' are None" % \
                        array_attribute
                elif first is None:
                    strg += "\n  The first array attribute \'%s\' is None" % \
                        array_attribute
                else:
                    strg += "\n  The second array attribute \'%s\' is None" % \
                        array_attribute              
                raise TypeError(strg)

    if tables:
        for table_attribute in tables:
            first = getattr(a, table_attribute)
            second = getattr(b, table_attribute)
            
            if first is not None and second is not None:
                try:
                    assert_recarray_equal(first,second)
                except Exception as e:
                    strg = "Data objects \'%s\' and \'%s\' have failed the test." % \
                        (a.__class__.__name__, b.__class__.__name__)
                    strg += "\n  \'%s\' table attribute has changed \n" % \
                        table_attribute
                    strg += "(%s)" % str(e)
                    strg += "\n" + str(first)
                    strg += "\n\t c.f.\n" + str(second)
                    raise TypeError(strg)   
            else:
                strg = "Data objects \'%s\' and \'%s\' have failed the test." % \
                        (a.__class__.__name__, b.__class__.__name__)
                if first is None and second is None:
                    strg += "\n  Both table attributes \'%s\' are None" % \
                        table_attribute
                elif first is None:
                    strg += "\n  The first table attribute \'%s\' is None" % \
                        table_attribute
                else:
                    strg += "\n  The second table attribute \'%s\' is None" % \
                        table_attribute
                raise TypeError(strg)   

def add_subarray_metadata(datamodel, subname=''):
    """
    
    Add missing subarray keywords to the metadata of a data model.
        
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be updated.
    subname: str (optional)
        The name of the subarray to be applied. It is only necessary to
        specify a new subarray if 'FULL' is not appropriate, or if the
        data model specifies an incorrect subarray mode.
        By default, the subarray mode is inferred from the data content.
        
    :Returns:
    
    subarray_modified: bool
        Returns True if the subarray metadata have been modified.
    
    """
    # Check subarray information.
    subarray_modified = False
    if hasattr(datamodel, 'meta') and hasattr(datamodel.meta, 'subarray'):
        if not subname:
            subname = 'GENERIC'
        if hasattr(datamodel.meta.subarray, 'name'):
            if datamodel.meta.subarray.name:
                subname = datamodel.meta.subarray.name
        # Remember the shape of the primary data array.
        paname = datamodel.get_primary_array_name()
        if paname == 'data':
            datashape = datamodel.data.shape
        elif paname == 'dq':
            datashape = datamodel.dq.shape
        else:
            parray = getattr(datamodel, paname)
            datashape = getattr(parray, 'shape')
        # Only add new subarray information if the existing information
        # is not defined.
        if datamodel.meta.subarray.xstart is None or \
           datamodel.meta.subarray.xstart <= 0 or \
           datamodel.meta.subarray.ystart is None or \
           datamodel.meta.subarray.ystart <= 0 or \
           datamodel.meta.subarray.xsize is None or \
           datamodel.meta.subarray.xsize <= 0 or \
           datamodel.meta.subarray.ysize is None or \
           datamodel.meta.subarray.ysize <= 0:
            # Subarray information is not present.
            if subname != 'GENERIC':
                datamodel.set_subarray_metadata( subname )
            else:
                subtuple = (1, 1, datashape[1], datashape[0])
                datamodel.set_subarray_metadata( subtuple )
            strg = "Subarray metadata updated to %s: " % \
                datamodel.meta.subarray.name
            strg += "%d " % datamodel.meta.subarray.xstart
            strg += "%d " % datamodel.meta.subarray.ystart
            strg += "%d " % datamodel.meta.subarray.xsize
            strg += "%d " % datamodel.meta.subarray.ysize
            print( strg )
            subarray_modified = True
            
        # Check the declared subarray size against the actual size of
        # the data array.
        if datamodel.meta.subarray.xsize != datashape[1] or \
           datamodel.meta.subarray.ysize != datashape[0]:
            strg = "\nDeclared size of subarray (%d,%d) " % \
                (datamodel.meta.subarray.xsize,
                 datamodel.meta.subarray.ysize)
            strg += "is not the same as the primary array size (%d,%d)." % \
                (datashape[1], datashape[0])
            warnings.warn(strg)
            
        # Ensure the fast and slow axes are declared properly
        if datamodel.meta.subarray.fastaxis is None or \
           datamodel.meta.subarray.fastaxis <= 0 or \
           datamodel.meta.subarray.slowaxis is None or \
           datamodel.meta.subarray.slowaxis <= 0:
            datamodel.meta.subarray.fastaxis = 1
            datamodel.meta.subarray.slowaxis = 2
            strg = "Fast/slow axis metadata updated to: %d/%d" % \
                (datamodel.meta.subarray.fastaxis,
                 datamodel.meta.subarray.slowaxis)
            print( strg )
            subarray_modified = True
    else:
        strg = "Subarray metadata attributes missing from data model %s" % \
            datamodel.__class__.__name_
        raise TypeError(strg)
    return subarray_modified

def verify_subarray_metadata(datamodel):
    """
    
    Verify that the subarray metadata contained in the given calibration
    data model is self-consistent.
    
    The function assumes that a basic metadata check verifying the existence
    of subarray metadata has already been made.
    
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be verified.
        
    :Returns:
    
    failure_string: str
        A string indicating problems with the metadata.
        If the metadata passes the test, an empty string is returned.
    
    """
    failure_string = ''
    
    if hasattr(datamodel.meta, 'subarray'):
        if hasattr(datamodel.meta.subarray, 'name') and \
           hasattr(datamodel.meta.subarray, 'xstart') and \
           hasattr(datamodel.meta.subarray, 'ystart') and \
           hasattr(datamodel.meta.subarray, 'xsize') and \
           hasattr(datamodel.meta.subarray, 'ysize'):
            if datamodel.meta.subarray.name and \
               str(datamodel.meta.subarray.name) in list(SUBARRAY.keys()):
                subname = str(datamodel.meta.subarray.name)
                if subname != 'GENERIC':
                    xstart = int(datamodel.meta.subarray.xstart)
                    ystart = int(datamodel.meta.subarray.ystart)
                    xsize = int(datamodel.meta.subarray.xsize)
                    ysize = int(datamodel.meta.subarray.ysize)
    
                    (expected_xstart, expected_ystart,
                     expected_xsize, expected_ysize) = SUBARRAY[subname]
                     
                    if (expected_xstart != xstart) or (expected_ystart != ystart) or \
                       (expected_xsize != xsize) or (expected_ysize != ysize):
                        failure_string = "  Wrong subarray coordinates for %s:" % subname
                        failure_string += " Expected [%d, %d, %d, %d];" % \
                            (expected_xstart, expected_ystart,
                             expected_xsize, expected_ysize)
                        failure_string += " Actual [%d, %d, %d, %d]" % \
                            (xstart, ystart, xsize, ysize)
            else:
                failure_string = "  Invalid subarray name (%s)" % \
                    str(datamodel.meta.subarray.name)
        else:
            failure_string = "  No subarray metadata!"
    else:
        failure_string = "  No subarray metadata!"
    return failure_string
    
def verify_metadata(datamodel):
    """
    
    Verify that a given calibration data model contains the compulsory
    keywords.
        
    :Parameters:
    
    datamodel: MiriDataModel
        The calibration data model whose metadata is to be verified.
        
    :Returns:
    
    failure_string: str
        A string indicating problems with the metadata.
        If the metadata passes the test, an empty string is returned.
    
    """
    failure_string = ''
    
    if not hasattr(datamodel, 'meta'):
        failure_string = "  Data model object %s has no metadata." % \
            datamodel.__class__.__name__
        return failure_string
    if not hasattr(datamodel.meta, 'reftype'):
        failure_string = "  Data model object %s has no reference file type keyword (REFTYPE)." % \
            datamodel.__class__.__name__
        return failure_string
    if not hasattr(datamodel.meta, 'useafter'):
        failure_string = "  Data model object %s has no reference file type keyword (USEAFTER)." % \
            datamodel.__class__.__name__
        return failure_string

    if not datamodel.meta.reftype:
        failure_string += "  Empty reference file type keyword (REFTYPE).\n"
    elif datamodel.meta.reftype not in CDP_DICT:
        failure_string += "  Reference file type of \'%s\' is not recognised.\n" % \
            datamodel.meta.reftype
            
    # A reference file must contain a USEAFTER keyword containing a
    # date and a time.

    if not datamodel.meta.useafter:
        failure_string += "  Empty reference file use after keyword (USEAFTER).\n"
    elif 'T' not in str(datamodel.meta.useafter):
        failure_string += "  USEAFTER string (%s) does " % datamodel.meta.useafter
        failure_string += "not include a time component (e.g. T00:00:00).\n"    
    
    # Check for compulsory keywords and default values.
    # Check all the compulsory metadata. Add extra subarray keywords if
    # a subarray is defined and is not 'GENERIC'. Make a copy so that
    # future tests made by this module are not affected.
    check_list = copy.copy(CDP_METADATA)
    subarray = datamodel.get_fits_keyword('SUBARRAY')
    if subarray is not None and 'GENERIC' not in subarray:
        check_list += CDP_SUBARRAY
        
    # The FLAT, LINEARITY, PSF and DISTORTION data must contain both
    # imager and MRS keywords
    if 'FLAT' in datamodel.meta.reftype or \
       datamodel.meta.reftype == 'LINEARITY' or \
       datamodel.meta.reftype == 'PSF' or \
       datamodel.meta.reftype == 'DISTORTION':
        check_list += CDP_METADATA_IMAGER + CDP_METADATA_MRS
    else:
        # Otherwise the remaining compulsory keywords depend on the detector.
        detector = datamodel.meta.instrument.detector
        if detector == 'MIRIMAGE':
            check_list += CDP_METADATA_IMAGER
        elif detector == 'MIRIFUSHORT' or detector == 'MIRIFULONG':
            check_list += CDP_METADATA_MRS
            
    for (name,test_values) in check_list:
#         print("Looking for %s with allowed values %s" % (name, str(test_values)))
        try:
            value = datamodel.get_fits_keyword(name)
#             print("%s = %s" % (name, str(value)))
            # Check for a null value, empty string or white space.
            if value is not None:
                if len(str(value).strip()) > 0:
                    if isinstance(test_values,(tuple,list)):
                        # Multiple possible values
                        if len(test_values) > 0:
                            if value not in test_values:
                                failure_string += "  Value of \'%s\' is " % name
                                failure_string += "%s instead of one of %s.\n" % \
                                    (str(value),str(test_values))
                    else:
                        # Single possible value
                        if value != test_values:
                            failure_string += "  Value of \'%s\' is " % name
                            failure_string += "%s instead of %s.\n" % \
                                (str(value),str(test_values))
                else:
                    failure_string += "  Compulsory metadata keyword "
                    failure_string += "\'%s\' is blank.\n" % name
            else:
                failure_string += "  Compulsory metadata keyword "
                failure_string += "\'%s\' is not defined.\n" % name
        except KeyError:
            failure_string += "  Compulsory metadata keyword "
            failure_string += "\'%s\' not found in the data model schema.\n" % name
            
    # Check for HISTORY records.
#     hlist = datamodel.find_fits_values('HISTORY')
    hlist = datamodel.get_history()
    if hlist:
        history_checked = len(CDP_HISTORY) * [False]
        for history in hlist:
            for ii in range(0,len(CDP_HISTORY)):
                if CDP_HISTORY[ii] in history:
                    history_checked[ii] = True
        not_found = []
        for ii in range(0,len(CDP_HISTORY)):
            if not history_checked[ii]:
                not_found.append("\'" + CDP_HISTORY[ii] + "\'")
        if len(not_found) > 0:
            strg = '  HISTORY records %s not found' % ','.join(not_found)
            failure_string += strg
    else:
        failure_string += '  No HISTORY records found.'

    # Check for some exceptional cases.
    # A SKYFLAT must, for the time being, have PEDIGREE defined as 'DUMMY'.
    if datamodel.meta.reftype == 'SKYFLAT':
        if hasattr(datamodel.meta, 'pedigree') and datamodel.meta.pedigree:
            if datamodel.meta.pedigree != 'DUMMY':
                failure_string = "  SKYFLAT file contains PEDIGREE=%s." % \
                    str(datamodel.meta.pedigree)
                failure_string += " PEDIGREE=DUMMY expected."
        else:
            failure_string = "  No PEDIGREE metadata found."
    return failure_string

def verify_cdp_file(filename, datatype=None, overwrite=False, keepfile=False):
    """
    
    Verify that a given file contains a valid MIRI calibration data
    product. The function can also be used to ensure a new software
    release has not affected the validity of a data product.
    
    Throws an exception if the CDP is not valid.
    
    NOTE: Part of the test includes writing the data product to a
    new file, which has the same name as the original file with
    "_copy.fits" appended to the name. The function will fail if
    a file of that name already exists, unless called with overwrite=True.
        
    :Parameters:
    
    filename: str
        The name of the FITS file containing the CDP.
    datatype: str (option)
        The data type expected. If not given, it will be obtained from
        the REFTYPE or TYPE keyword in the file header.
    overwrite: boolean (optional)
        Set True to overwrite an existing file when making a
        temporary copy (not the original file itself).
    keepfile: boolean (optional)
        Set to True to keep the temporary file instead of
        removing it.
        
    :Returns:
    
    datatype: str
        The data type as tested. This may have been obtained from the
        file header.
    
    """
    assert isinstance(filename, str)
    outputfile = filename + "_copy.fits"
    # Does the file exist?
    if not os.path.isfile(filename):
        strg = "verify_cdp_file: Specified file \'%s\' does not exist!" % filename
        raise IOError(strg)

    # Read the data model either using a class derived from the data type
    # or using the class specified in the datatype parameter.
    # Test 1 - the CDP must be readable.
    with open( init=filename, astype=datatype ) as datamodel:
        # Verify the metadata
        metadata_failure = False
        metadata_strg = ''
        failure_string = verify_metadata(datamodel)
        subarray_string = verify_subarray_metadata(datamodel)
        if failure_string or subarray_string:
            metadata_strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            metadata_strg += "\n  Incorrect metadata."
            if failure_string:
                metadata_strg += "\n" + failure_string
            if subarray_string:
                 metadata_strg += "\n" + subarray_string
               
            # For now only set a flag, so the data content may also be checked.
            # The flag will trigger a warning or an exception later.
            metadata_failure = True
                
        if isinstance(datamodel, MiriDataModel):
            dataarrays = datamodel.list_data_arrays()
            datatables = datamodel.list_data_tables()
#             print("List of data arrays:", dataarrays)
#             print("List of data tables:", datatables)
        else:
            if metadata_failure:
                warnings.warn(metadata_strg)   
            strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            strg += "\n  Not a recognised MIRI data model."
            raise TypeError(strg)
        
        # The CDP must contain either at least one data array or
        # at least one data table.
        if len(dataarrays) <= 0 and len(datatables) <= 0:
            if metadata_failure:
                warnings.warn(metadata_strg)   
            strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            strg += "\n  No data arrays and no data tables."
            raise TypeError(strg)

        # Warn about obsolete data tables (CDP-2 to CDP-3 only).
        if 'FIELD_DEF' in datatables or 'field_def' in datatables:
            strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            strg += "\n  It contains an obsolete FIELD_DEF table."
            raise TypeError(strg)
        
        if 'GROUP_DEF' in datatables or 'group_def' in datatables:
            strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            strg += "\n  It contains an obsolete GROUP_DEF table."
            raise TypeError(strg)

        # Update the data type
        if hasattr(datamodel.meta, 'reftype'):
            datatype = datamodel.meta.reftype
        elif hasattr(datamodel.meta, 'datatype'):
            datatype = datamodel.meta.datatype
            
        # Check the data model type. This keyword is compulsory.
        if hasattr(datamodel.meta, 'model_type'):
            model_type = str(datamodel.meta.model_type)
            class_name = str(datamodel.__class__.__name__)
            if not model_type:
                if metadata_failure:
                    warnings.warn(metadata_strg)   
                strg = "Data object \'%s\' has failed the test." % \
                    (datamodel.__class__.__name__)
                strg += "\n  It contains an empty DATAMODL keyword."
                raise TypeError(strg)
            else:
                # If a data model has no STScI equivalent, the DATAMODL keyword
                # is set to the MIRI class name.
                if model_type.strip().upper() == class_name.strip().upper():
                    strg = "\n***WARNING: Data object \'%s\' contains " % \
                        class_name
                    strg += "a non-STScI DATAMODL=\'%s\' in the metadata." % \
                        model_type
                    warnings.warn(strg)
        else:
            if metadata_failure:
                warnings.warn(metadata_strg)   
            strg = "Data object \'%s\' has failed the test." % \
                (datamodel.__class__.__name__)
            strg += "\n  No data model type (DATAMODL) defined."
            raise TypeError(strg)                
    
        # Test 2 - all data arrays and tables must be readable and not empty.
#         print("Test 2")
        for dataarray in dataarrays:
            isempty = False
            if datamodel[dataarray] is None:
                isempty = True
            else: #  or len(datamodel[dataarray]) < 1
                if datamodel[dataarray].size < 1:
                    isempty = True
            
            if isempty:
                if metadata_failure:
                    warnings.warn(metadata_strg)   
                strg = "Data object \'%s\' has failed the test." % \
                    (datamodel.__class__.__name__)
                strg += "\n  Data array %s is empty." % dataarray             
                raise TypeError(strg)
        
        for datatable in datatables:
            isempty = False
            if datamodel[datatable] is None:
                isempty = True
            else:
                if len(datamodel[datatable]) < 1:
                    isempty=True                
            if isempty:
                if metadata_failure:
                    warnings.warn(metadata_strg)   
                strg = "Data object \'%s\' has failed the test." % \
                    (datamodel.__class__.__name__)
                strg += "\n  Data table %s is empty." % datatable
                raise TypeError(strg)
                        
        # Check that the primary array of subarray data has the correct size.
        # The check is made when the subarray is not FULL and the primary
        # array exists and contains 2-D data.
        if hasattr(datamodel.meta, 'subarray'):
            if hasattr(datamodel.meta.subarray, 'name') and \
               hasattr(datamodel.meta.subarray, 'xsize') and \
               hasattr(datamodel.meta.subarray, 'ysize'):
                if datamodel.meta.subarray.name and \
                   datamodel.meta.subarray.name != 'FULL':
                    xsize = datamodel.meta.subarray.xsize
                    ysize = datamodel.meta.subarray.ysize
                    primary = datamodel.get_primary_array_name()
                    if (xsize is not None) and (ysize is not None) and \
                        hasattr(datamodel, primary):
                        data = getattr(datamodel, primary)
                        if data.ndim > 1 and \
                           (data.shape[-1] != xsize or data.shape[-2] != ysize):
                            if metadata_failure:
                                warnings.warn(metadata_strg)   
                            strg = "Data object \'%s\' has failed the test." % \
                                (datamodel.__class__.__name__)
                            strg += "\n  Primary data array shape (%d x %d)" % \
                                (data.shape[-1], data.shape[-2])
                            strg += " incompatible with subarray shape (%d x %d)" % \
                                (xsize, ysize)
                            raise TypeError(strg)
                    
        # Bail out here if there has been a metadata problem.
        if metadata_failure:
            raise TypeError(metadata_strg)

        # Test 3 - It must be possible to save the product to a new file
#         print("Test 3")
        try:
            datamodel.save( outputfile, overwrite=overwrite )
        except Exception as e:
            
            strg = str(e) + "\nNot possible to save data product to a file."
            if not keepfile:
                # Delete the data model to close the temporary file.
                del datamodel
                # Attempt to remove the temporary file.
                # Do not throw an exception if this fails.
#                 print("Attempt to remove old file after 3")
                try:
                    os.remove(outputfile)
                except Exception:
                    warnings.warn("\n***Could not remove temporary file: " + \
                                  outputfile)
            raise TypeError(strg)
        
        # Test 4 - The input and copied data products must be equal
#         print("Test 4")
        with open( init=outputfile, astype=datatype ) as newmodel:
             
            assert_products_equal(datamodel, newmodel, arrays=dataarrays,
                                  tables=datatables)
            del newmodel
        del datamodel
        
        if not keepfile:
            # Attempt to remove the temporary file.
            # Do not throw an exception if this fails.
#             print("Attempt to remove old file after 4")
            try:
                os.remove(outputfile)
            except Exception:
                warnings.warn("\n***Could not remove temporary file: " + \
                              outputfile)
        return datatype

def verify_fits_file(filename, cdp_checks=False, fitsverify_checks=True):
    """
    
    Verify that a given FITS file looks vaguely like a MIRI/JWST
    data model. The file must be structured so that the primary HDU
    is empty (apart from a valid header) and information is contained
    in extension HDUs.
    
    Throws an exception if the file is not valid.
    
    :Parameters:
    
    filename: str
        The name of the FITS file containing the CDP.
    cdp_checks: bool
        Set True for additional CDP checks. Defaults to False.
    fitsverify_checks: bool
        Set True for an additional fitsverify check. Defaults to True.
        
    :Returns:
    
    None.
    
    """
    assert isinstance(filename, str)
    # Does the file exist?
    if not os.path.isfile(filename):
        strg = "verify_fits_file: Specified file \'%s\' does not exist!" % filename
        raise IOError(strg)
    
    # First verify that the file obeys the rules of the FITS standard.
    # This step is only executed if the host system recognises the
    # "fitsverify" shell command.
    if fitsverify_checks:
        fitsverify_available = False
        try:
            if subprocess.call(["which", "fitsverify"]) == 0:
                fitsverify_available = True
        except:
            pass
        if fitsverify_available:
            nerrors = subprocess.call(["fitsverify", "-q", filename])
            if nerrors > 0:
                strg = "FITS standard violation (%d errors)." % nerrors
                strg += " Rerun \'fitsverify %s\' for details." % filename
                raise TypeError(strg)

    failure_strg = ''
    with pyfits.open( filename ) as hdulist:
        if len(hdulist) < 1:
            failure_strg += "Only one HDU. "
        if hdulist[0].size > 0:
            failure_strg += "Primary HDU is not empty. "
        if len(hdulist[0].header) < 1:
            failure_strg += "No primary FITS header. "
        else:
            if 'SIMPLE' not in hdulist[0].header or \
               'NAXIS' not in hdulist[0].header:
                failure_strg += "Invalid primary FITS header. "   
            if 'TELESCOP' not in hdulist[0].header or \
               'INSTRUME' not in hdulist[0].header:
                failure_strg += "No TELESCOP or INSTRUME keywords in primary FITS header. "
            if cdp_checks:   
                if 'REFTYPE' not in hdulist[0].header:
                    failure_strg +=  "No REFTYPE keyword in primary header. "
                    
                # Check whether that common data tables have correctly
                # named columns.
                if 'DQ_DEF' in hdulist:
                    dq_def_header = hdulist['DQ_DEF'].header
                    if 'TFIELDS' in dq_def_header:
                        tfields = int(dq_def_header['TFIELDS'])
                        if tfields == 4:
                            fieldnames = []
                            testnames = ['BIT', 'VALUE', 'NAME', 'DESCRIPTION']
                            names_ok = True
                            for index in range(0, tfields):
                                index1 = index + 1
                                ttypename = 'TTYPE%d' % index1
                                ttype = dq_def_header[ttypename]
                                fieldnames.append(ttype)
                                if ttype.strip() != testnames[index].strip():
                                    names_ok = False
                            if not names_ok:
                                strg = "Incorrect DQ_DEF column names:"
                                strg += " Expected %s;" % str(testnames)
                                strg += " Actual %s" % str(fieldnames)                
                                raise TypeError(strg)
                        else:
                            strg = "DQ_DEF binary table is the wrong size "
                            strg += "(expected 4 columns, got %d)." % tfields
                            raise TypeError(strg) 
                    else:
                        strg = "DQ_DEF HDU is not a binary table!"              
                        raise TypeError(strg)                    

        # NOTE: The following tidy up code now generates the astropy message:
        # WARNING:astropy:AstropyDeprecationWarning: Accessing an HDU after an HDUList is closed
        #hdulist.close()
        #del hdulist
    if failure_strg:
        raise TypeError("Not a valid JWST/MIRI data file: " + failure_strg)
        
# -----------------------------------------------------------------------------
#
# Data format conversion utilities. Many of these functions exist only to
# support old data formats.
#
def convert_detector(old_detector):
    """
    
    Convert an old detector string into a new string.
    
    NOTE: Only exists to support old data formats.
    
    """
    if old_detector:
        old_detector = old_detector.strip() # Strip off superfluous white space
    if old_detector == 'IM':
        return 'MIRIMAGE'
    elif old_detector == 'LW':
        return 'MIRIFULONG'
    elif old_detector == 'SW':
        return 'MIRIFUSHORT'
    else:
        # For any other value, return the same string
        return old_detector

def convert_band(old_band):
    """
    
    Convert an old band string into a new string.
    
    NOTE: Only exists to support old data formats.
     
    """
    if old_band:
        old_band = old_band.strip() # Strip off superfluous white space
    if old_band == 'A':
        return 'SHORT'
    elif old_band == 'B':
        return 'MEDIUM'
    elif old_band == 'C':
        return 'LONG'
    else:
        # For any other value, return the same string
        return old_band

def convert_cdp_2to3(infile, outfile, datatype=None, settings=None,
                     printmodel=False, version=None, 
                     update_filenames=False, overwrite=False):
    """
    
    Convert a MIRI calibration data product from CDP-2 format to CDP-3
    format.
    
    Throws an exception if the CDP cannot be converted.
    Issues a warnings message if some parts of the CDP are missing.
    
    NOTE: The function will fail if the output file already exists,
    unless called with overwrite=True.
    
    NOTE: Only exists to support old data formats.
         
    :Parameters:
    
    infile: str
        The name of the input file containing the olf format CDP.
    outfile: str
        The name of the output file to contain the converted CDP.
    datatype: str (option)
        The data type expected. If not given, it will be obtained from
        the TYPE or REFTYPE keyword in the file header.
    settings: str
        The detector settings used to create the CDP file.
        Allowed values are 'RAL1', 'JPL1' or 'ANY' or 'N/A'.
        If None, the detector settings are inferred from the MODELNAM
        and DETECTOR.
    printmodel: boolean (optional)
        Set to True to print the converted model before saving it.
        The default is False.
    version: string (optional)
        Set to a string to update the VERSION keyword in the output
        to this new input string.
    update_filenames: boolean (optional)
        Set to True to update the ORIGFILE and FILENAME keywords.
        The default is False.
    overwrite: boolean (optional)
        Set True to overwrite an existing file when making a
        temporary copy (not the original file itself).
        The default is False.
        
    :Returns:
    
    datatype: str
        The data type as converted. This may have been obtained from the
        file header.
    
    """
    assert isinstance(infile, str)
    assert isinstance(outfile, str)
    # Open the data model using the class derived from the data type.
    with open( init=infile, astype=datatype ) as datamodel:
         # Report the kind of model being processed
        if hasattr(datamodel.meta, 'reftype'):
            reftype = datamodel.meta.reftype
            strg = "The data model is of reference type \'%s\'" % \
                str(reftype)
        elif hasattr(datamodel.meta, 'model_type'):
            model_type = datamodel.meta.model_type
            strg = "The data model is of type \'%s\'" % \
                str(model_type)
        else:
            strg = "The data model class name is \'%s\'" % \
                datamodel.__class__.__name__
            
        detector = datamodel.meta.instrument.detector
        if detector:
            # Convert the detector to the new name
            new_detector = convert_detector( detector )
            datamodel.meta.instrument.detector = new_detector
            strg += " for detector \'%s\'-->\'%s\'" % \
                (str(detector), str(new_detector))
        mirifilter = datamodel.meta.instrument.filter
        if mirifilter:
            strg += " and filter \'%s\'" % str(mirifilter)
        
        # Convert to the new MRS band values
        if detector == 'SW' or detector == 'LW' or \
           detector == 'MIRIFUSHORT' or detector == 'MIRIFULONG':
            old_band = datamodel.meta.instrument.band
            new_band = convert_band( old_band )
            datamodel.meta.instrument.band = new_band
            strg += " and band \'%s\'-->\'%s\'" % \
                (str(old_band), str(new_band))
        print(strg)
        
        # If detector settings have not been given, infer them from
        # existing metadata.
        if not settings:
            model = datamodel.meta.instrument.model
            if model == 'JPL':
                settings = 'JPL1'
            elif model == 'VM':
                settings = 'RAL1'
            else:
                if detector == 'JPL':
                    settings = 'JPL1'
                else:
                    settings = 'N/A'
            print("Assumed detector settings: %s" % settings)
        else:
            print("Detector settings explicitly set to: %s" % settings)
            
        settings = settings.strip() # Strip off superflous whitespace
        datamodel.meta.instrument.detector_settings = settings
        if settings == 'ANY' or settings == 'N/A':
            # The detector settings do not matter.
            pass
        elif settings in DETECTOR_SETTINGS_DICT:
            # Set the FRMRSETS, ROWRSETS and RPCDELAY metadata.
            settings_tuple = DETECTOR_SETTINGS_DICT[settings]
            datamodel.meta.exposure.frame_resets = settings_tuple[0]
            datamodel.meta.exposure.reset_time = settings_tuple[1]
            datamodel.meta.exposure.rpc_delay = settings_tuple[2]
        else:
            strg = "\n***Unrecognised detector settings: %s" % settings
            warnings.warn(strg)
            
        # If necessary, ensure the data uses the new subarray keywords.
        try:
            subarray_name = datamodel.meta.subarray.name
            if subarray_name and subarray_name != 'FULL':
                print("Converting subarray keywords")
                datamodel.set_subarray_metadata(subarray_name)
        except:
            pass
        
        # If necessary, convert the data quality flags using a conversion table.
        if hasattr(datamodel,'_default_dq_def'):
            if hasattr(datamodel,'dq_def'):
                datamodel.dq_def = datamodel._default_dq_def
            elif hasattr(datamodel,'pixeldq_def'):
                datamodel.pixeldq_def = datamodel._default_dq_def
        datamodel.convert_flags( )
        
        # If a new 'version' has been provided, update VERSION
        if version:
            datamodel.meta.version = version
        
        if update_filenames:
            datamodel.meta.filename = os.path.basename(outfile)
            datamodel.meta.filename_original = os.path.basename(infile)
        
        if printmodel:
            print("\nConverted data model:")
            print(datamodel)
        
        # Save the converted data model
        print("Saving converted data object to %s." % outfile)
        datamodel.save( outfile, overwrite=overwrite )
        del datamodel
        return datatype

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the util module.")
    #import copy
    
    # Check that the data classes are looked up properly
    lookups = [['Bad'],
               ['Dark'],
               ['Flat'],
               ['Flat', 'ANY'],
               ['Flat', 'IM'],
               ['PixFlat'],
               ['Lin'],
               ['Flux'],
               ['Flux', 'ANY'],
               ['Flux', 'IM'],
               ['Flux', 'IM', 'F1800W'],
               ['Flux', 'IM', 'P750L'],
               ['Flux', 'IM', 'ANY'],
               ['Flux', 'SW'],
               ['Flux', 'LW'],
               ['Flux', 'LW', 'ANY'],
               ['AbsFlux'],
               ['AbsFlux', 'ANY'],
               ['AbsFlux', 'IM', 'ANY'],
               ['AbsFlux', 'IM', 'P750L'],
               ['SRF'],
               ['SRF', 'IM'],
               ['SRF', 'SW'],
               ['SRF', 'LW'],
               ['SRF', 'ANY'],
               ['PSF'],
               ['PSF', 'ANY'],
               ['PSF', 'IM'],
               ['PSF', 'IM', 'F1800W'],
               ['PSF', 'IM', 'P750L'],
               ['PSF', 'IM', 'ANY'],
               ['PSF', 'SW'],
               ['PSF', 'LW'],
               ['SILLYNAME']
               ]
    for lookup in lookups:
        strg = "Looking up " + str(lookup) + " - "
        gotclass = get_data_class( lookup )
        if gotclass is not None:
            strg += "got " + str(gotclass.__name__)
        else:
            strg += "got None (UNMATCHED) \t< ***"
        print(strg)
    
