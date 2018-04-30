#!/usr/bin/env python

"""

Module exposure_data - Contains the ExposureData class and
associated functions.

NOTE: This module now exists only for backwards compatibility.
It has been superceeded by the datamodels/miri_exposure_model.py

:History:
18 Jun 2010: Created.
28 Jun 2010: Reworked to recommended Python syntax standards.
05 Jul 2010: Changed name of test file.
21 Jul 2010: Added plotting methods. Corrected name space clash
             by using "intg" for integration instead of "int".
28 Jul 2010: Added plot_ramp and data size check. Corrected mistake
             in exposure data plot for large number of groups.
02 Aug 2010: Added ability to average groups and integrations.
06 Aug 2010: Documentation formatting problems corrected.
11 Aug 2010: FITS header updated.
01 Sep 2010: Improved attribute checking for public functions.
08 Sep 2010: write_file renamed to save to be compatible with other
             scasim classes. load_exposure_data function added.
             plotting functions extended.
24 Sep 2010: Check that plot clipping parameters are in the range
             0.0-1.0 and plot style is recognised.
07 Oct 2010: FITS header keywords updated after reading "Definition
             of the MIRI FM ILT Keywords" by T.W.Grundy.
             New statistics function added.
12 Oct 2010: FITSWriter option added and made the default.
             NFRAMES and GROUPGAP written properly to FITS header.
19 Oct 2010: Make sure exposure data is stored as floating point.
08 Nov 2010: Maximum data size increased to 160MB. Comments added.
09 Dec 2010: Maximum data size increased to 512MB.
26 Jan 2011: Added MIRI metadata to test.
26 Jan 2010: FITS header update moved from save function to new
             set_metadata function, so that header can be processed
             by a MiriMetadata object.
18 Feb 2011: Updated to match change of group name from DQ to DQL in
             mirikeyword.py.             
04 Mar 2011: Documentation tweaks to resolve bad formatting.
22 Mar 2011: Plotting functions modified to use matplotlib figures
             and axes more flexibly (unfortunately the colorbar
             function doesn't work properly).
01 Apr 2011: Documentation formatting problems corrected.
14 Jul 2011: Use of exceptions made more consistent and documented.
14 Sep 2011: FITS files written in 16-bit unsigned integer format
             (BITPIX=16) rather than floating point (BITPIX=-32).
             The latter causes problems for the DHAS miri_sloper
             V5.0.5 and later.
14 Sep 2011: Modified to keep up with the separation of
             miri.miritools.miridata into miri.miritools.metadata and
             miri.miritools.miricombdata
20 Oct 2011: Destructor added to remove large objects.
31 Oct 2011: Deprecated use of the .has_key() dictionary method removed.
15 Nov 2011: Worked around bug in pyfits V3.0 - all HDU creation arguments
             referenced by keyword rather than by name.
11 Jan 2012: Renamed symbols to avoid potential name clash:
             format-->file_format, max-->maxval, and min-->minval.
04 Mar 2012: add_to_header replaced with to_fits_header
22 Mar 2012: Metadata group import now specifies a wildcard.
02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
26 Apr 2012: Brought up to date with changes in Metadata class.
13 Nov 2012: Major restructuring of the package folder. Import statements
             updated.
14 Jan 2013: olddataproduct subpackage split from dataproduct.
05 Jun 2013: Removed verbose flag and rationalised attribute names.
13 Jun 2013: Changed metadata to match the new data model. New functions
             for accessing the FITS header added to be compatible with new
             data model.
18 Jun 2013: Create ExposureData objects with non-null header objects.
09 Jun 2014: Fixed pyfits deprecation warnings. Added ability to generate
             slope data and save in level 2 format.
19 Jun 2014: Use uint32 format to prevent data being truncated unnecessarily.
14 Jul 2014: Legacy MiriMetadata class moved to scasim, to remove
             dependence on the olddatamodel package.
06 Mar 2015: Fixed bug in FITSWriter format output.
27 May 2015: Replaced pyfits with astropy.io.fits
11 Jun 2015: Added get_history and set_history methods for forwards
             compatibility.
06 Aug 2015: Added set_exposure_times method for forwards compatibility.
02 Sep 2015: "slope" option added to output data shape.
             Removed dependence on old legacy MiriMetadata class.
08 Sep 2015: Metadata class and get_file_header function moved here.
             Made compatible with Python 3.
22 Sep 2015: get_exposure_times method added.
07 Oct 2015: Made exception catching Python 3 compatible.
03 Nov 2015: Changed to use floating point exposure times.
09 Mar 2016: Make the message about exceeding the data size limit a
             warning rather than an exception.
28 Sep 2016: miri.miritools.dataproduct renamed miri.datamodels and
             miri.miritools renamed miri.tools.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
27 Jul 2017: Added the ability to add DARK data in situ (for simulation).
08 Aug 2017: Log an error instead of raising an exception when there are
             missing header keywords.
11 Sep 2017: Ensure EXPSTART, EXPMID and EXPEND keywords are MJD.
30 Apr 2018: Replaced xrange with range for Python 3.

@author: Steven Beard

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.exposure_data") # Get a default parent logger

import numpy as np
import scipy.stats
from astropy.extern import six
import astropy.io.fits as pyfits

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt
    
# Maximum acceptable data size in the output file (in MB). A warning
# is issued if the user attempts to create a larger file. This
# prevents the simulator from getting choked creating huge files
# when inappropriate exposure time and readout mode combinations are
# used by accident. (NOTE: Ground based data is limited to around
# <300 reads, or ~1GB. In orbit data will be limited to around <30
# reads, or ~100MB.)
_MAX_DATA_SIZE_MB = 512.0

# The maximum number of plots that can be sensibly fitted
# into one matplotlib figure. Depends on monitor resolution
# and eyesight.
_MAX_PLOTS = 20

from miri.simulators.integrators import linear_regression


# Simple Metadata class used for legacy data models
class Metadata(object):
    """
    
    Class Metadata - a simple dictionary-like object which maintains a
    collection of keyword/value pairs plus a list of comment and history
    records.

    :Parameters:
        
    description: string, optional
        An optional description string for the metadata.
            
    """
    def __init__(self, description=''):
        """
        
        Create a new empty Metadata object.
        
        For parameters see class description.
       
        """
        # A string describing this set of metadata
        self.description = description
        
        # These attributes keep a record of the keywords and values
        # and keyword comments currently contained in the metadata.
        self._kwdict = {}
        self._kwcomment = {}
        self._keywords = []
        
        # Lists of comment and history records.    
        self.comments = []
        self.history = []

    def set_key_comment(self, keyword, comment):
        """
        
        Set the comment associated with a particular keyword.
        The keyword must already exist and have a value.
        
        :Parameters:
        
        keywords: string
            The keyword
        comment: string
            Comment to be associated with the keyword.
        
        """
        if keyword in self:
            self._kwcomment[keyword] = str(comment)
        else:
            strg = "Keyword %s not present in metadata" % keyword
            if self.description:
                strg += " \'%s\'" % self.description
            raise KeyError(strg)

    def __setitem__(self, keyword, value):
        """
        
        Set a metadata keyword item with the operation::
        
            metadata[keyword] = properties
            
        Properties can be a tuple containing (value,comment) for
        the keyword or a single item containing the value.
        
        """
        self._kwdict[keyword] = value
        self._keywords.append(keyword)
            
    def __contains__(self, keyword):
        """
        
        Returns True if the given keyword is contained in the metadata.
        Called in response to the operation::
        
            keyword in metadata
        
        """
        return (keyword in self._kwdict)
    
    def __delitem__(self, keyword):
        """
        
        Remove the given keyword entry from the metadata with the operation::
        
            del metadata[keyword]
        
        """
        if keyword in self:
            del self._kwdict[keyword]
            del self._kwcomment[keyword]
            self._keywords.remove(keyword)

    def __getitem__(self, keyword):
        """
        
        Return the value associated with the given keyword with
        the operation::
        
            value = metadata[keyword]
            
        This is a wrapper for the get_value() method.
        
        """
        if keyword in self:
            return self._kwdict[keyword]
        else:
            strg = "Keyword %s not present in metadata" % keyword
            if self.description:
                strg += " \'%s\'" % self.description
            raise KeyError(strg)

    def keys(self):
        """
        
        Return the list of all active keywords contained in the metadata.
        Unlike a regular Python dictionary, the list is returned in
        the same order in which the metadata was defined.

        :Returns:
        
        kwlist: list of keywords
            A list of all the active keywords present in the metadata, in
            the same order in which they were added.
         
        """
        return self._keywords

    def add_comment(self, comment):
        """
        
        Add a comment to the list of comments stored in the metadata.
        The comments are added to the FITS header when the
        to_fits_header method is called.
                 
        :Parameters:
    
        comment: string
            A comment string to be added to the list (ignored if
            a null string is given).
        
        """
        if comment:
            self.comments.append(comment)

    def get_comments(self):
        """
        
        Return a list of strings containing all the comments added.
        
        """
        return self.comments

    def has_comment(self, comment):
        """
        
        Check whether a give comment is (word for word) already
        present in the metadata.
        
        :Parameters:
        
        comment: string
            Comment string.
            
        :Returns:
        
        present: bool
            True if the comment string is present.
            A null comment string will always return True
        
        """
        if comment:
            return comment in self.comments
        else:
            return True

    def add_history(self, history):
        """
        
        Add a history record to the list stored in the metadata.
                 
        :Parameters:
    
        history, str
            A history description to be added to the list (ignored if
            a null string is given).
        
        """
        if history:
            self.history.append(history)

    def get_history(self):
        """
        
        Return a list of strings containing all the history records added.
        
        """
        return self.history

    def has_history(self, history):
        """
        
        Check whether a give history string is (word for word) already
        present in the metadata.
        
        :Parameters:
        
        history: string
            History string.
            
        :Returns:
        
        present: bool
            True if the history string is present.
            A null history string will always return True
        
        """
        if history:
            return history in self.history
        else:
            return True

    def _comment_valid(self, comment):
        """
        
        Helper function which determines whether a string is a valid
        comment.
        
        """
        # A comment must be a string.
        if isinstance(comment,six.string_types):
            return True
        else:
            return False

    def _extract_comment(self, keywords, comment):
        """
        
        Helper function which extracts the relevant parts of a
        comment string, missing out the initial keyword.
        
        """
        # Convert the comment into a string and split it into words.
        comment = str(comment)
        comwords = comment.split()
        # Check whether a single keyword or a list of keywords has
        # been provided.
        if isinstance(keywords, (tuple,list)):
            for keyword in keywords:
            # If the first word matches the given keyword, remove it.
                if (len(comwords) > 0) and (comwords[0] == keyword):
                    comment = " ".join(comwords[1:])
                    comwords = comment.split()
        else:
            # Assume a single keyword
            # If the first word matches the given keyword, remove it.
            if (len(comwords) > 0) and (comwords[0] == keywords):
                comment = " ".join(comwords[1:])

        # Return a copy of the comment with leading and trailing space removed.
        rscomment = comment.strip()
        return rscomment

    def _key_name_valid(self, key, maxlength=1024):
        """
        
        Helper function which determines whether a string is a valid
        keyword name.
        
        """
        # A keyword must be a string.
        if isinstance(key,six.string_types):
            # It cannot be white space or a null string or be longer
            # than the maximum length.
            if key.isspace() or (len(key) == 0) or (len(key) > maxlength):
                return False
            else:
                return True            
        else:
            return False

    FITS_MAX_KEYWORD_LENGTH = 8
    FITS_MAX_COMMENT_LENGTH = 72

    def _key_name_valid_fits(self, key):
        """
        
        Helper function which determines whether a string is a valid
        keyword name for use in a FITS file.
        
        """
        # It cannot be white space and must be shorter than
        # the maximum length limit.
        valid = self._key_name_valid(key, maxlength=self.FITS_MAX_KEYWORD_LENGTH)
        if not valid:
            return False
        else:
            # The first character must be alphabetic.
            if key[0].isalpha():
                return True
            else:
                return False

    def from_data_object(self, data_object, hdu_name='PRIMARY'):
        """
        
        Create the metadata from the metadata found in a new MIRI
        data object.
        
        :Parameters:
        
        data_object: MIRI data object
            The data object from which to extract the metadata.
        hdu_name: string, optional
            The name of the HDU the metadata should be associated with.
            
        :Returns:
        
        keys_added: int
            The number of keywords added to the metadata (not including
            COMMENT and HISTORY keywords).

        """
        keys_added = 0
        if not hasattr(data_object, 'fits_metadata_dict'):
            strg = "Attribute \'%s\' is not a MIRI data object" % \
                data_object.__class__.__name__
            raise TypeError(strg)

        metadict = data_object.fits_metadata_dict()
        metakeys = metadict.keys()
        for keyw in metakeys:
            (fitshdu, fitskey, fitscomment) = metadict[keyw]
            value = data_object[keyw]
            if value is not None and fitshdu == hdu_name:
                keys_added += 1
                self[fitskey] = value
                if fitscomment:
                    self.set_key_comment(fitskey, fitscomment)
                        
        comments = data_object.get_fits_comments(hdu_name=hdu_name)
        for comment in comments:
            self.add_comment( comment )

        histories = data_object.get_history()
        for history in histories:
            self.add_history( history )

        return keys_added

    def from_fits_header(self, fitsheader):
        """
        
        Create the metadata from the given pyFits FITS header.
        
        :Parameters:
        
        fitsheader: pyfits Header object
            The FITS header from which the metadata are to be read.
            
        :Returns:
        
        keys_added: int
            The number of keywords added to the metadata (not including
            COMMENT and HISTORY keywords).

        """
        keys_added = 0
        # Populate metadata from the given FITS header. NOTE: There
        # can be multiple COMMENT and HISTORY records, so these need
        # to be copied differently. Skip blank records.
        for kw in fitsheader:
            if kw == 'COMMENT' or kw == 'HISTORY' or \
               (len(kw) == 0) or kw.isspace():
                pass
            else:
                keys_added += 1
                self[kw] = fitsheader[kw]
                # TODO: Extract FITS comment field?
        
        # Copy all COMMENT and HISTORY records with trailing white
        # space and initial labels (and superflous '=' signs) removed.
        # Empty records are ignored.
        if 'COMMENT' in fitsheader:
#             comments = fitsheader.get_comment()
            comments = fitsheader['COMMENT']
            for comment in comments:
                rscomment = self._extract_comment(('COMMENT','='), comment)
                if rscomment:
                    self.add_comment(rscomment)

        if 'HISTORY' in fitsheader:
#             histories = fitsheader.get_history()
            histories = fitsheader['HISTORY']
            for history in histories:
                rshistory = self._extract_comment(('HISTORY','='), history)
                if rshistory:
                    self.add_history(rshistory)
        
        return keys_added
 
    def to_fits_header(self, fitsheader, overwrite=False, add_description=True):
        """
        
        Add metadata to the given pyfits FITS header object
        in the correct keyword order. Only keywords which have
        defined values are written, plus any COMMENT and HISTORY
        descriptions stored within the metadata.
        
        :Parameters:
        
        fitsheader: pyfits Header object
            The FITS header to which metadata should be
            added. If this is None, a new header is created.
        overwrite: bool, optional, default=False
            If True, if any keywords within the metadata match
            existing keyword in the FITS header the old keywords
            will be overwritten.
            If False, existing keywords within the FITS header
            will be protected.
        add_description: bool, optional, default=True
            If True, the description associated with this metadata
            will be added to the FITS header as a comment keyword.
        
        :Returns:
        
        fitsheader: pyfits Header object
            The FITS header with added metadata.        
        
        """
        # Create a new header if an existing one has not been provided.
        if fitsheader is None:
            fitsheader = pyfits.Header()
        
        if add_description:
            fitsheader.add_comment(self.description)
        
        for key in self.keys():

            # Check the keyword name is a valid FITS keyword.
            if not self._key_name_valid_fits(key):
                strg = "String \'%s\' is not a valid FITS keyword." % \
                    str(key)
                raise ValueError(strg)
           
            # Skip an existing keyword if overwrite is False.
            if not overwrite:
                if key in fitsheader:
                    continue
            
            value = self._kwdict[key]
            fitsheader[key] = value
            
        # Add comments and history records to the FITS header.
        # Split multi-line records into individual lines to avoid
        # an "unprintable string" problem with pyfits. Also detect
        # multi-line records that have a CONTINUE keyword but no
        # newline.
        if self.comments:
            for comment in self.comments:
                comlines = comment.splitlines()
                for comline in comlines:
                    comparts = comline.split('CONTINUE')
                    for compart in comparts:
                        if compart:
                            # Remove superflous labels and '=' statements.
                            newline = \
                                self._extract_comment(('CONTINUE','='), compart)
                            fitsheader.add_comment(newline)
                
        if self.history:
            for history in self.history:
                histlines = history.splitlines()
                for histline in histlines:
                    histparts = histline.split('CONTINUE')
                    for histpart in histparts:
                        if histpart:
                            # Remove superflous labels and '=' statements.
                            newline = \
                                self._extract_comment(('CONTINUE','='), histpart)
                            fitsheader.add_history(newline)
               
        return fitsheader

    def __str__(self):
        """
        
        Return a string description of a set of metadata. This function
        is called when an object is printed or when a statement
        contains str(object).
            
        :Returns:
        
        strg: string
            A readable string describing the object.
        
        """
        strg = 'Metadata'
        if self.description:
            strg += " \'%s\'" % self.description
        strg += ":\n"
        # Add a description for each known keyword
        for key in self.keys():
            value = self._kwdict[key]
            strg += "%8s = %s\n" % (key, str(value))
        # Add any comment and history records contained in the
        # metadata.
        if self.comments:
            for comment in self.comments:
                strg += "COMMENT: %s\n" % comment
        if self.history:
            for history in self.history:
                strg += "HISTORY: %s\n" % history
        return strg

    def __repr__(self):
        """
        
        Return a representation of the object. This function is called
        when the name of an object is entered interactively, or when
        a statement contains repr(object).
            
        :Returns:
        
        strg: string
            An unambiguous string describing the type and content of the object.
        
        """
        strg = type(self).__name__ + '('
        strg += 'kwdict=%s' % repr(self._kwdict)
        strg += ')'
        return strg


class ExposureData(object):
    """
    
    Class ExposureData - Manages the simulated exposure data and the
    level 1 FITS file output.
    
    LEGACY DATA FORMAT WHICH IS COMPATIBLE WITH DHAS.
                 
    :Parameters:
    
    rows: int
        The number of rows making up the data.
    columns: int
        The number of columns making up the data.
    ngroups: int
        The number of readout groups making up each integration.
    nints: int
        The number of integrations making up each exposure.
    readpatt: string
        Name of detector readout pattern.
    grpavg: int, optional, default=1
        The number of groups to be averaged to reduce the file size.
        *NOTE: An error will be reported if ngroups does not divide
        exactly by grpavg.*
    intavg: int, optional, default=1
        The number of integrations to be averaged to reduce the file
        size. 
        *NOTE: An error will be reported if nints does not divide
        exactly by intavg.*
    nframes: int, optional, default=1
        The number of frames per group. Normally 1 for MIRI data.
    groupgap: int, optional, default=0
        The number of dropped frames in between groups.
        
    :Requires:
    
    pyfits
        
    :Raises:

    A warning is issued if the combination of initialisation parameters
    will generate too much data.
    
    ValueError
        Raised if any of the initialisation parameters are out of range.

    """
    def __init__(self, rows, columns, ngroups, nints, readpatt, grpavg=1,
                 intavg=1, nframes=1, groupgap=0):
        """
        
        Constructor for class ExposureData.
        
        Parameters: See class doc string.
        
        """
        if int(rows) <= 0 or int(columns) <= 0:
            strg = "The exposure data must have a non-zero size."
            raise ValueError(strg)
        self.rows = int(rows)
        self.columns = int(columns)
        
        self.ngroups = int(ngroups)
        if self.ngroups <= 0:
            strg = "There must be at least one group."
            raise ValueError(strg)
        self.grpavg = int(grpavg)
        self.ngroups_file = self.ngroups // self.grpavg
        
        self.nints = int(nints)
        if self.nints <= 0:
            strg = "There must be at least one integration."
            raise ValueError(strg)
        self.intavg = int(intavg)
        self.nints_file = self.nints // self.intavg
        
        self.nframes = int(nframes)
        self.groupgap = int(groupgap)
        self._written = False
        self._primary_header = pyfits.Header()

        self._data_header = pyfits.Header()
        self.readpatt = readpatt
        
        # Check the integrity of the requested data averaging. It is
        # assumed the caller will provide nicely rounded values.
        if self.ngroups % self.grpavg != 0:
            strg = "%d groups does not divide exactly by %d. " % \
                (self.ngroups, self.grpavg)
            strg += "Please adjust the number of groups."
            raise ValueError(strg)    
        if self.nints % self.intavg != 0:
            strg = "%d integrations does not divide exactly by %d. " % \
                (self.nints, self.intavg)
            strg += "Please adjust the number of integrations."
            raise ValueError(strg)    
        
        # Check the file data size is manageable.
        data_size_mb = 4 * self.columns * self.rows * \
            self.ngroups_file * self.nints_file / 1.024E6
        if data_size_mb > _MAX_DATA_SIZE_MB:
            strg = "This simulation will generate %.1fMB of data "  % \
                data_size_mb
            strg += "(> %.1fMB)." % _MAX_DATA_SIZE_MB
            strg += "\nConsider rerunning with a reduced exposure time or " \
                "adjust the readout mode (e.g. FAST --> SLOW or FASTGRPAVG)."
            LOGGER.warning(strg)

        # Create the 4-D data to hold all the data readings.
        # NOTE: The ordering of the dimensions of the array is transposed
        # when reading or writing a FITS file, so the FITS data will have
        # dimensions: NAXIS1=columns, NAXIS2=rows, NAXIS3=ngroups and
        # NAXIS4=nints.
        # NOTE: The new MIRI exposure model uses floating point data.
        # MEMORY MANAGEMENT: Don't create the array until needed.
        self.datashape = [self.nints, self.ngroups, self.rows, self.columns]
        self.datasize = self.nints * self.ngroups * self.rows * self.columns
        self.data = None
#         self.data = np.zeros(self.datashape, dtype=np.uint32)
        
        # Create a place holder to store averaged SCI data, if needed.
        self._data_averaged = None
        
    def __del__(self):
        """
        
        Destructor for class ExposureData.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        try:
            # Objects created in the constructor
            del self.data
            
            # Objects created in class methods
            if self._data_averaged is not None:
                del self._data_averaged
        except Exception:
            pass

    def set_group(self, data, group, intg):
        """
        
        Add a new readout to the SCI data array.
        
        :Parameters:
        
        data: array_like (converted to uint32)
            An array containing the readout data to be added to exposure
            SCI data. It must have the same shape as one frame of the
            exposure SCI data (i.e. a 2-D array with the same number of
            rows and columns specified in the ExposureData constructor).
        group: int
             The group to which the data readout belongs.
        intg: int
            The integration to which the data readout belongs.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        # Copy the input data to a 2-D slice for this group/intg combination.
        # NOTE: This only works if data array is broadcastable so the shape
        # of the data array is checked.
        data = np.asarray(data, dtype=np.uint32)
        expected_shape = (self.rows, self.columns)
        if data.shape == expected_shape:
            # MEMORY MANAGEMENT: Don't create the array until needed.
            if self.data is None:
                self.data = np.zeros(self.datashape, dtype=np.uint32)    
            self.data[intg,group,:,:] = data
        else:
            strg = "Data array has the wrong shape: " \
                "(%d, %d) instead of (%d, %d)." % (data.shape[0],
                                                   data.shape[1],
                                                   expected_shape[0],
                                                   expected_shape[1])
            raise TypeError(strg)
        
        # There is more data waiting to be written.
        self._written = False

    def set_exposure(self, data):
        """
        
        Define the entire SCI data array in one go.
        
        :Parameters:
        
        data: array_like (converted to uint32)
            An array containing the SCI data. It's size must match the
            number of rows, columns, groups and integrations defined
            when the ExposureData object was created.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        # The given data must be the same size as the existing SCI_data array.
        # MEMORY MANAGEMENT: Don't create the array until needed.
        data = np.asarray(data, dtype=np.uint32)
#         if data.size == self.data.size:
        if data.size == self.datasize:
            if self.data is not None:
                del self.data
            self.data = data
            # Reshape the array, in case it has been stored as a cube.
            self.data.shape = (self.nints, self.ngroups, self.rows,
                               self.columns)
        else:
            strg = "Data array has the wrong size: " \
                "%d instead of %d." % (data.size, self.data.size)
            raise TypeError(strg)
        
        # There is more data waiting to be written.
        self._written = False

    def set_exposure_times(self, exposure_time=None, duration=None,
                           start_time='NOW', mid_time=None, end_time=None):
        """
        
        Convenience function to define exposure times in the metadata.
        If called with no arguments, the function will attempt to set
        the exposure end time from the existing start time and duration.
        
        Exposure times are stored as MJD in days.
        
        :Parameters:
        
        exposure_time: number, optional
            Exposure time.
        duration: number, optional
            The exposure duration. Defaults to the same as the exposure time.
        start_time: str or float, optional
            Date/time of start of exposure (MJD days). Strings other than
            'NOW' are converted to floating point. If set to 'NOW', the
            current date-time is used. By default this is not set.
        mid_time: str or float, optional
            Date/time of start of exposure (MJD days). Strings other than
            'NOW' are converted to floating point. If set to 'NOW', the
            current date-time is used. By default this is not set.
        end_time: str or float, optional
            Date/time of start of exposure (MJD days). Strings other than
            'NOW' are converted to floating point. If set to 'NOW', the
            current date-time is used. If not defined, the end time is set
            to start_time + duration, or failing that is not set at all.
        
        """
        import time, datetime
        # Modified Julian date of the "zero epoch" of the time library (1/1/70)
        MJD_ZEROPOINT = 40587.0
        # Number of seconds per day.
        SECONDS_PER_DAY = 86400.0
        if exposure_time is not None:
            self._primary_header['EXPTIME'] = exposure_time
            self._primary_header['EFFEXPTM'] = exposure_time
        if duration is not None:
            self._primary_header['DURATION'] = duration
        elif exposure_time is not None:
            self._primary_header['DURATION'] = exposure_time
                
        if start_time == 'NOW':
#             start_time = time.strftime('%Y-%m-%dT%H:%M:%S')
            start_time = MJD_ZEROPOINT + (time.time()/SECONDS_PER_DAY)
        if start_time is not None:
            self._primary_header['EXPSTART'] = start_time
                
        if mid_time == 'NOW':
#             mid_time = time.strftime('%Y-%m-%dT%H:%M:%S')
            mid_time = MJD_ZEROPOINT + (time.time()/SECONDS_PER_DAY)
        if mid_time is not None:
            self._primary_header['EXPMID'] = mid_time
                
        if end_time == 'NOW':
#             end_time = time.strftime('%Y-%m-%dT%H:%M:%S')
            end_time = MJD_ZEROPOINT + (time.time()/SECONDS_PER_DAY)
        elif self._primary_header['EXPSTART'] is not None and \
             self._primary_header['DURATION'] is not None and end_time is None:
            # Set the end time to start_time + duration
            # FIXME: Does not work. Low priority, since data model is deprecated.
#             end_time = datetime.datetime(self._primary_header['EXPSTART']) + \
#                 datetime.timedelta(seconds=self._primary_header['DURATION'])
                # Set the end time to start_time + duration
            end_time = self._primary_header['EXPSTART'] + \
                            self._primary_header['DURATION']
        if end_time is not None:
            self._primary_header['EXPEND'] = end_time

    def get_exposure_times(self):
        """
        
        Return the exposure time metadata to the caller
        
        :Returns:
        
        (exposure_time, duration, start_time, mid_time, end_time): tuple of 5 floats
            The exposure time in seconds, clock duration in seconds, exposure
            start time, mid time and end time in MJD days. Undefined values are
            returned None.
        
        """
        if 'EFFEXPTM' in self._primary_header:
            exposure_time = self._primary_header['EFFEXPTM']
        else:
            exposure_time = None
        if 'DURATION' in self._primary_header:
            duration = self._primary_header['DURATION']
        else:
            duration = None
        if 'EXPSTART' in self._primary_header:
            start_time = self._primary_header['EXPSTART']
        else:
            start_time = None
        if 'EXPMID' in self._primary_header:
            mid_time = self._primary_header['EXPMID']
        else:
            mid_time = None
        if 'EXPEND' in self._primary_header:
            end_time = self._primary_header['EXPEND']
        else:
            end_time = None
        return (exposure_time, duration, start_time, mid_time, end_time)
   
    def set_header(self, header):
        """
        
        Set the primary FITS header.
        Any existing header is deleted.
        
        :Parameters:
        
        header: pyfits header object
            The primary FITS header.
            
        """
        if self._primary_header is not None:
            del self._primary_header
        self._primary_header = header
        
    def set_data_header(self, header):
        """
        
        Set the header for the SCI data extensions.
        Any existing header is deleted.
        
        :Parameters:
        
        header: pyfits header object
            The SCI data extension FITS header.
            
        """
        if self._data_header is not None:
            del self._data_header
        self._data_header = header

    def add_fits_comment(self, comment, hdu_name='PRIMARY'):
        """
        
        Add a comment to the comment records contained in the metadata
        associated with a given FITS HDU.
        
        :Parameters:
        
        comment: string
            A string containing a comment.
        hdu_name: string, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
            If None, any HDU is matched, but the keyword must be
            unique within the data structure.
                        
        """
        # Convert unusual entries into a plain string.
        if not isinstance(comment, str):
            comment = str(comment)
        if hdu_name == 'DATA':
            self._data_header.add_comment(comment)
        else:
            self._primary_header.add_comment(comment)

    def add_comment(self, comment):
        """
        
        General purpose add_comment. Simply calls add_fits_comment.
        
        """
        self.add_fits_comment(comment)

    def get_fits_comments(self):
        """
        
        Locates all the COMMENT metadata and returns them as a list.
        Comments found in the PRIMARY metadata are listed first.
        
        :Parameters:
        
        None
        
        :Returns:
        
        comments: list of list of strings
            All the comments found in the metadata. Each item in the top
            level list corresponds to one HDU and each item in the second
            level list corresponds to a separate COMMMENT entry.
                
        """
        comments = self._primary_header.comments()
        comments += self._data_header.comments()
        return comments
    
    def get_fits_comments_str(self):
        """
        
        Return a structured string containing the comments
        associated with the data product.
        
        """
        strg = ''
        commentlist = self.get_fits_comments()
        for comment in commentlist:
            strg += str(comment[0]) + "\n"
            for ii in range(1,len(comment)):
                strg += "\t" + str(comment[ii]) + "\n"
        return strg

    # Forwards compatibility
    def add_history(self, history):
        return self.add_fits_history(history)

    def add_fits_history(self, history, hdu_name='PRIMARY'):
        """
        
        Add a history record to the metadata
        associated with a given FITS HDU.
        
        :Parameters:
        
        history: string
            A string containing a history record.
        hdu_name: string, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
            If None, any HDU is matched, but the keyword must be
            unique within the data structure.
                        
        """
        # Convert unusual entries into a plain string.
        if not isinstance(history, str):
            history = str(history)
        if hdu_name == 'DATA':
            self._data_header.add_history(history)
        else:
            self._primary_header.add_history(history)

    # Forwards compatibility
    def get_history(self):
        return self.get_fits_history()

    def get_fits_history(self):
        """
        
        Locates all the HISTORY metadata and returns them as a list.
        History records found in the PRIMARY metadata are listed first.
        
        :Parameters:
        
        None
        
        :Returns:
        
        history: list of list of strings
            All the history found in the metadata. Each item in the top
            level list corresponds to one HDU and each item in the second
            level list corresponds to a separate HISTORY entry.
                
        """
        history = self._primary_header.history()
        history += self._data_header.history()
        return history

    # Forwards compatibility
    def get_history_str(self):
        return self.get_fits_history_str()

    def get_fits_history_str(self):
        """
        
        Return a structured string containing the history records
        associated with the data product.
        
        """
        strg = ''
        historylist = self.get_fits_history()
        for history in historylist:
            strg += str(history[0]) + "\n"
            for ii in range(1,len(history)):
                strg += "\t" + str(history[ii]) + "\n"
        return strg

    def get_fits_keyword(self, keyword, hdu_name='PRIMARY'):
        """
  
        Return the value within the metadata associated with
        a FITS keyword. Only one value is returned.
     
        :Parameters:
        
        keyword: string
            The keyword to be obtained.
        hdu_name: string, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
        
        :Returns:
        
        value: object
            An object containing the current value associated with the keyword.
            If multiple occurences are matched, only the first will be returned.
        
        """
        if hdu_name == 'DATA':
            if keyword in self._data_header:
                return self._data_header[keyword]
            else:
                strg = "Keyword %s not found in %s header" % (keyword, hdu_name)
                raise KeyError(strg)
        else:
            if keyword in self._primary_header:
                return self._primary_header[keyword]
            else:
                strg = "Keyword %s not found in %s header" % (keyword, hdu_name)
                raise KeyError(strg)

    def set_fits_keyword(self, keyword, value, hdu_name='PRIMARY'):
        """
  
        Set the metadata associated with a FITS keyword to a given value.
     
        :Parameters:
        
        keyword: string
            The keyword to be updated. The keyword must be unique.
        value: object
            The value to set the metadata to. 
        hdu_name: string, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
        
        """
        # Convert unusual values into a string
        if six.PY2:
            usual_types = (str,float,int)
        else:
            usual_types = (str,float,int)
        if not isinstance(value, usual_types):
            value = str(value)
        if hdu_name == 'DATA':
            self._data_header[keyword] = value
        else:
            self._primary_header[keyword] = value
    
    def add_dark(self, darkarray ):
        """
        
        Add a DARK array to exposure data in-situ (for simulation purposes).
        The DARK array is assumed to have been read from a MIRI DARK CDP file.
        The exposure data and DARK are assumed both to be in DN units.
        
        NOTE: The DARK CDP often contains negative values, especially
        in bad pixel zones. The result is clipped to ensure these
        negative values don't generate a zero or negative result.
        
        :Parameters:
        
        darkarray: numpy array
            An array containing the DARK data. Its size must match the
            number of rows, columns, groups and integrations defined
            when the ExposureDataProduct object was created.
            A 4-D array of shape (2, ngroups, columns, rows) is
            expected, where the first dimension provides DARK data for
            the first and subsequent integrations.
            A 3-D array of shape (groups, columns, rows) is also acceptable
            but will be applied to all integrations.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        darkarray = np.asarray(darkarray)
        LOGGER.debug("Adding DARK array of shape %s" % str(darkarray.shape) + \
              " to data array of shape %s" % str(self.data.shape))

        ngroups = self.data.shape[1]
        if darkarray.ndim == 4:
            # 4-D data has been provided.
            # The DARK must have at least as many groups as the exposure.
            if darkarray.shape[1] < ngroups:
                raise ValueError("DARK data has insufficient groups.")

            # The dark is expected to be smaller than the data because
            # of the reference rows.
            darkrows = darkarray.shape[-2]
            darkcols = darkarray.shape[-1]

            # Add the DARK to the first integration
            # Skip the reference rows (assuming they are at the top)
            self.data[0, :, :darkrows, :darkcols] = \
                    self.data[0, :, :darkrows, :darkcols] + \
                    darkarray[0, :ngroups, :, :]
            # Add the DARK to the second and subsequent integrations.
            self.data[1:, :, :darkrows, :darkcols] = \
                    self.data[1:, :, :darkrows, :darkcols] + \
                    darkarray[1, :ngroups, :, :]

        elif darkarray.ndim == 3:
            # 3-D data has been provided.
            # The DARK must have at least as many groups as the exposure.
            if darkarray.shape[0] < ngroups:
                raise ValueError("DARK data has insufficient groups.")

            # The dark is expected to be smaller than the data because
            # of the reference rows.
            darkrows = darkarray.shape[-2]
            darkcols = darkarray.shape[-1]

            # Add the DARK to the all groups.
            # Skip the reference rows (assuming they are at the top)
            self.data[:, :, :darkrows, :darkcols] = \
                    self.data[:, :, :darkrows, :darkcols] + \
                    darkarray[:ngroups, :, :]
 
        else:
            strg = "DARK data array has the wrong shape "
            strg += "(%s compared with %s)." % (str(darkarray.shape),
                                                str(self.data.shape))
            raise TypeError(strg)
        
        # Adding the dark must not allow the resulting data to go negative
        # or contain zeros. It must also not allow the data to go above
        # the maximum value for 16-bit telemetry data.
        whereneg = np.where( self.data < 1.0 )
        if whereneg and (len(whereneg[0]) > 0):
            strg = "%d negative pixels after adding DARK." % len(whereneg[0])
            LOGGER.debug(strg)
            self.data[whereneg] = 1.0
        whereclipped = np.where( self.data > 65535.0 )
        if whereclipped and (len(whereclipped[0]) > 0):
            strg = "%d saturated pixels after adding DARK." % len(whereclipped[0])
            LOGGER.debug(strg)
            self.data[whereclipped] = 65535.0 

    def save(self, filename, fileformat='level1', datashape='cube',
             overwrite=False, recalculate=False):
        """
        
        Creates and writes a level 1 FITS file with the exposure data.
        
        :Parameters:
        
        filename: string
            The name of the file to be created.
        fileformat: string, optional, default='FITSWriter'
            The kind of file format to be written.
            
            * 'FITSWriter' - emulate the format written by the
              FITSWriter during MIRI VM and FM tests and read by the
              DHAS.
            * 'level1' - emulate the level 1 FITS format used by
              the JWST DMS data pipeline.
            * 'level2' - write simple slope data.
              
        datashape: string, optional, default='cube'
            The SCI data format to be written.
            
            * 'hypercube' - write the SCI data to a 4 dimensional FITS
              image with separate columns x rows x groups x
              integrations dimensions.
            * 'cube' - append the groups and integrations to make a 3
              dimensional FITS image with columns x rows x (groups and
              integrations) dimensions.
            * 'slope' - write fitted slope data.
              
        overwrite: bool, optional, default=False
            Parameter passed to pyfits.HDUlist.writeto
        recalculate: boolean, optional, default=False
            If True the averaged data are recalculated even if already
            available.
            
        :Raises:
    
        IOError
            Raised if the output file could not be written.
            
        """
        if self.data is None:
            raise TypeError("No exposure data to be saved.")
        
        # Create a primary FITS header if one doesn't already exist
        # and append the file name, copying any previous file name to
        # ORIGFILE.
        if self._primary_header is None:
            header = pyfits.Header()
        else:
            header = self._primary_header

        if self._data_header is None:
            self._data_header = pyfits.Header()
        
        try:
            origfile = header['FILENAME']
            header["ORIGFILE"] = (origfile, "Original file")
        except (KeyError, NameError):
            pass
        header["FILENAME"] = (filename, "File name")
        
        if 'READPATT' not in header:
            header["READPATT"] = (self.readpatt, "Detector readout pattern")
        
        if fileformat == 'level1' or fileformat == 'level2':
            primary_hdu = pyfits.PrimaryHDU(data=None, header=header)

            # Add some statistics to the SCI data HDU header.
            self._data_header["BUNIT"] = ("DN", "Data units")
            self._data_header["MIN"] = (self.data.min(), "Minimum value")
            self._data_header["MAX"] = (self.data.max(), "Maximum value")

        # Decide whether the SCI_data array contents need to be
        # averaged before writing to the file.
        if self.grpavg <= 1 and self.intavg <= 1:
            
            # No averaging
            if datashape == "cube":
                # Take a copy of the SCI data and reshape it into a cube
                datacube = self.data.copy()
                datacube.shape = (self.nints*self.ngroups, self.rows,
                                  self.columns)
                # TODO: Add World Coordinates?
                if fileformat == 'level1':
                    sci_hdu = pyfits.ImageHDU(data=datacube,
                                              header=self._data_header,
                                              name='SCI')
                else:
                    # Ensure the data units are preserved when saving to
                    # FITSWriter format.
                    if self._data_header is not None:
                        if 'BUNIT' in self._data_header:
                            bunit = self._data_header['BUNIT']
                            header["BUNIT"] = (bunit, "Data units")
                    primary_hdu = pyfits.PrimaryHDU(data=datacube,
                                                    header=header)
                
                # Tidy up, since datacube could be quite large.
                del datacube
            elif datashape == "slope":
                # Take a copy of the SCI data and reshape it into a cube
                dataslope = self.slope_data()
                # TODO: Add World Coordinates?
                if fileformat == 'level1':
                    sci_hdu = pyfits.ImageHDU(data=dataslope,
                                              header=self._data_header,
                                              name='SCI')
                else:
                    # Ensure the data units are preserved when saving to
                    # FITSWriter format.
                    if self._data_header is not None:
                        if 'BUNIT' in self._data_header:
                            bunit = self._data_header['BUNIT']
                            header["BUNIT"] = (bunit, "Data units")
                    primary_hdu = pyfits.PrimaryHDU(data=dataslope,
                                                    header=header)
                
                # Tidy up, since dataslope could be quite large.
                del dataslope
            else:
                # Assume hypercube.
                # TODO: Add World Coordinates?
                if fileformat == 'level1':
                    sci_hdu = pyfits.ImageHDU(data=self.data,
                                              header=self._data_header,
                                              name='SCI')
                else:
                    # Ensure the data units are preserved when saving to
                    # FITSWriter format.
                    if self._data_header is not None:
                        if 'BUNIT' in self._data_header:
                            bunit = self._data_header['BUNIT']
                            header["BUNIT"] = (bunit, "Data units")
                        else:
                            header["BUNIT"] = ('DN', "Data units")
                    primary_hdu = pyfits.PrimaryHDU(data=self.data,
                                                    header=header)
        elif fileformat == 'level2':
            # Make a copy of the SCI array which is converted to slope
            # data.
            datacopy = self.slope_data(diff_only=True)
            sci_hdu = pyfits.ImageHDU(data=datacopy,
                                      header=self._data_header,
                                      name='SCI')
        else:
            # Make a copy of the SCI_array which is reduced in size by
            # averaging the groups and/or integrations.
            datacopy = self._average_data(recalculate=recalculate)         
            if datashape == "cube":
                # Reshape the new array into a cube
                datacopy.shape = (self.nints_file*self.ngroups_file, 
                                  self.rows, self.columns)
                # TODO: Add World Coordinates?
                if fileformat == 'level1':
                    sci_hdu = pyfits.ImageHDU(data=datacopy,
                                              header=self._data_header,
                                              name='SCI')
                else:
                    # Ensure the data units are preserved when saving to
                    # FITSWriter format.
                    if self._data_header is not None:
                        if 'BUNIT' in self._data_header:
                            bunit = self._data_header['BUNIT']
                            header["BUNIT"] = (bunit, "Data units")
                    primary_hdu = pyfits.PrimaryHDU(data=datacopy,
                                                    header=header)
            else:
                # Assume hypercube.
                # TODO: Add World Coordinates?
                if fileformat == 'level1':
                    sci_hdu = pyfits.ImageHDU(data=datacopy,
                                              header=self._data_header,
                                              name='SCI')
                else:
                    # Ensure the data units are preserved when saving to
                    # FITSWriter format.
                    if self._data_header is not None:
                        if 'BUNIT' in self._data_header:
                            bunit = self._data_header['BUNIT']
                            header["BUNIT"] = (bunit, "Data units")
                    primary_hdu = pyfits.PrimaryHDU(data=datacopy,
                                                    header=header)
            
        if fileformat == 'level1' or fileformat == 'level2':
            hdulist = pyfits.HDUList([primary_hdu, sci_hdu])
        else:
            hdulist = pyfits.HDUList([primary_hdu])
        
        try:
            hdulist.writeto(filename, overwrite=overwrite)
            self._written = True
        except Exception as e:
            # If the file could not be written re-raise the exception
            # with a more meaningful error message.
            strg = "%s: Could not write exposure data file.\n   %s" % \
                (e.__class__.__name__, e)
            raise IOError(strg)
        finally:
            # Always close the file.
            hdulist.close()

    def slope_data(self, grptime=1.0, diff_only=False):
        """
        
        Generate a copy of the SCIdata array where all groups have
        been combined to make slope data.

        NOTE: This function is very inefficient!
        
        :Parameters:
        
        grptime: float, optional, default=2.785s
            The time interval between groups (which determines the
            time axis for the slope calculation).
        diff_only: bool, optional, default=False
            If True, implement a quick and dirty estimate of the
            slope by subtracting the last frame from the first.
           
        :Returned:
        
        slope_data: array_like float32
            The SCIdata array converted to slopes.
            The is a slope image for each integration.
            It can be 2-D or 3-D.
            
        """
        assert self.ngroups > 1
        output = np.zeros([self.nints, self.rows, self.columns],
                            dtype=np.float32)
        
        if diff_only:
            # Quick and dirty estimate which subtracts the last
            # ramp from the first
            timediff = grptime * (self.ngroups - 1)
            output = (self.data[:,-1,:,:] - self.data[:,0,:,:]) / float(timediff)
        else:
            # Full straight line fit
            timearray = grptime * np.array( range(0, self.ngroups) )
            for intg in range(0, self.nints):
                for row in range(0, self.rows):
                    for column in range(0, self.columns):
                        # Ouch! Slope calculated one (row,column) at a time.
                        (slope, ic) = linear_regression( timearray,
                                        self.data[intg,:,row,column])
                        output[intg,row,column] = slope
        return output

    def _average_data(self, recalculate=False):
        """
        
        Generate a copy of the SCIdata array which has been reduced
        in size by averaging integrations and groups, as prescribed
        by the intavg and grpavg parameters.
        
        :Parameters:
        
        recalculate: boolean, optional, default=False
            If True, recalculate and overwrite any existing set of
            averaged data.
            
        :Returned:
        
        average_data: array_like (float32)
            The averaged SCIdata array.
            
        """
        if self.data is None:
            raise TypeError("No exposure data to be averaged.")
        output = np.zeros([self.nints_file, self.ngroups_file,
                            self.rows, self.columns],
                            dtype=np.float32)
        count = np.zeros([self.nints_file, self.ngroups_file,
                            1, 1], dtype=np.float32)
        
        for intg in range(0, self.nints):
            intg_file = intg // self.intavg
            for grp in range(0, self.ngroups):
                grp_file = grp // self.grpavg
                output[intg_file,grp_file,:,:] += \
                    self.data[intg,grp,:,:]
                count[intg_file,grp_file,0,0] += 1.0

        # Avoid divide by zero.
        iszero = np.where(count < 1.0)
        if iszero:
            count[iszero] = 1.0
        output = output / count
        del count
        # NOTE: The new MIRI exposure model uses floating point data, so the
        # average is no longer converted to integer.
        #return output.astype(np.uint32)
        return output

    @property
    def data_averaged(self):
        # Generate the averaged data on the fly. This ensures the
        # averaging is always up to date with the latest data array.
        if self.data is not None and (self.grpavg > 1 or self.intavg > 1):
            # Generate the averaged data if it is not available.
            if self._data_averaged is None:
                self._data_averaged = self._average_data()
            return self._data_averaged
        else:
            # No averaging. Just return a copy of the data.
            return self.data

    def __str__(self):
        """
        
        Returns a string describing the exposure data.
        
        """
        strg = "Exposure data containing %d rows, %d columns, " % \
            (self.rows, self.columns)
        strg += "%d groups and %d integrations" % (self.ngroups, self.nints)
        if self.grpavg > 1 or self.intavg > 1:
            strg += "\n   with"
        if self.grpavg > 1:
            strg += " every %d groups averaged to make %d" % \
                (self.grpavg, self.ngroups_file)
        if self.intavg > 1:
            if self.grpavg > 1:
                strg += " and"
            strg += " every %d integrations averaged to make %d" % \
                (self.intavg, self.nints_file)
        strg += "."
        if self.data is not None:
            strg += "\n   The SCI_data array range is %d to %d." % \
                (self.data.min(), self.data.max())
        else:
            strg += "\n   The SCI_data array is not defined."
        if self._written:
            strg += "\n   Data have been written to a file"
        else:
            strg += "\n   Data have not yet been written to a file."
            
        return strg
    
    def statistics(self, integration=None, group=-1,
                   rowstart=None, rowstop=None, colstart=None, colstop=None,
                   clip=3.0, verbose=False):
        """
        
        Returns a string describing the mean, standard deviation and other
        useful statistics for a portion of the exposure data.
        
        :Parameters:
        
        integration: int, optional, default=all integrationsn
            The integration to be analysed. If None, all integrations
            are analysed. -1 means last integration.
        group: int, optional, default=last group
            The group to be analysed. If None, all groups are analyzed.
            -1 means last group.
        rowstart: int, optional, default=all rows
            The first row to be analysed.
        rowstop: int, optional, default=all rows
            The last row to be analysed.
        colstart: int, optional, default=all columns
            The first column to be analysed.
        colstop: int, optional, default=all columns
            The last column to be analysed.
        clip: float, optional, default=3.0
            The number of standard deviations outside of which to clip
            the data to calculate a second pass mean.
        verbose: bool
            Set True for verbose output.
            
        :Returned:
        
        stats: string
            String describing the exposure data statistics.
            
        """
        if self.data is None:
            return "No exposure data has been defined."
        if rowstart is None and rowstop is None and \
           colstart is None and colstop is None:
            if integration is None and group is None:
                data_frame = self.data
            elif integration is None:
                data_frame = self.data[:, group, :,:]
            elif group is None:
                data_frame = self.data[integration, :, :,:]
            else:
                data_frame = self.data[integration, group, :,:]
        else:
            if integration is None and group is None:
                data_frame = self.data[:, :, \
                                       rowstart:rowstop,colstart:colstop]
            elif integration is None:
                data_frame = self.data[:, group, \
                                       rowstart:rowstop,colstart:colstop]
            elif group is None:
                data_frame = self.data[integration, :, \
                                       rowstart:rowstop,colstart:colstop]
            else:
                data_frame = self.data[integration, group, \
                                       rowstart:rowstop,colstart:colstop]
        
        mean = np.mean(data_frame)
        std = np.std(data_frame)
        
        if std > 0.0 and clip > 0.1:
            lower = mean - (clip * std)
            upper = mean + (clip * std)
            cmean = scipy.stats.tmean(data_frame, limits=(lower,upper))
            cstd = scipy.stats.tstd(data_frame, limits=(lower,upper))
        else:
            cmean = mean
            cstd = std
        
        strg = "mean=%f, stdev=%f; %.1f-sigma clipped mean=%f, stdev=%f" % \
            (mean, std, clip, cmean, cstd)
            
        # Determine the variation of the mean across the data
        # in a ZONES x ZONES grid and report the extreme values.
        if verbose:
            ZONES = 5
            minima = np.array(list(range(0,ZONES))) / float(ZONES)
            maxima = np.array(list(range(1,ZONES+1))) / float(ZONES)
            rowminima = (minima * data_frame.shape[0]).astype(int)
            colminima = (minima * data_frame.shape[1]).astype(int)
            rowmaxima = (maxima * data_frame.shape[0]).astype(int)
            colmaxima = (maxima * data_frame.shape[1]).astype(int)
            minval = data_frame.max()
            minstd = 0.0
            maxval = data_frame.min()
            maxstd = 0.0
        
            for rowzone in range(0,ZONES):
                for colzone in range(0,ZONES):
                    rmin = rowminima[rowzone]
                    rmax = rowmaxima[rowzone]
                    cmin = colminima[colzone]
                    cmax = colmaxima[colzone]
                    mean = np.mean(data_frame[rmin:rmax, cmin:cmax])
                    std = np.std(data_frame[rmin:rmax, cmin:cmax])
                
                    if std > 0.0:
                        lower = mean - (3.0 * std)
                        upper = mean + (3.0 * std) 
                        cmean = scipy.stats.tmean(data_frame,
                                                  limits=(lower,upper))
                        cstd = scipy.stats.tstd(data_frame,
                                                limits=(lower,upper))
                    else:
                        cmean = mean
                        cstd = std

                    if cmean < minval:
                        minval = cmean
                        minstd = cstd
                    elif cmean > maxval:
                        maxval = cmean
                        maxstd = cstd
                    
            strg += "\nMinimum clipped mean found in " \
                "%dx%d zones is %f with stddev=%f" % \
                (ZONES, ZONES, minval, minstd)
            strg += "\nMaximum clipped mean found in " \
                "%dx%d zones is %f with stddev=%f" % \
                (ZONES, ZONES, maxval, maxstd)
        return strg

    def plot(self, plotstyle='IMAGE', integrations=None, groups=None,
                averaged=False, recalculate=False, clip=(0.0,1.0),
                plotsperpage=_MAX_PLOTS, description=''):
        """
        
        Plot the simulated exposure data as a series of frames for each
        integration and group. A new plot figure is created for each
        integration. One or more frames can be plotted on each figure.
        
        :Parameters:
        
        plotstyle: string, optional
            The plot style, which may be:
            
            * IMAGE - Plot each frame as an image.
            * HISTOGRAM - Plot each frame as a histogram.
            
            The default is 'IMAGE'. 
        integrations: tuple of ints, optional
            A list of the integrations to be plotted. If not specified,
            all integrations are plotted.
        groups: tuple of ints, optional
            A list of the groups to be plotted within each integration.
            If not specified, all groups are plotted.
        averaged: boolean, optional, default=False
            If True, plot what the data would look like after averaging
            groups and integrations according to grpavg and intavg.
        recalculate: boolean, optional, default=False
            If averaged=True this flag controls whether the averaged
            data are recalculated if already available.
        clip: tuple of 2 floats, optional, default=(0.0, 1.0)
            The lower and upper clipping limits expressed as a fraction
            of the data range. The default of (0.0, 1.0) displays the
            full range of the data.
        plotsperpage: int, optional
            The number of plots to be displayed on each figure.
            If not specified, the default is the maximum number
            which can sensibly be fitted. Must be greater than 0.
        description: string, optional
            Additional description to be shown on the plot, if required.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        :Raises:
     
        ValueError
            Raised if any of the parameters are out of range.
        IndexError
            Raised if an attempt is made to plot outside the
            exposure data area.
            
        """
        if self.data is None:
            raise TypeError("No exposure data to be plotted.")
        # Generate averaged data if needed. The numbers of groups
        # and integrations are reduced to ngroups_file and nints_file
        if averaged:
            nints = self.nints_file
            ngroups = self.ngroups_file
            plotdata = self._average_data(recalculate=recalculate)
        else:
            nints = self.nints
            ngroups = self.ngroups
            plotdata = self.data

        # Generate the integration and group lists. Single values are
        # converted into lists.
        if integrations is None:
            integrations = list(range(0,nints))
        elif isinstance(integrations,(int,float)):
            integrations = [int(integrations)]
        if groups is None:
            groups = list(range(0,ngroups))
        elif isinstance(groups,(int,float)):
            groups = [int(integrations)]
                            
        # Define the page layout in an almost square grid up to
        # the maximum number of plots per page.
        pltrows, pltcols, nfigs = mplt.subdivide(len(groups),
                                              maxplots=int(plotsperpage))
        # Add a colour bar if there is only one plot per page.
        withbar = (pltrows * pltcols) < 2
            
        # Step through the integrations and create a new figure
        # for each one.
        page = 0
        for intg in integrations:
            # New page - create a new pyplot figure (plt is a
            # global variable).
            subno = 1
            page += 1
            fig = mplt.new_figure(page, stitle=description)

            # One subplot per group up to a maximum of the given number
            # of plots per page.
            for grp in groups:
                if subno > plotsperpage:
                    # New page - create a new pyplot figure (plt is
                    # a global variable).
                    subno = 1
                    page += 1
                    fig = mplt.new_figure(page, stitle=description)
                # New subplot
                ax = mplt.add_subplot(fig, pltrows, pltcols, subno)
                    
                # Extract the data frame to be plotted and trap when the
                # integration and group requested are outside the data.
                try:
                    data_frame = plotdata[intg,grp,:,:]
                except IndexError as e:
                    strg = "Integration %d group %d is " % (intg+1,grp+1)
                    strg += "outside the range of the data."
                    strg += "\n   %s" % e
                    raise IndexError(strg)
                        
                if plotstyle == 'IMAGE':
                    self._plot_image(fig, ax, intg, grp, data_frame,
                                     averaged=averaged, clip=clip,
                                     withbar=withbar)
                elif plotstyle == 'HISTOGRAM':
                    self._plot_histogram(fig, ax, intg, grp, data_frame,
                                         averaged=averaged, clip=clip)
                else:
                    strg = "Unrecognised plot style: %s" % plotstyle
                    raise ValueError(strg)
                subno += 1

            strg = "Plotting simulated data for integration %d." % \
                    (intg+1)
            strg += " Close the plot window(s) to continue."
            mplt.show_plot(prompt=strg)

    def _plot_image(self, plotfig, plotaxis, intg, grp, data_frame,
                    averaged=False, clip=(0.0,1.0), cmap='hot', labels=True,
                    withbar=False):
        """
        
        Plot an image of a data frame which is labelled with the
        integration and group numbers derived (by adding 1) from
        the array incices intg and grp. The "averaged" flag indicates
        that averaged data are being plotted, so each plot shows the
        average of a range of integrations or groups. The data are
        clipped to the range specified in clip. A colour map is defined
        from cmap and a olour bar is added when withbar is True.
        
        :Parameters:
        
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
        intg: int
            Integration index.
        grp: int
            Group index.
        data_frame: array-like
            The 2-D data to be plotted.
        averaged: bool, optional
            Set True if averaged exposure data is being plotted.
        clip: tuple of 2 floats, optional
            Any clipping to be applied to the data in the form
            (lower, upper) scaled to the range 0.0-1.0. The
            default value of (0.0,1.0) will not clip the data.
        cmap: string, optional
            The matplotlib colour map to be used. Defaults to 'hot'.
        labels: bool, optional
            Set True to add axis labels. The default is True.
        withbar: bool, optional
            Set to True to add a colour bar. The default is False.

        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
        
        """                   
        # Label the data with integration and group numbers.
        if averaged:
            int1 = (intg * self.intavg) + 1
            int2 = int1 + self.intavg - 1
            grp1 = (grp * self.grpavg) + 1
            grp2 = grp1 + self.grpavg - 1
            strg = "Int %d-%d, Grp %d-%d" % (int1,int2,grp1,grp2)
        else:
            strg = "Int %d, Grp %d" % (intg+1, grp+1)

        # Plot the data as an image.
        mplt.plot_image(data_frame, plotfig=plotfig, plotaxis=plotaxis,
                        withbar=withbar, clip=clip, xlabel='Columns',
                        ylabel='Rows', title=strg)
         
    def _plot_histogram(self, plotfig, plotaxis, intg, grp, data_frame,
                        averaged=False,
                    clip=(0.0,1.0), nbins=256, labels=True):
        """
        
        Plot a histogram of a data frame which is labelled with the
        integration and group numbers derived (by adding 1) from the
        array incices intg and grp. The "averaged" flag indicates
        that averaged data are being plotted, so each plot shows the
        average of a range of integrations or groups. The histogram
        has nbins bins and is restricted to the data range specified
        in clip.
        
        :Parameters:
        
        pyplt: matplotlib pyplot object
            A top level matplotlib pyplot object.          
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
        intg: int
            Integration index.
        grp: int
            Group index.
        data_frame: array-like
            The 2-D data to be plotted.
        averaged: bool, optional
            Set True if averaged exposure data is being plotted.
        clip: tuple of 2 floats, optional
            Any clipping to be applied to the data in the form
            (lower, upper) scaled to the range 0.0-1.0. The
            default value of (0.0,1.0) will not clip the data.
        nbins: int, optional
            The number of histogram bins to be plotted. The
            default is 256.
        labels: bool, optional
            Set True to add axis labels. The default is True.

        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
                
        """
        # Label the data with integration and group numbers.
        if averaged:
            int1 = (intg * self.intavg) + 1
            int2 = int1 + self.intavg - 1
            grp1 = (grp * self.grpavg) + 1
            grp2 = grp1 + self.grpavg - 1
            strg = "Int %d-%d, Grp %d-%d" % (int1,int2,grp1,grp2)
        else:
            strg = "Int %d, Grp %d" % (intg+1, grp+1)
            
        # Plot the data as a histogram.
        mplt.plot_hist(data_frame.flatten(), nbins, plotfig=plotfig,
                       plotaxis=plotaxis, clip=clip, xlabel='Value',
                       ylabel='Frequency', title=strg, facecolor='b',
                       edgecolor='k')

    def plot_ramp(self, rows, columns, stime=1.0, tunit='', averaged=False,
                  show_ints=False, description=''):
        """
        
        Plot a ramp showing how the signal changes as a function of
        integration and group at the specified location in the data.
        
        :Parameters:
        
        rows: int or tuple of ints
            Either: The row at which the ramp is to be plotted.
            Or: A list of rows at which ramps are to be plotted.
        columns: int or tuple of its.
            Either: The column at which the ramp is to be plotted.  
            Or: A list of columns at which ramps are to be plotted.
        stime: float, optional, default=1.0
            The sample time per group in time units.
        tunit: string, optional
            The time unit.
        averaged: boolean, optional, default=False
            If True, plot what the data would look like after averaging
            groups and integrations according to grpavg and intavg.
        description: string, optional
            Additional description to be shown on the plot, if required. 
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        if self.data is None:
            raise TypeError("No exposure data to be plotted.")
        # Generate averaged data if needed. The numbers of groups
        # and integrations are reduced to ngroups_file and nints_file
        if averaged:
            nints = self.nints_file
            ngroups = self.ngroups_file
            plotdata = self._average_data()
        else:
            nints = self.nints
            ngroups = self.ngroups
            plotdata = self.data
            
        # Convert the data into a cube
#         datacube = plotdata[:,:,:,:]
        npoints = nints*ngroups
#         datacube.shape = (npoints, self.rows, self.columns)
        datacube = plotdata.reshape(npoints, self.rows, self.columns)
            
        # Generate an approximate timeline from the sample time per group.
        tmin = 0.0
        tmax = npoints * stime
        time = np.linspace(tmin, tmax, npoints)            

        # Define fixed labels
        xlabel = "Time"
        if tunit:
            xlabel += " (%s)" % tunit
        elif stime == 1.0:
            xlabel += " (groups)"
        ylabel = 'DN'

        if (isinstance(rows, (int,float)) and \
            isinstance(columns, (int,float))):
            
            # A single row and column has been provided.
            tstrg = "Ramp at (%d,%d) for %d groups and %d integrations" % \
                (rows, columns, ngroups, nints)

            # Extract the ramp from the data cube.
            ramp = datacube[:,rows,columns]
            
            # Plot the ramp as an XY plot.
            mplt.plot_xy(time, ramp, linefmt='bo', xlabel=xlabel,
                            ylabel=ylabel, title=tstrg)
        else:
            # Multiple rows and columns have been provided. The shortest
            # list will be the reference.
            if len(rows) < len(columns):
                nloc = len(rows)
            else:
                nloc = len(columns)
            
            tstrg = "Ramps for %d groups and %d integrations at: " % \
                (ngroups, nints)
            for loc in range(0, nloc):
                tstrg += "(%d,%d), " % (columns[loc], rows[loc])
            
            # Create a matplotlib figure and axis.
            fig = mplt.new_figure(1, stitle=tstrg)
            ax = mplt.add_subplot(fig, 1, 1, 1)

            linefmtlist = ['bo', 'rx', 'g+', 'ko', \
                           'bx', 'r+', 'go', 'kx', \
                           'b+', 'ro', 'gx', 'k+' ]
            
            for loc in range(0, nloc):

                # Extract the ramp from the data cube.
                ramp = datacube[:,rows[loc],columns[loc]]

                # The plots will cycle around the defined line formats,
                ifmt = loc % len(linefmtlist)
                lfmt = linefmtlist[ifmt] 

                # Plot each ramp as an XY plot overlaid on the same axis.
                mplt.plot_xy(time, ramp, plotfig=fig, plotaxis=ax,
                             linefmt=lfmt, xlabel=xlabel,
                             ylabel=ylabel, title='')
                
        strg = "Plotting ramp data for rows %s and columns %s." % \
                (str(rows), str(columns))
        strg += " Close the plot window to continue."
        mplt.show_plot(prompt=strg)


def load_exposure_data(filename):
    """
    
    Create an ExposureData object from an existing level 1 FITS file.
    NOTE: Data that has been averaged before being written to a file
    cannot be unaveraged when read back.
                 
    :Parameters:
    
    filename: string or tuple of strings
        The name(s) of the file(s) to be opened.
        
    :Raises:
    
    ValueError
        Raised if any of the parameters are out of range.
    IOError
        Raised if there is an error opening, reading or interpreting
        the data from the input file.
        
    :Returns:
    
        A new ExposureData object.
        
    """
    # Read a level 1 FITS file containing exposure data. The file must
    # contain a primary FITS header and a SCI extension.
    try:
        hdulist = pyfits.open(filename)
 
        # Read the contents of the file       
        header = hdulist[0].header
        kw = str('SCI')
        if kw in hdulist:
            # Assume level 1 FITS data
            sci_data = hdulist[str('SCI')].data
            data_header = hdulist[str('SCI')].header
        else:
            # Assume FITSWriter data
            sci_data = hdulist[0].data
            data_header = None
    except Exception as e:
        # If the file could not be read re-raise the exception
        # with a more meaningful error message.
        strg = "%s: Could not read exposure data file.\n   %s" % \
            (   e.__class__.__name__, e)
        raise IOError(strg)

    hdulist.close()

    # Bail out if the data has not been read successfully.
    if sci_data is None: return None

    rows = sci_data.shape[-2]
    columns = sci_data.shape[-1]

    missing_items = []
    if 'INTAVG' in header:
        intavg = int(header['INTAVG'])
    else:
        intavg = 1
        missing_items.append('INTAVG')
    if 'NINTS' in header:
        nints = int(header['NINTS']) / intavg
    elif 'NINT' in header:
        nints = int(header['NINT']) / intavg
    else:
        nints = 1
        missing_items.append('NINTS')
    if 'GRPAVG' in header:
        grpavg = int(header['GRPAVG'])
    else:
        grpavg = 1
        missing_items.append('GRPAVG')
    if 'NGROUPS' in header:
        ngroups = int(header['NGROUPS']) / grpavg
    elif 'NGROUP' in header:
        ngroups = int(header['NGROUP']) / grpavg
    else:
        ngroups = 1
        missing_items.append('NGROUPS')

    if 'READPATT' in header:
        readpatt = header['READPATT']
    elif 'READOUT' in header:
        readpatt = header['READOUT']
    else:
        readpatt = 'UNKNOWN'
        missing_items.append('READPATT')
    
    if missing_items:
        strg = "The following header keywords are missing: "
        for item in missing_items:
            strg += item + " "
        LOGGER.error(strg)

    if intavg > 1 or grpavg > 1:
        LOGGER.warning("The ExposureData object read from the file " \
            "will contain already averaged data.")

    # Finally, create and return an ExposureData object
    expdata = ExposureData(rows, columns, ngroups, nints, readpatt,
                           grpavg=1, intavg=1)
    expdata.set_header(header)
    expdata.set_data_header(data_header)
    expdata.set_exposure(sci_data.astype(np.uint32))
    return expdata


def get_file_header(filename):
    """
    
    Return the primary header contained in a given illumination
    map FITS file.
    
    :Parameters:
    
    filename: string or tuple of strings
        The name(s) of the file(s) to be opened.
        
    :Returns:
    
    header: pyfits Header object
        Either the header contained in the file, or None if the
        file could not be opened.
    
    """
    # Attempt to read the FITS file and obtain the primary header.
    try:
        hdulist = pyfits.open(filename)
        header = hdulist[0].header
        
    except Exception as e:
        # If the file cannot be read, return a null identification.
        strg = "Illumination map file %s cannot be read:" % filename
        strg += "\n   (%s)." % e
        LOGGER.error( strg )
        return None
    
    return header

def set_metadata(exposure_data, metadata):
    """
        
    Add exposure data keywords to the given MIRI metadata.
        
    :Parameters:
        
    metadata: dictionary-like object
        The MIRI metadata to which exposure data keywords should be
        added.
            
    :Returns:
        
    metadata: dictionary-like object
        An updated MIRI metadata object.
            
    """
    metadata["NINTS"] = exposure_data.nints
    metadata["INTAVG"] = exposure_data.intavg
    metadata["NGROUPS"] = exposure_data.ngroups
    metadata["NFRAMES"] = exposure_data.nframes
    metadata["GROUPGAP"] = exposure_data.groupgap
    metadata["GRPAVG"] = exposure_data.grpavg
    metadata['READPATT'] = exposure_data.readpatt
        
    return metadata


#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_exposure_data module.
#
if __name__ == '__main__':
    print( "Testing the LEGACY ExposureData class" )

    PLOTTING = True        # Set to False to turn off plotting.
    test_file_name = './data/SCAExposureDataTest.fits'
    
    # Create an ExposureData object for a 6 row x 8 column detector
    # and display its initial state.
    rows = 6
    columns = 8
    nints = 2
    ngroups = 2
    readpatt = 'FAST'
    exposure = ExposureData(rows, columns, ngroups, nints, readpatt)
    print( exposure )
    print( exposure.data )
    print( "Initial data: " + exposure.statistics() )

    # Generate a series of integrations and groups with predictable data values
    # and build up the exposure data.    
    for intg in range(0,nints):
        for group in range(0,ngroups):
            data = 1 * np.indices([rows,columns])[0] * (group+1) * (intg+2)
            exposure.set_group(data, group, intg)
            del data

    print( exposure )
    print( "After adding data: " + exposure.statistics() )
    if PLOTTING:
        exposure.plot(description='Test Exposure Data')
        exposure.plot(plotstyle='HISTOGRAM', description='Test Exposure Data')
        exposure.plot_ramp(3, 4, stime=2.75, tunit='seconds',
                           description='(test 1)')
    else:
        print( "*** matplotlib plotting library not available." )

    # Write the data to the output file. Allow the file to be overwritten.
    print( "Writing", test_file_name )
    exposure.save(test_file_name, overwrite=True)
    del exposure
    
    # Create another ExposureData object, only this time with a larger
    # number of nints and groups which are expected to be averaged.
    nints = 6
    intavg = 2
    ngroups = 8
    grpavg = 4
    exposure = ExposureData(rows, columns, ngroups, nints, readpatt,
                            grpavg=grpavg, intavg=intavg)
    print( exposure )
    print( "Initial larger data: " + exposure.statistics() )

    # Generate a series of integrations and groups with predictable data values
    # and build up the exposure data.    
    for intg in range(0,nints):
        for group in range(0,ngroups):
            data = 1 * np.indices([rows,columns])[0] * (group+1) * (intg+2)
            exposure.set_group(data, group, intg)
            del data

    print( exposure )
    print( "After adding data: " + exposure.statistics() )
    if PLOTTING:
        exposure.plot(averaged=True,
                         description='Test Exposure Data - Averaged')
        exposure.plot_ramp(list(range(0,6)), list(range(0,6)),
                           stime=2.785, tunit='seconds', averaged=True,
                           description='(test 2 - averaged)')
        slope_data = exposure.slope_data(grptime=2.785)
        mplt.plot_image(slope_data, xlabel='Columns', ylabel='Rows',
                        title='Test Exposure Data - Slopes')
    else:
        print( "*** matplotlib plotting library not available." )
    
    # Add MIRI metadata to the header
    header = pyfits.Header()
    metadata = Metadata()    
    metadata = set_metadata(exposure, metadata)
    metadata.to_fits_header(header)
    exposure.set_header(header)
    
    # Write the data to the output file. Allow the file to be overwritten.
    print( "Writing", test_file_name )
    exposure.save(test_file_name, overwrite=True)

    # Read the same data back from the file.
    newexp = load_exposure_data(test_file_name)
    print( "After reading back from file: " + newexp.statistics() )
    print( newexp )
    if PLOTTING:
        newexp.plot(description='Test Exposure Data - read back from file')
    else:
        print( "*** matplotlib plotting library not available." )     

    print( "Test finished." )
