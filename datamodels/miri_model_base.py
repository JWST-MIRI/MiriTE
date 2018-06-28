#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

The base model for MIRI data. An extension to the standard STScI data 
model.

:Reference:

The STScI jwst.datamodel library documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

20 Nov 2012: First created as miri_measured_model
05 Feb 2013: Features common to all MIRI data moved here, plus some
             compatibility utilities.
07 Feb 2013: Conversion of metadata into a readable string improved
             with the aid of the new fits_metadata_dict function.
05 Jun 2013: Added get_fits_comments, get_fits_history and
             find_fits_values methods.
10 Jun 2013: All FITS-specific methods renamed so their names include
             the word "fits". Methods reordered into sensible blocks.
             New functions for getting and setting FITS keywords and
             getting and setting FITS comment and history records.
13 Jun 2013: FITS query functions modified to look for a specific HDU
             when requested.
03 Jul 2013: Metadata defaults moved to the schema. maskable()
             function added.
05 Jul 2013: Check for data objects without a type.
29 Jul 2013: Corrected typo that prevented data tables from being
             listed in the __str__() function. stats() and
             get_data_stats() function added.
22 Aug 2013: get_field_names function added.
13 Sep 2013: If a data product contains only one table, the stats()
             method will list it.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
23 Jan 2014: Modified for jsonschema draft 4: Convenience functions
             in MiriDataModel modified to use the "walk_schema"
             functions so they are independent of schema structure.
04 Mar 2014: New set_housekeeping_metadata function. MODELNAM added
             to set_instrument_metadata and nints and groups checked
             for None in set_exposure_metadata.
11 Apr 2014: Modified list_data_arrays, list_data_tables, dataname_to_hduname,
             hduname_to_dataname and get_field_names functions to be
             more independent of schema structure (so they work with
             jsonschema draft 4 schemas). Added set_data_units method.
18 Jun 2014: Removed some surplus code. Some formatting corrections.
             Corrected a mix-up between detector frame time and sample time.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
30 Sep 2014: Added the convert_flags method, to aid CDP-3 delivery.
10 Oct 2014: set_table_units method added.
16 Oct 2014: Improved the logic of the set_table_units method. Some print
             statements added and commented out. (Use for debugging.)
07 Nov 2014: Change the unit comparison logic so that a units setting of
             'None' is equivalent to a null string.
09 Dec 2014: Added ability to set detector settings and reference file
             metadata.
16 Jan 2015: SUBARRAY_DICT modified so that SUBSTRTi and SUBSIZEi keywords
             are written even for full-frame data.
06 Feb 2015: Define DATE-OBS and TIME-OBS metadata separately.
04 Mar 2015: Improve the error message from set_fits_keyword.
11 Mar 2015: Metadata names bought up to date with changes to JWST core schema.
             reffile.pedigree added to housekeeping metadata.
             observation.id changed to observation.obs_id.
             group_integration_time changed to group_time.
             All metadata values are skipped if they are None.
10 Jun 2015: Use 'GENERIC' for a generic subarray rather than 'T' or 'ANY'.
             Obsolete set_fits_history and get_fits_history removed and a
             problem corrected which prevented the first COMMENT record being
             created.
03 Jul 2015: Modified get_meta_str() method so that table unit metadata can
             be excluded.
09 Jul 2015: Removed redundant set_table_units function. axes label set
             automatically to 'records' when the data object is a table.
             Automatically list the column names of a table.
06 Aug 2015: Do not automatically set the end time when exposure metadata
             are added to a data model.
31 Aug 2015: Exposure start, middle and end time are now floating point
             values rather than date-time strings.
02 Sep 2015: Metadata copying function now copies every header item
             defined in the schema.
10 Sep 2015: Corrected a bug in add_fits_comment function.
             Add a history record stating the data model used to create
             each file.
11 Sep 2015: Moved the "plot" method here, to reduce the number of
             duplications of this same method. It will fail gracefully
             if given an unplottable data model 
07 Oct 2015: Made exception catching Python 3 compatible.
22 Oct 2015: Missing metadata attributes will now raise an AttributeError
             exception.
08 Dec 2015: Correction to subarray parameters to match Mike Ressler's PASP
             paper.
25 Jan 2016: Changes for compatibility with ASDF version of jwst.datamodels:
             Try imports from models.fits and models.fits_support.
             Date objects converted to string before writing to metadata.
28 Jan 2016: Modified to meet the new requirements of the ASDF release of
             the jwst.datamodels package. The API of the callback function for
             walk_schema has changed from (results, schema, path) to
             (subschema, path, combiner, ctx, recurse). Schema paths are
             managed as string lists rather than dot-separated strings.
29 Jan 2016: Renamed "title" attribute to "_subtitle" to avoid clashing
             with schema components.
10 Mar 2016: copy_metadata method now uses the source data model's metadata
             as a baseline rather than the destination data model.
             Ignore string list added to copy_metadata method so that
             irrelevant metadata can be skipped without generating warnings.
05 Apr 2016: Work around jwst.datamodels problem by copying self.meta.observation.date
             to self.meta.date. Renamed 'dtype' to 'datatype in query functions.
             Search for data entries by looking for type=None rather than
             type='data'.
06 Apr 2016: Change for compatibility with ASDF version of jwst.datamodels: 
             the "add_history" method now wraps the string in an 
             asdf.tags.core.HistoryEntry object, including a timestamp, before
             adding this HistoryEntry to the list.
09 May 2016: Moved get_exposure_type and get_exp_type to here from
             miri_exposure_model.py.
08 Jun 2016: Corrected a typo in the setting of the exposure type.
             Also changed the logic so that 'MIR_DARK' can be defined even
             for an unknown detector.
19 Jul 2016: Added set_pointing and set_wcs_metadata functions.
31 Aug 2016: Added add_referencefile_history function. Added options to
             set the USEAFTER function to today's date or to a default date.
02 Sep 2016: Added add_schema_mapping function and used it to obtain a
             resolver for accessing the JWST data model schemas from
             within the MIRI data model schemas.
20 Sep 2016: Resolver parameter converted to ASDF extension.
29 Sep 2016: ASDF extension added for accessing the MIRI data models.
30 Sep 2016: Extra World Coordinates keywords added.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. Added "date" parameter to
             housekeeping metadata.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
12 Sep 2017: World Coordinates extended to 4 axes. Added warnings when an
             attempt is made to write more axes than defined in the data model.
15 Sep 2017: World Coordinates setting functions split into smaller parts.
             "get_" functions added to make accessing the metadata easier.
22 Sep 2017: Added more WCS-related metadata to the set_wcs_metadata function.
20 Nov 2017: Corrected a bug in the set_data_units method.
29 Nov 2017: Updated get_history and get_history_str methods to cope with
             a change to the data type of the .history attribute in the
             underlying JWST data model.
04 Jan 2018: Fill a masked array before calculating statistics in
             get_data_stats() to prevent a numpy warning.
19 Jan 2018: Added set_telescope function for creating test data.
24 Jan 2018: Added association metadata to set_observation_metadata.
11 Apr 2018: FLENS added to filter list and "N/A" added to some options in
             data model metadata.
17 May 2018: Python 3: Converted dictionary keys return into a list.
24 May 2018: Modified to work with history records stored either as a list
             of strings or a dictionary (to cope with a change to the JWST
             library).
28 Jun 2018: Modified to work with history records stored either as a dictionary
             or a HistoryList object (to cope with another change to the JWST
             library). Added the get_title_and_metadata() function to display
             the history records along with metadata.

@author: Steven Beard (UKATC), Vincent Geers (UKATC)

"""
# This module is now converted to Python 3.


import os
import datetime
import warnings
import numpy as np
from astropy.extern import six
from astropy.time import Time
import asdf.util
from asdf.tags.core import HistoryEntry
from asdf.extension import AsdfExtension
#import numpy.ma as ma

# Import the STScI base data model and utilities
import jwst.datamodels
import jwst.datamodels.fits_support as mfits
import jwst.datamodels.schema as mschema
from jwst.datamodels.model_base import DataModel
from jwst.datamodels.history import HistoryList

# Import the MIRI data models and the data model plotter.
import miri.datamodels
from miri.datamodels.plotting import DataModelPlotVisitor


# List all classes and global functions here.
__all__ = ['get_exp_type', 'MiriDataModel']

#
# Dictionary of available subarray options. Each tuple contains
# (first column, first row, number of columns, number of rows).
# This dictionary is used by the set_subarray_metadata convenience
# function.
# TODO: Should this table be moved to a MIRI parameters file?
#
SUBARRAY_DICT = {}
SUBARRAY_DICT['FULL'] =          (   1,   1, 1032, 1024 )
SUBARRAY_DICT['GENERIC'] =       (   1,   1, 1032, 1024 )
SUBARRAY_DICT['MASK1065'] =      (   1,  19,  288,  224 )
SUBARRAY_DICT['MASK1140'] =      (   1, 245,  288,  224 )
SUBARRAY_DICT['MASK1550'] =      (   1, 467,  288,  224 )
SUBARRAY_DICT['MASKLYOT'] =      (   1, 717,  320,  304 )
SUBARRAY_DICT['BRIGHTSKY'] =     ( 457,  51,  512,  512 )
SUBARRAY_DICT['SUB256'] =        ( 413,  51,  256,  256 )
SUBARRAY_DICT['SUB128'] =        (   1, 889,  136,  128 )
SUBARRAY_DICT['SUB64'] =         (   1, 779,   72,   64 )
SUBARRAY_DICT['SLITLESSPRISM'] = (   1, 529,   72,  416 )

# Pre-BURST-mode versions of the subarrays
#SUBARRAY_DICT['BRIGHTSKY'] =     (   1,  37,  968,  512 )
#SUBARRAY_DICT['SUB256'] =        (   1,  37,  668,  256 )

#
# Public global function.
#
def get_exp_type(detector, filter, subarray='FULL', datatype='SCIENCE'):
    """
    
    Helper function which derives an exposure type given a combination
    of detector name, filter, subarray and data type
    
    NOTE: This function needs to know the data type to return the
    correct exposure type for MIRI calibration data.
    
    :Parameters:
    
    detector: str
        Name of detector (MIRIMAGE, MIRIFUSHORT or MIRIFULONG)
    filter: str
        Name of filter
    subarray: str, optional
        Name of subarray (defaults to 'FULL')
    datatype: str, optional
        Known data type (SCIENCE, DARK, FLAT or TARGET)
        Defaults to SCIENCE.
        
    :Returns:
    
    exp_type: str
        The exposure type string (EXP_TYPE)
    
    """
    assert isinstance(detector, str)
    assert isinstance(filter, str)
    
    if 'MIRIFU' in detector:
        if (filter and 'OPAQUE' in filter) or (datatype == 'DARK'):
            exp_type = 'MIR_DARK'
        elif datatype == 'FLAT':
            exp_type = 'MIR_FLAT-MRS'
        else:
            exp_type = 'MIR_MRS'
    elif 'IM' in detector:
        if (filter and 'OPAQUE' in filter) or (datatype == 'DARK'):
            exp_type = 'MIR_DARK'
        elif datatype == 'FLAT':
            exp_type = 'MIR_FLAT-IMAGE'
        elif filter and 'P750L' in filter:
            if 'SLITLESS' in subarray:
                exp_type = 'MIR_LRS-SLITLESS'
            else:
                exp_type = 'MIR_LRS-FIXEDSLIT'
        elif subarray and 'LYOT' in subarray:
            if datatype == 'SCIENCE':
                exp_type = 'MIR_LYOT'
            elif datatype == 'TARGET':
                exp_type = 'MIR_TACQ'
            else:
                exp_type = 'MIR_CORONCAL'
        elif subarray and 'MASK' in subarray:
            if datatype == 'SCIENCE':
                exp_type = 'MIR_4QPM'
            elif datatype == 'TARGET':
                exp_type = 'MIR_TACQ'
            else:
                exp_type = 'MIR_CORONCAL'
        else:
            exp_type = 'MIR_IMAGE'
    elif datatype == 'DARK':
        exp_type = 'MIR_DARK'
    else:
        # Exposure type cannot be derived
        exp_type = ''
    return exp_type

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

def _shape_to_string(shape, axes=None):
    """
    
    Return an annotated string describing the given data shape.
    
    For example, if shape=(3,4) and axes=('columns', 'rows')
    the string "3 rows x 4 columns" will be returned. If axes
    is not given, the string "3 x 4" will be returned.
    
    If shape=(2,3,4) and axes=('hypercubes','cubes','rows','columns')
    the string "2 cubes x 3 rows x 4 columns" will be returned. The
    axes array is truncated on the left to match.
        
    :Attributes:
        
    shape: tuple of ints
        The shape to be described.
    axes: tuple of str, optional
        The names of the axis labels corresponding to the given shape.
        If the axes list is longer that the shape it will be truncated
        on the left.
        If the axes list is not given, or is too short, the shape will
        be described without labels.
            
    :Returns:
        
    description: str
        A human readable description of the given shape.
        For example, "3 cubes x 4 rows x 6 columns".
    
    """
    if axes is not None:
        if len(shape) == 1:
            strg = "%d %s" % (int(shape[0]), axes[-1])  
        elif len(axes) >= len(shape):
            strg = ''
            iishift = len(axes) - len(shape)
            for ii in range(0,len(shape)):
                jj = ii + iishift
                if ii < len(shape) - 1:
                    strg += "%d %s x " % (int(shape[ii]), axes[jj])
                else:
                    strg += "%d %s" % (int(shape[ii]), axes[jj])
        else:
            # A much simpler implementation when there are inadequate
            # axis labels
            strg = ' x '.join(str(s) for s in shape)
    else:
        # A much simpler implementation when there are no axis labels
        strg = ' x '.join(str(s) for s in shape)
    return strg

# THIS CLASS HAS NOW BEEN IMPLEMENTED IN JWST.DATAMODELS
# class JWSTExtension(AsdfExtension):
#     """
#     
#     A class which defines an extension to the ASDF URI mapping
#     which recognises the URI 'http://jwst.stsci.edu/schemas/'
#     as a reference to the schemas belonging to the jwst.datamodels
#     package.
#     
#     """
#     @property
#     def types(self):
#         return []
# 
#     @property
#     def tag_mapping(self):
#         return []
#
#     @property
#     def url_mapping(self):
#         # Define the uri and directory containing the additional schemas
#         schema_uri_base = 'http://jwst.stsci.edu/schemas/'
#         schema_path = os.path.abspath(os.path.join( \
#                             os.path.dirname(jwst.datamodels.__file__),
#                             'schemas'))
#         return [(schema_uri_base,
#                  asdf.util.filepath_to_url(schema_path) +
#                  '/{url_suffix}')]


class MIRIExtension(AsdfExtension):
    """
    
    A class which defines an extension to the ASDF URI mapping
    which recognises the URI 'http://miri.stsci.edu/schemas/'
    as a reference to the schemas belonging to the miri.datamodels
    package.
    
    """
    @property
    def types(self):
        return []

    @property
    def tag_mapping(self):
        return []

    @property
    def url_mapping(self):
        # Define the uri and directory containing the additional schemas
        schema_uri_base = 'http://miri.stsci.edu/schemas/'
        schema_path = os.path.abspath(os.path.join( \
                            os.path.dirname(miri.datamodels.__file__),
                            'schemas'))
        return [(schema_uri_base,
                 asdf.util.filepath_to_url(schema_path) +
                 '/{url_suffix}')]


class MiriDataModel(DataModel):
    """
    
    A base data model for MIRI data, with MIRI-specific utilities and
    add-ons, based on the STScI base model, DataModel. The MiriDataModel
    also includes MIRI-specific additions to the metadata.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided later, the shape is derived from the data array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.

    title: str, optional
        An additional descriptive title for the data structure (if the
        default title contained in the schema is too vague).
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_core.schema.yaml"

    def __init__(self, init=None, title='', extensions=None, **kwargs):
        """
        
        Initialises the MiriDataModel class.
        
        Parameters: See class doc string.

        """
        # It is necessary to handle the "extensions" argument because it
        # will be defined when a data model is copied.
        # When a data model is first created, extensions will be None.
        # in which case the MIRI extension needs to be defined.
        # When a data model is copied, extensions will contain the current
        # list of extensions, in which case the list needs to be passed
        # unchanged to the parent class.
        if extensions is None:
#            extensions = [MIRIExtension(), JWSTExtension()]
           extensions = [MIRIExtension()]

        # Initialise the underlying STScI data model.
        super(MiriDataModel, self).__init__(init=init, extensions=extensions,
                                            **kwargs)

        # Initialise the observation date if not already defined.
        if hasattr(self, 'meta'):
            if hasattr(self.meta, 'observation') and hasattr(self.meta, 'date'):
                if self.meta.date and not self.meta.observation.date:
                    self.meta.observation.date = self.meta.date
        
        # Initialise the metadata cache used by the schema search functions.
        self._metadata_cache = {}
        
        # Set the data product subtitle and history, if provided.
        self._subtitle = title
        created_strg = "Created from: %s" % self.__class__.__name__
        self.add_unique_history(created_strg)

    #
    # Convenience functions for setting commonly associated blocks
    # of MIRI metadata. All these functions assume the metadata is
    # as defined in the "miri_core.schema.yaml" file.
    #
    def set_housekeeping_metadata(self, origin, author=None, pedigree=None,
                                  version=None, date=None, useafter=None,
                                  description=None):
        """
        
        Convenience function to define housekeeping metadata.
        
        Useful for setting up test data.

        NOTE: The author, pedigree, useafter and descroption parameters are
        only valid for data models which define the reference file metadata
        in their data model.
                
        :Parameters:
        
        origin: str
            Organisation responsible for creating the file.
        author: str, optional
            Author of the data.
        pedigree: str, optional
            Pedigree of a reference file.
        version: str, optional
            Version number of data in file
        date: str, optional
            Date when a data model was created.
            'TODAY' is translated into 'YYYY-MM-DDT00:00:00', where
            YY-MM-DD is today's date.
        useafter: str, optional
            Date after which a reference file is valid.
            'DEFAULT' is translated into '2000-01-01T00:00:00'.
            'TODAY' is translated into 'YYYY-MM-DDT00:00:00', where
            YY-MM-DD is today's date.
        description: str, optional
            Description of reference file.
            
        """
        if hasattr(self, 'meta'):
            missing = ''
            if origin is not None:
                self.meta.origin = origin
            if version is not None:
                self.meta.version = version

            if date == 'TODAY':
                self.meta.date = str(datetime.date.today()) + \
                    'T00:00:00'
            elif date is not None:
                self.meta.date = date
                    
            if pedigree is not None or author is not None or \
               useafter is not None or description is not None:
                self.set_referencefile_metadata(reftype=None, author=author,
                                                pedigree=pedigree,
                                                useafter=useafter,
                                                description=description)
        else:
            strg = "Metadata attributes missing from data model"
            raise AttributeError(strg)

    def set_referencefile_metadata(self, reftype=None, author=None,
                                   pedigree=None, useafter=None,
                                   description=None):
        """
        
        Convenience function to define reference file metadata.
        
        NOTE: This function is only valid for data models which define
        the reference file metadata in their data model.
        
        :Parameters:

        reftype: str, optional
            Type of reference data contained.
        author: str, optional
            Author of the reference data.
        pedigree: str, optional
            Pedigree of a reference file.
        useafter: str, optional
            Date after which a reference file is valid
            'DEFAULT' is translated into '2000-01-01T00:00:00'.
            'TODAY' is translated into 'YYYY-MM-DDT00:00:00', where
            YY-MM-DD is today's date.
        description: str, optional
            Description of reference file.
            
        """
        if hasattr(self, 'meta'):
            if reftype is not None:
                self.meta.type = reftype
            if pedigree is not None:
                self.meta.pedigree = pedigree
            if author is not None:
                self.meta.author = author
            if useafter == 'DEFAULT':
                self.meta.useafter = '2000-01-01T00:00:00'
            elif useafter == 'TODAY':
                self.meta.useafter = str(datetime.date.today()) + \
                    'T00:00:00'
            elif useafter is not None:
                self.meta.useafter = useafter
            if description is not None:
                self.meta.description = description
        else:
            strg = "Reference file metadata attributes missing from data model"
            raise AttributeError(strg)

    def add_referencefile_history(self, document=None, software=None,
                                   dataused=None, differences=None,
                                   history=None):
        """
        
        Convenience function to add required HISTORY strings to
        reference file metadata.
        
        :Parameters:

        document: str, optional
            DOCUMENT string describing name of document.
        software: str, optional
            SOFTWARE string decsribing name of software.
        dataused: str, optional
            DATA USED string describing data used.
        differences: str, optional
            DIFFERENCES string describing differences.
        history: list of str, optional
            Optional list of additional history strings.
            
        """
        if document is not None:
            document_strg = 'DOCUMENT: ' + str(document)
            self.add_history(document_strg)
        if software is not None:
            software_strg = 'SOFTWARE: ' + str(software)
            self.add_history(software_strg)
        if dataused is not None:
            dataused_strg = 'DATAUSED: ' + str(dataused)
            self.add_history(dataused_strg)
        if differences is not None:
            differences_strg = 'DIFFERENCES: ' + str(differences)
            self.add_history(differences_strg)
        if history is not None:
            if isinstance(history, (tuple,list)):
                for hstrg in history:
                    self.add_history(str(hstrg))
            else:
                self.add_history(str(history))
                    
    def set_observation_metadata(self, dateobs=None, timeobs=None, obsid=None,
                                 obsnumber=None, programnumber=None,
                                 visitid=None, visitnumber=None, visitgroup=None,
                                 activityid=None, exposurenumber=None):
        """
        
        Convenience function to define observation and association metadata.
        By default it sets the observation date to the current date and time
        and leaves the other fields empty.
        Useful for setting up test data.
        
        :Parameters:
        
        dateobs: date object, optional
            Observation date, formatted into a Python date object or string.
            If not given the current date will be used.
        timeobs: time object, optional
            Observation time, formatted into a Python time object or string.
            If not given the current time will be used.
        obsid: str, optional
            Observation ID.
        obsnumber: number, optional
            Observation number
        programnumber: str, optional
            Program number
        visitid: str, optional
            Visit ID
        visitnumber: str, optional
            Visit number
        visitgroup: str, optional
            Visit group
        activityid: str, optional
            Activity ID
        exposurenumber: str, optional
            Exposure number
            
        """
        from datetime import datetime
        if hasattr(self, 'meta') and hasattr(self.meta, 'observation'):
            if dateobs is None:
                dt = datetime.now()
                dateobs = dt.date()
                self.meta.observation.date = str(dateobs)
            else:
                self.meta.observation.date = str(dateobs)
            if timeobs is None:
                dt = datetime.now()
                timeobs = dt.time()
                self.meta.observation.time = str(timeobs)
            else:
                self.meta.observation.time = str(timeobs)
            if obsid is not None:
                self.meta.observation.obs_id = str(obsid)
            if obsnumber is not None:
                self.meta.observation.observation_number = str(obsnumber)
            if programnumber is not None:
                self.meta.observation.program_number  = str(programnumber)
            if visitid is not None:
                self.meta.observation.visit_id = str(visitid)
            if visitnumber is not None:
                self.meta.observation.visit_number = str(visitnumber)
            if visitgroup is not None:
                self.meta.observation.visit_group = str(visitgroup)
            if activityid is not None:
                self.meta.observation.activity_id  = str(activityid)
            if exposurenumber is not None:
                self.meta.observation.exposure_number  = str(exposurenumber)
        else:
            strg = "Observation metadata attributes missing from data model"
            raise AttributeError(strg)

    def set_target_metadata(self, ra, dec):
        """
        
        Convenience function to define target metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        ra: number
            Target Right Ascension.
        dec: number
            Target declination.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'target'):
            if ra is not None:
                self.meta.target.ra = ra
            if dec is not None:
                self.meta.target.dec = dec
        else:
            strg = "Target metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_target_metadata(self):
        """
        
        A trivial helper function which returns the target metadata.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (ra, dec):
            Target Right Ascension and Declination.
            (None, None) if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'target'):
            return (self.meta.target.ra, self.meta.target.dec)
        else:
            return (None, None)

    def set_pointing_metadata(self, ra_v1, dec_v1, pa_v3):
        """
        
        Convenience function to define telescope pointing metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        ra_v1: number
            Telescope pointing V1 axis Right Ascension.
        dec_v1: number
            Telescope pointing V1 axis declination.
        pa_v3: number
            Telescope pointing V3 axis position angle.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'pointing'):
            if ra_v1 is not None:
                self.meta.pointing.ra_v1 = ra_v1
            if dec_v1 is not None:
                self.meta.pointing.dec_v1 = dec_v1
            if pa_v3 is not None:
                self.meta.pointing.pa_v3 = pa_v3                
        else:
            strg = "Pointing metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_pointing_metadata(self):
        """
        
        A trivial helper function which returns the pointing metadata.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (ra_v1, dec_v1, pa_v3):
            Telescope pointing V1 axis Right Ascension and Declination
            and Telescope pointing V3 axis position angle.
            (None, None, None) if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'pointing'):
            return (self.meta.pointing.ra_v1, self.meta.pointing.dec_v1,
                    self.meta.pointing.pa_v3)
        else:
            return (None, None, None)

    def set_telescope(self):
        """
        
        A convenience function to ensure the telescope metadata is set
        correctly. The telescope should always be 'JWST'.
        
        :Parameters:
        
        None
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'telescope'):
            self.meta.telescope = 'JWST'

    def set_instrument_metadata(self, detector, modelnam='FM', detsetng='ANY',
                                filt='', channel='', band='',
                                ccc_pos='OPEN', deck_temperature=None,
                                detector_temperature=None):
        """
        
        Convenience function to define instrument metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        detector: str
            Name of detector.
        modelnam: str
            Name of instrument model (if not 'FM')
        detsetng: str
            Detector settings used ('RAL1', 'JPL1' or 'ANY')
        filt: str, optional
            Name of instrument filter (if any).
        channel: str, optional
            Name of instrument channel (if any)
        band: str, optional
            Name of instrument band (if any)
        ccc_pos: str, optional
            MIRI CCC position (if not 'OPEN')
        deck_temperature: number, optional
            MIRI deck temperature in K.
        detector_temperature: number, optional
            MIRI detector temperature in K.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'instrument'):
            self.meta.instrument.name = 'MIRI'
            if detector is not None:
                self.meta.instrument.detector = detector
            if modelnam:
                self.meta.instrument.model = modelnam
            if detsetng:
                self.meta.instrument.detector_settings = detsetng
            if filt:
                self.meta.instrument.filter = filt
            if channel:
                self.meta.instrument.channel = channel
            if band:
                self.meta.instrument.band = band
            if ccc_pos:
                self.meta.instrument.ccc_pos = ccc_pos
            if deck_temperature is not None:
                self.meta.instrument.deck_temperature = deck_temperature
            if detector_temperature is not None:
                self.meta.instrument.detector_temperature = detector_temperature
        else:
            strg = "Instrument metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_instrument_detector(self):
        """
        
        A trivial helper function which returns the instrument detector name,
        model and settings.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (detector, modelnam, detsetng):
            Instrument detector name, model and settings.
            (None, None, None) if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'instrument'):
            return (self.meta.instrument.detector, self.meta.instrument.model,
                    self.meta.instrument.detector_settings)
        else:
            return (None, None, None)

    def get_instrument_filter(self):
        """
        
        A trivial helper function which returns the instrument filter name.
        Only valid for MIRI imager data. None will be returned for MRS
        data.
        
        :Parameters:
        
        None
        
        :Returns:
        
        filt: str
            Name of instrument filter. None if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'instrument'):
            if hasattr(self.meta.instrument, 'filter'):
                return self.meta.instrument.filter
            else:
                return None
        else:
            return None

    def get_instrument_channel_band(self):
        """
        
        A trivial helper function which returns the instrument channel
        and band. Only valid for MIRI MRS data. None will be returned
        for imager data.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (channel, band)
            Name of instrument channel and band.
            (None, None) if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'instrument'):
            if hasattr(self.meta.instrument, 'channel') and hasattr(self.meta.instrument, 'band'):
                return (self.meta.instrument.channel, self.meta.instrument.band)
            else:
                return (None, None)
        else:
                return (None, None)

    def set_exposure_metadata(self, readpatt, nints, ngroups, nframes=1,
                              sample_time=10.0, frame_time=None,
                              integration_time=None, exposure_time=None,
                              group_time=None, groupgap=0, grpavg=1, intavg=1,
                              reset_time=0, frame_resets=3, start_time='NOW',
                              end_time=None):
        """
        
        Convenience function to define exposure metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        readpatt: str
            Name of detector readout pattern.
        nints: int
            Number of integrations per exposure (missed out if set to None).
        ngroups: int
            Number of groups per integration (missed out if set to None).
        nframes: int, optional
            Number of frames coadded in group (always 1 for MIRI data)
        sample_time: number, optional
            Detector sample time in microseconds (default 10).
        frame_time: number, optional
            Detector frame time in seconds
        integration_time: number, optional
            Integration time.
        exposure_time: number, optional
            Exposure time.
        group_time: number, optional
            Time between groups.
        groupgap: int, optional
            Number of dropped frames between groups (default 0).
        grpavg: int, optional
            Number of groups averaged (default 1).
        intavg: int, optional.
            Number of integrations averaged (default 1).
        reset_time: number, optional
            Detector reset time.
        frame_resets: int, optional
            Number of detector frame resets (if not 3).
        start_time: str or float, optional
            Date/time of start of exposure. Strings other than 'NOW' are
            converted to floating point. If set to 'NOW', the current
            date-time is used. Skipped if set to None. By default
            this is set to the current date and time.
        end_time: str or float, optional
            Date/time of start of exposure. Strings other than 'NOW' are
            converted to floating point. If set to 'NOW', the current
            date-time is used. By default this is not set.
        
        """
        import time
        if hasattr(self, 'meta') and hasattr(self.meta, 'exposure'):
            if readpatt is not None:
                self.meta.exposure.readpatt = readpatt
            if nints is not None:
                self.meta.exposure.nints = nints
            if ngroups is not None:
                self.meta.exposure.ngroups = ngroups
            if nframes is not None:
                self.meta.exposure.nframes = nframes
            if sample_time is not None:
                self.meta.exposure.sample_time = sample_time
            if frame_time is not None:
                self.meta.exposure.frame_time = frame_time
            if integration_time is not None:
                self.meta.exposure.integration_time = integration_time
            if exposure_time is not None:
                self.meta.exposure.exposure_time = exposure_time
            if group_time is not None:
                self.meta.exposure.group_time = \
                    group_time
            if groupgap is not None:
                self.meta.exposure.groupgap = groupgap
            if grpavg is not None:
                self.meta.exposure.groups_averaged = grpavg
            if intavg is not None:
                self.meta.exposure.integrations_averaged = intavg
            if reset_time is not None:
                self.meta.exposure.reset_time = reset_time
            if frame_resets is not None:
                self.meta.exposure.frame_resets = frame_resets
            if start_time == 'NOW':
                start_time = time.time()
            if end_time == 'NOW':
                end_time = time.time()
            if start_time is not None:
                self.meta.exposure.start_time = float(start_time)
            if end_time is not None:
                self.meta.exposure.end_time = float(end_time)
        else:
            strg = "Exposure metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_exposure_readout(self):
        """
        
        A trivial helper function which returns the exposure readout
        metadata.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (readpatt, nints, ngroups, nframes):
            Detector readout pattern and the number of integrations,
            groups and frames.
            (None, None, None, None) if not defined.  
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'exposure'):
            return (self.meta.exposure.readpatt, self.meta.exposure.nints,
                    self.meta.exposure.ngroups, self.meta.exposure.nframes)
        else:
            return (None, None, None, None)

    def set_exposure_type(self, detector=None, filter=None, subarray=None,
                          datatype='SCIENCE'):
        """
    
        Convenience function which derives an exposure type given a
        combination of detector name, filter, subarray and data type
    
        :Parameters:
    
        detector: str, optional
            Name of detector (MIRIMAGE, MIRIFUSHORT or MIRIFULONG).
            If not specified, derived from the existing metadata.
        filter: str, optional
            Name of filter
            If not specified, derived from the existing metadata.
        subarray: str, optional
            Name of subarray
            If not specified, derived from the existing metadata.
        datatype: str, optional
            Known data type (SCIENCE, DARK, FLAT or TARGET)
            Defaults to SCIENCE.
            
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'exposure'):
            if hasattr(self.meta, 'instrument'):
                if detector is None:
                    detector = self.meta.instrument.detector
                if filter is None:
                    filter = self.meta.instrument.filter
            if hasattr(self.meta, 'subarray'):
                if subarray is None:
                    subarray = self.meta.subarray.name   
            if not detector:
                detector = 'UNKNOWN'
            if not filter:
                filter = 'UNKNOWN'
            if not subarray:
                subarray = 'FULL'
            exp_type = get_exp_type(detector, filter=filter, subarray=subarray,
                                    datatype=datatype)
            if exp_type:
                self.meta.exposure.type = exp_type
            else:
                # Warning removed, since it cluttered up test scripts.
                pass
#                 strg = "Insufficient information to define exposure type."
#                 warnings.warn(strg)
        else:
            strg = "Exposure metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_exposure_type(self):
        """
        
        A trivial helper function which returns the exposure type.
        
        :Parameters:
        
        None
        
        :Returns:
        
        type: str
            Exposure type. None if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'exposure'):
            if hasattr(self.meta.exposure, 'type'):
                return self.meta.exposure.type
            else:
                return None
        else:
            return None

    def set_subarray_metadata(self, subarray):
        """
        
        Convenience function to define subarray metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        subarray: str or tuple
            If a string: The name of a known subarray mode
            If a tuple: The (xstart, ystart, xsize, ysize) of a
            user-defined subarray mode, where:
            
            * xstart is the starting pixel in the X direction.
            * ystart is the starting pixel in the Y direction.
            * xsize is the number of pixels in the X direction (columns).
            * ysize is the number of pixels in the Y direction (rows).
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'subarray'):
            if isinstance(subarray, six.string_types):
                subarray = subarray.strip() # Strip off superflous white space
                self.meta.subarray.name = subarray
                subtuple = SUBARRAY_DICT[subarray]
                if subtuple is not None:
                    self.meta.subarray.xstart = subtuple[0]
                    self.meta.subarray.ystart = subtuple[1]
                    self.meta.subarray.xsize = subtuple[2]
                    self.meta.subarray.ysize = subtuple[3]   
            elif isinstance(subarray, (list,tuple)):
                self.meta.subarray.name = "GENERIC"
                self.meta.subarray.xstart = subarray[0]
                self.meta.subarray.ystart = subarray[1]
                self.meta.subarray.xsize = subarray[2]
                self.meta.subarray.ysize = subarray[3]
            else:
                raise TypeError("subarray specification must be a string or tuple")
        else:
            strg = "Subarray metadata attributes missing from data model"
            raise AttributeError(strg)


    def get_subarray_metadata(self):
        """
        
        A trivial helper function which returns the subarray metadata.
        
        :Parameters:
        
        None
        
        :Returns:
        
        (name, xstart, ystart, xsize, ysize)
            Name of subarray mode, X and Y starting pixel and X and Y size.
            (None, None, None, None, None) if not defined.
        
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'subarray'):
            return (self.meta.subarray.name,
                    self.meta.subarray.xstart, self.meta.subarray.ystart,
                    self.meta.subarray.xsize, self.meta.subarray.ysize)
        else:
            return (None, None, None, None, None)

    def set_wcs_metadata(self, wcsaxes=3, crpix=[0,0,0], crval=[0.0,0.0,0.0],
                         ctype=['','',''], cunit=['','',''],
                         cdelt=[1.0,1.0,1.0], pc=None, s_region=None,
                         waverange_start=None, waverange_end=None,
                         spectral_order=None, v2_ref=None, v3_ref=None,
                         vparity=None, v3yangle=None,
                         ra_ref=None, dec_ref=None, roll_ref=None):
        """
        
        Convenience function to define world coordinates metadata.
        Useful for setting up test data.
        
        NOTE: Data models may define different numbers of WCS axes.
        For example, the core data model defines 3 axes, but ramp data
        contains 4 axes. Warnings will be issued if a data model
        has not defined a sufficient number of axes.
        
        :Parameters:
        
        wcsaxes: int
            Maximum number of WCS axes. Defaults to 3. Maximum 4.
        crpix: list of int
            Reference pixel for each axis. Defaults to list of 0.0.
        crval: list of float
            Reference value for each axis. Defaults to list of 0.0.
        ctype: list of str
            Coordinate type for each axis. Defaults to list of ''.
        cunit: list of str
            Coordinate unit for each axis. Defaults to list of ''.
        cdelt: list of float
            Delta coordinate for each axis. Defaults to list of 1.0.
        pc: array-like, optional
            Linear transformation matrix. Defaults to None.
        s_region: string, optional
            Spatial extent of the observation. Defaults to None.
        waverange_start: float, optional
            Lower bound of default wavelength range. Defaults to None.
        waverange_end: float, optional
            Upper bound of default wavelength range. Defaults to None.
        spectral_order: int, optional
            Default spectral order. Defaults to None.
        v2_ref: float, optional
            Telescope v2 coordinate of the reference point (arcmin).
            Defaults to None.
        v3_ref: float, optional
            Telescope v3 coordinate of the reference point (arcmin).
            Defaults to None.
        vparity: int, optional
            Relative sense of rotation between Ideal xy and V2V3.
            Defaults to None.
        v3yangle: float, optional
            Angle from V3 axis to Ideal y axis (deg)
            Defaults to None.
        ra_ref: float, optional
            Right ascension of the reference point (deg). Defaults to None.
        dec_ref: float, optional
            Declination of the reference point (deg). Defaults to None.
        roll_ref, float, optional.
            Telescope roll angle of V3 measured from North over East
            at the ref. point (deg). Defaults to None.
       
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'wcsinfo'):
            # Set the number of WCS to the supplied value, to 3, or to
            # the length of the shortest of the supplied (non-empty) lists,
            # whichever is smaller.
            if wcsaxes is not None:
                naxes = wcsaxes
            else:
                naxes = 3
            if naxes > 4:
                strg = "WCSAXES=%d is larger than the maximum of 4" % naxes
                strg += " and has been truncated."
                warnings.warn(strg)
                naxes = 4
            if len(crpix) > 0 and len(crpix) < naxes:
                naxes = len(crpix)
            if len(crval) > 0 and len(crval) < naxes:
                naxes = len(crval)
            if len(ctype) > 0 and len(ctype) < naxes:
                naxes = len(ctype)
            if len(cunit) > 0 and len(cunit) < naxes:
                naxes = len(cunit)
            if len(cdelt) > 0 and len(cdelt) < naxes:
                naxes = len(cdelt)

            self.meta.wcsinfo.wcsaxes = naxes
            if naxes > 0:  
                for axis in range(0, naxes):
                    wcsmetadata = self.meta.wcsinfo
                    axis1 = axis + 1
                    crpix_name = 'crpix%d' % axis1
                    crval_name = 'crval%d' % axis1
                    ctype_name = 'ctype%d' % axis1
                    cunit_name = 'cunit%d' % axis1
                    cdelt_name = 'cdelt%d' % axis1
                    if hasattr(wcsmetadata, crpix_name):
                        setattr(wcsmetadata, crpix_name, crpix[axis])
                    else:
                        strg = "%s.%s does not exist! Too many axes?" % \
                            (wcsmetadata, crpix_name)
                        warnings.warn(strg)
                    if hasattr(wcsmetadata, crval_name):
                        setattr(wcsmetadata, crval_name, crval[axis])
                    else:
                        strg = "%s.%s does not exist! Too many axes?" % \
                            (wcsmetadata, crval_name)
                        warnings.warn(strg)
                    if hasattr(wcsmetadata, ctype_name):
                        setattr(wcsmetadata, ctype_name, ctype[axis])
                    else:
                        strg = "%s.%s does not exist! Too many axes?" % \
                            (wcsmetadata, ctype_name)
                        warnings.warn(strg)
                    if hasattr(wcsmetadata, cunit_name):
                        setattr(wcsmetadata, cunit_name, cunit[axis])
                    else:
                        strg = "%s.%s does not exist! Too many axes?" % \
                            (wcsmetadata, cunit_name)
                        warnings.warn(strg)
                    if hasattr(wcsmetadata, cdelt_name):
                        setattr(wcsmetadata, cdelt_name, cdelt[axis])
                    else:
                        strg = "%s.%s does not exist! Too many axes?" % \
                            (wcsmetadata, cdelt_name)
                        warnings.warn(strg)
            
                if pc is not None and len(pc) > 0:
                    minaxes = min( naxes, len(pc) )
                    if len(pc[0]) > 1:
                        for axis in range(0, minaxes):
                            axis1 = axis + 1
                            pc1_name = 'pc%d_1' % axis1
                            pc2_name = 'pc%d_2' % axis1
                            if hasattr(wcsmetadata, pc1_name):
                                setattr(wcsmetadata, pc1_name, pc[axis][0])
                            if hasattr(wcsmetadata, pc2_name):
                                setattr(wcsmetadata, pc2_name, pc[axis][1])

            self.set_wcs_metadata_ranges(s_region=s_region,
                                waverange_start=waverange_start,
                                waverange_end=waverange_end,
                                spectral_order=spectral_order)
            self.set_wcs_metadata_refs(v2_ref=v2_ref, v3_ref=v3_ref,
                                vparity=vparity, v3yangle=v3yangle,
                                ra_ref=ra_ref, dec_ref=dec_ref,
                                roll_ref=roll_ref)
        else:
            strg = "World coordinates metadata attributes missing from data model"
            raise AttributeError(strg)

    def set_wcs_metadata_ranges(self, s_region=None, waverange_start=None,
                                waverange_end=None, spectral_order=None):
        """
        
        Convenience function to define the optional range subset of the
        world coordinates metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        s_region: string, optional
            Spatial extent of the observation. Defaults to None.
        waverange_start: float, optional
            Lower bound of default wavelength range. Defaults to None.
        waverange_end: float, optional
            Upper bound of default wavelength range. Defaults to None.
        spectral_order: int, optional
            Default spectral order. Defaults to None.
       
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'wcsinfo'):                
            if s_region is not None:
                self.meta.wcsinfo.s_region = s_region
            if waverange_start is not None:
                self.meta.wcsinfo.waverange_start = waverange_start
            if waverange_end is not None:
                self.meta.wcsinfo.waverange_end = waverange_end
            if spectral_order is not None:
                self.meta.wcsinfo.spectral_order = spectral_order
        else:
            strg = "World coordinates metadata attributes missing from data model"
            raise AttributeError(strg)

    def set_wcs_metadata_refs(self, v2_ref=None, v3_ref=None, vparity=None,
                              v3yangle=None, ra_ref=None, dec_ref=None,
                              roll_ref=None):
        """
        
        Convenience function to define the optional reference subset of the
        world coordinates metadata.
        Useful for setting up test data.
        
        :Parameters:
        
        v2_ref: float, optional
            Telescope v2 coordinate of the reference point (arcmin).
            Defaults to None.
        v3_ref: float, optional
            Telescope v3 coordinate of the reference point (arcmin).
            Defaults to None.
        vparity: int, optional
            Relative sense of rotation between Ideal xy and V2V3.
            Defaults to None.
        v3yangle: float, optional
            Angle from V3 axis to Ideal y axis (deg)
            Defaults to None.
        ra_ref: float, optional
            Right ascension of the reference point (deg). Defaults to None.
        dec_ref: float, optional
            Declination of the reference point (deg). Defaults to None.
        roll_ref, float, optional.
            Telescope roll angle of V3 measured from North over East
            at the ref. point (deg). Defaults to None.
       
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'wcsinfo'):
            if v2_ref is not None:
                self.meta.wcsinfo.v2_ref = v2_ref
            if v3_ref is not None:
                self.meta.wcsinfo.v3_ref = v3_ref
            if vparity is not None:
                self.meta.wcsinfo.vparity = vparity
            if v3yangle is not None:
                self.meta.wcsinfo.v3yangle = v3yangle
            if ra_ref is not None:
                self.meta.wcsinfo.ra_ref = ra_ref
            if dec_ref is not None:
                self.meta.wcsinfo.dec_ref = dec_ref
            if roll_ref is not None:
                self.meta.wcsinfo.roll_ref = roll_ref       
        else:
            strg = "World coordinates metadata attributes missing from data model"
            raise AttributeError(strg)
            
    def copy_metadata(self, other, ignore=[]):
        """
        
        Copies metadata from another data model. MIRI metadata is copied
        from the other data model object to this data model.  
        Useful when copying MIRI data products.
        
        NOTE: This method will only copy the metadata defined in the
        schemas for both the current and "other" data models. Other
        metadata will be ignored (and a warning issued). Use the 'ignore'
        parameter to specify keywords that are know to be missing or
        irrelevant.
        
        :Parameters:
        
        other: MiriDataModel
            The data model whose metadata is to be copied to this one.
        ignore: list of strings, optional
            If a metadata keyword contains any of these strings it is
            not copied.
        
        """
        assert isinstance(other, MiriDataModel)
        if hasattr(other, 'meta'):
            
            # Copy all metadata apart from keyword matches specified
            # in the ignore list.
            fitskwdict = other.fits_metadata_dict()
            metakeys = list(fitskwdict.keys())
            notcopied = 0
            names = []
            for key in metakeys:
                tobecopied = True
                for ignorekw in ignore:
                    if ignorekw and (ignorekw in key):
                        # Skip to the next keyword
                        tobecopied = False
                        break
                if tobecopied:
#                     print("Copying metadata key:", key)
                    # Ignore KeyError or AttributeError exceptions
                    # when metadata is not defined in both data models.
                    try:
                        if other[key] is not None:
                            self[key] = other[key]
                    except (KeyError, AttributeError):
                        notcopied += 1
                        names.append(key)
            if notcopied > 0:
                strg = "%d metadata items could not be copied." % notcopied
                strg += "\nMissed items: "
                for name in names:
                    strg += str(name) + " "
                warnings.warn(strg)
        else:
            strg = "There is no metadata to be copied."
            warnings.warn(strg)

    #
    # Convenience functions for deducing the data arrays and data tables
    # belonging to a data model.
    # NOTE: These are expected to be called by classes derived from
    # MiriDataModel. MiriDataModel itself has no data arrays, so the
    # functions below will return an empty list.
    #
    def list_data_arrays(self):
        """
        
        Return a list of data arrays contained in the product.
        
        :Parameters:
        
        None
        
        :Returns:
        
        arraylist: tuple of str
            A list the names of the data arrays found.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def find_data_arrays(subschema, path, combiner, ctx, recurse):
            # Data objects have a datatype and do not have a type
            type = subschema.get('type')
            subdtype = subschema.get('datatype')
            if type is None and subdtype is not None:
                if not isinstance(subdtype, (tuple,list)):
                    # Ensure each entry is made only once.
                    pathstr = '.'.join(path)
                    if pathstr not in results:
#                         print("Found data array at", pathstr)
                        results.append(pathstr)

        # Walk through the schema and apply the search function.
        results = []
        mschema.walk_schema(self.schema, find_data_arrays, results)
        return results

    def list_data_tables(self):
        """
        
        Return a list of data tables contained in the product.
        
        :Parameters:
        
        None
        
        :Returns:
        
        tablelist: tuple of str
            A list the names of the data tables found.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def find_data_tables(subschema, path, combiner, ctx, recurse):
            # Data objects have a datatype and do not have a type
            type = subschema.get('type')
            subdtype = subschema.get('datatype')
            if type is None and subdtype is not None:
                if isinstance(subdtype, (tuple,list)):
                    # Skip date fields
                    if 'date' not in path:
                        # Ensure each entry is made only once.
                        pathstr = '.'.join(path)
                        if pathstr not in results:
#                             print("Found data table at", pathstr)
                            results.append(pathstr)

        # Walk through the schema and apply the search function.
        results = []
        mschema.walk_schema(self.schema, find_data_tables, results)
        return results

    def has_dataarray(self, name, matchfits=True):
        """
        
        Return True if this data product contains the named data array
        within its schema. The function can also search for a match
        against the FITS HDU names.
        
        :Parameters:
        
        name: str
            Name of data array (e.g. 'data', 'err' or 'dq').
        matchfits: bool, optional
            If True (the default) also attempt to match the given
            name to the FITS HDU names.
            
        :Returns:
        
        True or False
        
        """
        # Obtain a list of known data arrays
        list_of_arrays = self.list_data_arrays()
        
        # Check for a match in the list of data arrays
        if name in list_of_arrays:
            return True
        
        # If there isn't a match in the data array names, try looking
        # for a match in the FITS HDU name.
        if matchfits:
            for dataname in list_of_arrays:
                hduname = self.dataname_to_hduname(dataname)
                if name == hduname:
                    return True
            
        # Match not found
        return False

    def has_datatable(self, name, matchfits=True):
        """
        
        Return True if this data product contains the named data table
        within its schema. The function can also search for a match
        against the FITS HDU names.
        
        :Parameters:
        
        name: str
            Name of data table (e.g. 'data', 'err' or 'dq').
        matchfits: bool, optional
            If True (the default) also attempt to match the given
            name to the FITS HDU names.
            
        :Returns:
        
        True or False
        
        """
        # Obtain a list of known data arrays
        list_of_tables = self.list_data_tables()
        
        # Check for a match in the list of data arrays
        if name in list_of_tables:
            return True
        
        # If there isn't a match in the data array names, try looking
        # for a match in the FITS HDU name.
        for tablename in list_of_tables:
            hduname = self.dataname_to_hduname(tablename)
            if name == hduname:
                return True
            
        # Match not found
        return False

    def maskable(self):
        """
        
        Returns True if the data model contains a valid data quality
        array with at least one non-zero value, meaning the array can
        be used to mask other arrays.
        
        Note: This function doesn't check the data quality array is
        broadcastable onto other arrays, only that it exists and
        contains valid data.

        :Returns:
        
        True or False
        
        """
        result = False
        # There must be a data quality array (dq).
        if hasattr(self, 'dq'):
            # The data quality array must not be null.
            if self.dq is not None:
                # The data quality array must be an array with
                # at least one element.
                try:
                    if len(self.dq) > 0:
                        # There must be at least one non-zero value.
                        if hasattr(self.dq, 'max') and self.dq.max() > 0:
                            result = True
                except TypeError:
                    pass
        return result

    def get_field_names(self, name):
        """
        
        Get the field names associated with the named data table.

        :Parameters:
        
        name: str
            The name of the data table.
            
        :Returns:
        
        field_names: list of str
            A list of column names found in the table. An empty list
            is returned if the data element is not a table.
        
        """
        field_names = []
        dtypes = self.find_schema_components('datatype')
#         print("Looking for", name, "in table components", dtypes)
        name_dt = name + ".datatype"
        if name in dtypes:
            dtype = dtypes[name]
#             print("Found table component:", dtype, "for", name)
            if isinstance(dtype,(tuple,list)):
                for element in dtype:
#                     print("Testing element", element, "for", name)
                    if 'name' in element:
                        field_names.append(element['name'])
                return field_names
            else:
#                 warnings.warn("Data element " + name + " is not a table.")
                return []
        elif name_dt in dtypes:
            dtype = dtypes[name_dt]
#             print("Found table component:", dtype, "for", name)
            if isinstance(dtype,(tuple,list)):
                for element in dtype:
#                     print("Testing element", element, "for", name)
                    if 'name' in element:
                        field_names.append(element['name'])
                return field_names
            else:
#                 warnings.warn("Data element " + name + " is not a table.")
                return []
        else:
#             warnings.error("No data element named " + name)
            return []

    def add_comment(self, comment):
        """
        
        Add a comment to the metadata associated with the data product.
        
        :Parameters:
        
        comment: str
            A string containing a comment.
                        
        """
        # The FITS comment function no longer works.
        # Add the comment as a history item instead.
        #self.add_fits_comment(comment)
        self.add_history(comment)

    def get_comments_str(self):
        """
        
        Return a structured string containing the comments
        associated with the data product.
        
        """
        # For now, just return the FITS comments.
        return self.get_fits_comments_str()

    def add_history(self, history):
        """
        
        Add a history record to the metadata associated with the
        data product.
        
        :Parameters:
        
        history: str
            A string containing a history record.
                        
        """
        # The data model can store history records in a list-like object or
        # in a dictionary, depending on which version of the jwst library is
        # being used.
        if hasattr(self, 'history'):
            if isinstance(self.history, (tuple, list, HistoryList)):
                self.history.append( history )
            elif isinstance(self.history, dict):  
                # Wrap history string in an ASDF HistoryEntry object 
                # and tag it with the current time.
                history_item = HistoryEntry({'description': history,
                'time': Time(datetime.datetime.now())})
                if 'entries' in self.history:
                    self.history['entries'].append(history_item)
                else:
                    self.history['entries'] = [history_item]
            else:
                warnings.warn("Unrecognised self.history format! HISTORY record not added.")
        else:
            warnings.warn("Data model contains no history attribute! HISTORY record not added.")

    def add_unique_history(self, history):
        """
        
        Add a history record to the metadata associated with the
        data product, but only if the same string doesn't already
        exist.
        
        :Parameters:
        
        history: str
            A string containing a history record.
                        
        """
        history_list = self.get_history()
        if history_list:
            if history not in history_list:
                # History entry does not already exist.
                self.add_history(history)
        else:
            # This is the first history entry.
            self.add_history(history)
            
    def get_history(self):
        """
        
        Return the history list as a list of strings.
        Exists for backwards compatibility.
                
        """
        # The data model can store history records in a list-like object or
        # in a dictionary, depending on which version of the jwst library is
        # being used.
        history_list = []
        if hasattr(self, 'history'):
            if isinstance(self.history, (tuple, list, HistoryList)):
                history_list = []
                for hitem in self.history:
                    if 'description' in list(hitem.keys()):
                        history_list.append( str(hitem['description']) )
            elif isinstance(self.history, dict):
                if 'entries' in self.history:
                    hlist = list(self.history['entries'])
        return history_list
        
    def get_history_str(self):
        """
        
        Return a structured string containing the history records
        associated with the data product.
        
        """
        strg = ''
        history_list = self.get_history()
        for hitem in history_list:
            if isinstance(hitem, dict):
                if 'description' in list(hitem.keys()):
                    hstrg = str(hitem['description'])
                else:
                    hstrg = ''
                if 'time' in list(hitem.keys()):
                    tstrg = str(hitem['time'])
                else:
                    tstrg = ''
                if hstrg:
                    strg += "HISTORY = \'" + hstrg + "\'"
                if tstrg:
                    strg += "; TIME = \'" + tstrg + "\'"
            else:
                strg += "HISTORY = \'" + str(hitem) + "\'"
            strg += "\n"
        return strg

    #
    # Convenience functions for locating, getting and setting items
    # associated with FITS keywords or FITS HDUs.
    #
    # NOTE: These FITS functions will become obsolete when the underlying
    # data model stops using the FITS storage object.
    #
    def fits_metadata_dict(self, include_comments=False, include_history=False,
                           include_builtin=False):
        """
        
        Locates all the metadata items associated with a FITS keyword
        and returns a dictionary which identifies the HDU name, FITS
        keyword name and comment associated with each item.
        
        The dictionary can be used to translate each metadata tree element
        into its associated FITS HDU, keyword and comment (the opposite
        of the find_fits_keyword convenience function provided with the
        STScI data model).

        NOTE: The to_flat_dict function included in the STScI data model
        provides a similar function, but it returns a dictionary of
        metadata values only. This function is more FITS-specific.
        
        :Parameters:
        
        include_comments: bool, optional, default=False
            Set True to include COMMENT and HISTORY records in the dictionary
        include_history: bool, optional, default=False
            Set True to include COMMENT and HISTORY records in the dictionary
        include_builtin: bool, optional, default=False
            Set True to include builtin FITS keywords (such as NAXIS) in
            the dictionary.
        
        :Returns:
        
        metadict = dictionary of {schemakw: (hduname, fitskw, fitscomment)}
            A dictionary whose keywords correspond to dot-separated
            metadata tree entries, such as `meta.observation.date`.
            Each element in the dictionary is a tuple of 3 strings,
            containing the HDU name, the FITS keyword and the title
            (or FITS comment string) for the metadata element.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def find_fits_elements(subschema, path, combiner, ctx, recurse):
            hdu = subschema.get('fits_hdu')
            if not hdu:
                # If there isn't a 'fits_hdu' declaration an item
                # is assumed to be in the PRIMARY HDU.
                hdu = 'PRIMARY'
            keyw = subschema.get('fits_keyword')
            # Reject a null path or a null or blank keyword
            if not path:
                return
            if not keyw or keyw == 'BLANK':
                return
            # Only include comment, history and builtin keywords if requested.
            if not include_comments and keyw == 'COMMENT':
                return
            if not include_history and keyw == 'HISTORY':
                return
            if not include_builtin and mfits._is_builtin_fits_keyword(keyw):
                return
            # Work around bug in mfits
            if not include_builtin and (keyw == 'SIMPLE' or keyw == 'EXTEND'):
                return
            comment = subschema.get('title')
            kw = '.'.join(path)
            results[kw] = (hdu, keyw, comment)

        # Walk through the schema and apply the search function.
        results = {}
        mschema.walk_schema(self.schema, find_fits_elements, results)
        return results

    def dataname_to_hduname(self, name):
        """
        
        Return the name of the FITS HDU associated with the named data item.
        
        :Parameters:
        
        name: str
            Name of data item (e.g. 'data', 'err' or 'dq').
            
        :Returns:
        
        hduname: str
            name of the FITS HDU associated with the named item.
            Returns None if there is no match or there is no FITS HDU 
            name.
        
        """
        # TODO: Modify to use walk_schema
        if 'properties' in self.schema:
            # The schema has objects defined at the top level.
            if name in self.schema['properties']:
                if 'fits_hdu' in self.schema['properties'][name]:
                    hduname = self.schema['properties'][name]['fits_hdu']
                    return hduname
        elif 'allOf' in self.schema:
            # The schema has as "allOf" list at the top level.
            # The schema objects can be found at the next level down.
            for sublevel in self.schema['allOf']:
                if 'properties' in sublevel:
                    if name in sublevel['properties']:
                        if 'fits_hdu' in sublevel['properties'][name]:
                            hduname = sublevel['properties'][name]['fits_hdu']
                            return hduname
        else:
            raise AttributeError("MiriDataModel object has unrecognised schema format") 

        return None

    def hduname_to_dataname(self, hduname):
        """
        
        Obtain the name of the data item associated with the given
        FITS HDU.
        
        :Parameters:
        
        hduname: str
            Name of the FITS HDU (e.g. 'SCI', 'ERR' or 'DQ').
            
        :Returns:

        name: str
            Name of the data item. Returns None if there is no match.
        
        """
        # TODO: Modify to use walk_schema
        if 'properties' in self.schema:
            # The schema has objects defined at the top level.
            for subkey in list(self.schema['properties'].keys()):
                if self.schema['properties'][subkey]['type'] == 'data':
                    fits_hdu = self.schema['properties'][subkey]['fits_hdu']
                    if hduname == fits_hdu:
                        return subkey
        elif 'allOf' in self.schema:
            # The schema has as "allOf" list at the top level.
            # The schema objects can be found at the next level down.
            for sublevel in self.schema['allOf']:
                if 'properties' in sublevel:
                    for subkey in list(self.schema['properties'].keys()):
                        if self.schema['properties'][subkey]['type'] == 'data':
                            fits_hdu = \
                                self.schema['properties'][subkey]['fits_hdu']
                            if hduname == fits_hdu:
                                return subkey
        else:
            raise AttributeError("MiriDataModel object has unrecognised schema format") 

        return None

    def add_fits_comment(self, comment, hdu_name='PRIMARY'):
        """
        
        Add a comment to the comment records contained in the metadata
        associated with a given FITS HDU.
        
        NOTE: This function is deprecated. The add_comment function
        will be more reliable.
        
        :Parameters:
        
        comment: str
            A string containing a comment.
        hdu_name: str, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
            If None, any HDU is matched, but the keyword must be
            unique within the data structure.
                        
        """
        # Define a function to be applied at each level of the schema.
        def find_fits_elements(subschema, path, combiner, ctx, recurse):
            keyw = subschema.get('fits_keyword')
            hdu = subschema.get('fits_hdu')
            if not hdu:
                # If there isn't a 'fits_hdu' declaration an item
                # is assumed to be in the PRIMARY HDU.
                hdu = 'PRIMARY'
            if hdu == hdu_name and keyw == 'COMMENT':
                results.append('.'.join(path))

        # Walk through the schema and apply the search function.
        # Keep trying until the comment has been added or something
        # has gone wrong.
        count = 0
        action_pending = True
        while action_pending:
            count += 1
            if count > 2:
                # There should only be up to 2 passes through the while loop.
                strg = "\n***WARNING: Too many attempts to set "
                strg += "COMMENT keyword (%s):\n" % comment
                warnings.warn(strg)
                action_pending = False
            results = []
            mschema.walk_schema(self.schema, find_fits_elements, results)
            if len(results) > 0:
                # COMMENT keyword found. Append to the first occurence.
                clist = self[results[0]]
                clist.append(comment)
                self[results[0]] = clist
                action_pending = False
            else:
                # No existing comment keywords. A new one needs to be
                # created in the _extra_fits area. Warn if this doesn't
                # work.
                try:
                    kwdef = {hdu_name: [('COMMENT', comment, '')]}
                    # FIXME: THIS FUNCTION NO LONGER EXISTS
                    mfits.extend_schema_with_fits_keywords(self, kwdef)
                except Exception as e:
                    strg = "\n***WARNING: Could not create "
                    strg += "a first COMMENT keyword (%s):\n" % comment
                    strg += "   %s" % str(e)
                    warnings.warn(strg)
                    action_pending = False

    def get_fits_comments(self, hdu_name='ALL'):
        """
        
        Locates all the COMMENT metadata from the given HDU and
        returns them as a list.
        
        :Parameters:
        
        hdu_name: str, optional, default='ALL'
            The name of the FITS HDU from whic to extract the comments.
            If 'ALL', all comments will be included, regardless of HDU.
        
        :Returns:
        
        comments: list of list of strings
            All the comments found in the metadata. Each item in the top
            level list corresponds to one HDU and each item in the second
            level list corresponds to a separate COMMMENT entry.
                
        """
        found = self.find_fits_values('COMMENT')
        # Extract the comments and check their HDU.
        comments = []
        for item in found:
            if hdu_name == 'ALL':
                comments.append(['Comments from ' + item[0] + " HDU:"] + item[1])
            elif item[0] == hdu_name:
                comments.append(item[1])
        return comments
    
    def get_fits_comments_str(self):
        """
        
        Return a structured string containing the comments
        associated with the data product.
        
        """
        strg = ''
        commentlist = self.get_fits_comments(hdu_name='ALL')
        for comment in commentlist:
            strg += str(comment[0]) + "\n"
            for ii in range(1,len(comment)):
                strg += "\t" + str(comment[ii]) + "\n"
        return strg

    def find_fits_values(self, keyword):
        """
        
        Locates all the references to a particular FITS keyword
        in the schema and return a list of all the values found,
        identified by HDU.
        
        NOTE: Unlike 'find_fits_keyword', this function also searches
        the "_extra_fits" tree  to match keywords not included in
        the schema file. This allows the function to track down
        ad-hoc comment and history records.
        
        :Parameters:
        
        keyword: str
            The keyword to be searched for.
        
        :Returns:
        
        values: list of (string, object)
            A list of the metadata items found. Each element in the
            list is a string containing name of the HDU plus an
            object containing the value.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def find_fits_elements(subschema, path, combiner, ctx, recurse):
            keyw = subschema.get('fits_keyword')
            hdu = subschema.get('fits_hdu')
            if not hdu:
                # If there isn't a 'fits_hdu' declaration an item
                # is assumed to be in the PRIMARY HDU.
                hdu = 'PRIMARY'
            if (keyw == keyword) and path:
                kw = '.'.join(path)
                value = self[kw]
                results[kw] = (hdu, value)

        # Walk through the schema and apply the search function.
        results = {}
        mschema.walk_schema(self.schema, find_fits_elements, results)
        vlist = []
        for key in results:
            vlist.append(results[key])
        return vlist

    def get_elements_for_fits_hdu(self, schema, hdu_name='PRIMARY'):
        """
        
        Returns a list of metadata element names that are stored in a
        given FITS HDU.

        :Parameters:
        
        schema : Data model schema fragment
            Data model schema fragment

        hdu_name : str, optional
            The name of the HDU to extract.  Defaults to ``'PRIMARY'``.

        :Returns:

        elements : dict
            The keys are FITS keywords, and the values are dot-separated
            paths to the metadata elements.
            
        """
        def find_fits_elements(subschema, path, combiner, ctx, recurse):
            hdu = subschema.get('fits_hdu')
            if not hdu:
                # If there isn't a 'fits_hdu' declaration an item
                # is assumed to be in the PRIMARY HDU.
                hdu = 'PRIMARY'
            if hdu == hdu_name:
                results[subschema.get('fits_keyword')] = '.'.join(path)

        results = {}
        mschema.walk_schema(schema, find_fits_elements, results)

        return results

    def get_fits_keyword(self, keyword, hdu_name='PRIMARY'):
        """
  
        Return the value within the metadata associated with
        a FITS keyword. Only one value is returned.
     
        :Parameters:
        
        keyword: str
            The keyword to be searched for.
        hdu_name: str, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
            If None, any HDU is matched, but the keyword must be
            unique within the data structure.
        
        :Returns:
        
        value: object
            An object containing the current value associated with the keyword.
            If multiple occurences are matched, only the first will be returned.
        
        """
        if hdu_name is None:
            # Find the keyword in any HDU
            location = self.find_fits_keyword(keyword, return_result=True)
            ntimes = len(location)
            if ntimes == 1:
                return self[location[0]]
            elif ntimes < 1:
                strg = "Keyword %s not found." % keyword
                strg += " All keywords must be defined in the schema."
                raise KeyError(strg)
            else:
                strg = "\n***WARNING: Keyword %s found %d times." % (keyword, ntimes)
                strg += " Returning the first occurrence only"
                warnings.warn(strg)
                return self[location[0]]
        else:
            # Cache the search results so the same HDU isn't
            # searched for more than once.
            if hdu_name in self._metadata_cache:
                metadict = self._metadata_cache[hdu_name]
            else:
                metadict = self.get_elements_for_fits_hdu(self.schema,
                                                          hdu_name=hdu_name)
                self._metadata_cache[hdu_name] = metadict
                
            if keyword in metadict:
                return self[metadict[keyword]]
            else:
                strg = "Keyword %s not found." % keyword
                strg += " All keywords must be defined in the schema."
                raise KeyError(strg)

    def set_fits_keyword(self, keyword, value, hdu_name='PRIMARY'):
        """
  
        Set the metadata associated with a FITS keyword to a given value.
     
        :Parameters:
        
        keyword: str
            The keyword to be searched for.
        value: object
            The value to set the metadata to. The keyword must be unique.
        hdu_name: str, optional, default='PRIMARY'
            The name of the FITS HDU the keyword is associated with.
            If None, any HDU is matched.
        
        """
        if hdu_name is None:
            # Find the keyword in any HDU
            location = self.find_fits_keyword(keyword, return_result=True)
            ntimes = len(location)
            if ntimes == 1:
                self[location[0]] = value
            elif ntimes < 1:
                strg = "Keyword %s not found." % keyword
                strg += " All keywords must be defined in the schema."
                raise KeyError(strg)
            else:
                strg = "Keyword %s found %d times." % (keyword, ntimes)
                strg += " Narrow down the search by specifying the HDU name."
                raise KeyError(strg)
        else:
            # Cache the search results so the same HDU isn't
            # searched for more than once.
            if hdu_name in self._metadata_cache:
                metadict = self._metadata_cache[hdu_name]
            else:
                metadict = self.get_elements_for_fits_hdu(self.schema,
                                                          hdu_name=hdu_name)
                self._metadata_cache[hdu_name] = metadict

            if keyword in metadict:
                self[metadict[keyword]] = value
            else:
                strg = "Keyword %s not found." % keyword
                strg += " All keywords must be defined in the schema."
                raise KeyError(strg)

    #
    # Convenience functions for returning a description of a data model,
    # or a part of a data model, as as string.
    #
    def find_schema_components(self, component):
        """
        
        Return the location of a specific data array within the schema.
        
        :Parameters:
        
        component: str
            Name of the schema components to be located.
        
        :Returns:
        
        componentList: dict
            A dictionary of component values, looked up by name.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def locate_component(subschema, path, combiner, ctx, recurse):
            if component in subschema:
                value = subschema.get(component)
#                 print("Analysing component:", component, "at path", path, "which has value", value)
                kw = '.'.join(path) + "." + str(component)
                results[kw] = value

        # Walk through the schema and apply the search function.
        results = {}
        mschema.walk_schema(self.schema, locate_component, results)
        return results

    def get_title(self, underline=False, underchar="-", maxlen=None):
        """
        
        Return a top level title string describing the data product,
        which may be underlined.
        
        :Parameters:
        
        underline: bool, optional
            Set to True if the title should be underlined. The default
            is False.
        underchar: str, optional
            If underline is True, the character to be used (default '-').
        maxlen: int, optional
            The maximum title length acceptable. If not specified, the
            length is unlimited.
        
        :Returns:
        
        title: str
            A title string
        
        """
        # The schema title will either be found at the top
        # level or be inside an "allOf" structure. If the latter,
        # only the last title is returned.
        strg = self.__class__.__name__ + ": "
        if 'title' in self.schema:
            strg += self.schema['title']
        elif 'allOf' in self.schema:
            if 'title' in self.schema['allOf'][-1]:
                strg += self.schema['allOf'][-1]['title']
            else:
                strg += "product"
        # Add the subtitle defined when the data product was created, if any.
        if self._subtitle:
            strg += " (%s)" % self._subtitle

        # Truncate the title if needed.
        if maxlen is not None:
            strg = _truncate_string_left( strg, maxlen )
        if underline and strg is not None:
            len1 = len(strg)
            strg += "\n"
            strg += underchar * len1
        return strg
           
    def get_meta_title(self, underline=False, underchar="-", maxlen=None):
        """
        
        Return the metadata title.
     
        :Parameters:
        
        underline: bool, optional
            Set to True if the title should be underlined. The default
            is False.
        underchar: str, optional
            If underline is True, the character to be used (default '-').
        maxlen: int, optional
            The maximum title length acceptable. If not specified, the
            length is unlimited.
        
        :Returns:
        
        title: str
            A title string.
        
        """
        titles = self.find_schema_components('title')
        if 'meta' in titles:
            title = titles['meta']
        else:
            title = "Metadata"

        # Truncate and/or underline the title if needed.
        if maxlen is not None:
            title = _truncate_string_left( title, maxlen )
        if underline:
            len1 = len(title)
            title += "\n"
            title += underchar * len1
        return title

    def get_meta_str(self, includenulls=False, includeunits=False,
                     underline=False, underchar="-"):
        """
        
        Get a readable string describing the metadata.
     
        :Parameters:
        
        includenulls: bool, optional
            If True, the description will include all the metadata items
            with null values. The default is False.
        includeunits: bool, optional
            If True, the description will include all table units.
            The default is False.
        underline: bool, optional
            Set to True if the metadata title should be underlined.
            The default is False.
        underchar: str, optional
            If underline is True, the character to be used (default '-').
        
        :Returns:
        
        strg: str
            A string describing the metadata.
        
        """
        strg = self.get_meta_title(underline=underline, underchar=underchar) + "\n"
        # Locate all the metadata associated with FITS keywords and
        # sort the locations into order (so related keywords are grouped
        # together).
        fitskwdict = self.fits_metadata_dict()
        metakeys = list(fitskwdict.keys())
        metakeys.sort()
        # Format the metadata into a human readable table.
        strgfmt1 = "%-36s   %-8s   %-20s   %s\n"
        strgfmt2 = "%-36s | %-8s = %-20s / %s\n"
        strg += strgfmt1 % ("Data model location", "FITS key", "Value", "Comment")
        strg += strgfmt1 % ("~~~~~~~~~~~~~~~~~~~", "~~~~~~~~", "~~~~~", "~~~~~~~")
        for keyw in metakeys:
            (fitshdu, fitskey, fitscomment) = fitskwdict[keyw]
            if includeunits or ('UNIT' not in fitskey):
                value = self[keyw]
                if includenulls or (value is not None):
                    strg += strgfmt2 % (keyw, fitskey, str(value), fitscomment)
        return strg
 
    def get_data_title(self, name, underline=False, underchar="-", maxlen=None):
        """
        
        Return the title of the named data item, which may be underlined.
        
        :Parameters:
        
        name: str
            The name of the data array
        underline: bool, optional
            Set to True if the title should be underlined. The default
            is False.
        underchar: str, optional
            If underline is True, the character to be used (default '-').
        maxlen: int, optional
            The maximum title length acceptable. If not specified, the
            length is unlimited.
        
        :Returns:
        
        title: str
            A title string.
        
        """
        titles = self.find_schema_components('title')
        if name in titles:
            title = titles[name]
        else:
            title = "Data array " + name
   
        # Truncate and/or underline the title if needed.
        if maxlen is not None:
            title = _truncate_string_left( title, maxlen )
        if underline:
            len1 = len(title)
            title += "\n"
            title += underchar * len1
        return title

    def get_data_axes(self, name):
        """
        
        Return the axes of the named data item.
        
        :Parameters:
        
        name: str
            The name of the data array or table.
        
        :Returns:
        
        axes: tuple of str
            A list of axis names. Returns ['records'] if the item is a table.
            Returns None if the data name does not match, or if it has no axes.
        
        """
        titles = self.find_schema_components('axes')
        if name in titles:
            axes = titles[name]
        else:
            axes = None
            dtypes = self.find_schema_components('datatype')
            if name in dtypes:
                dtype = dtypes[name]
                if isinstance(dtype,(tuple,list)):
                    axes = ["records"]
        return axes

    def get_data_labels(self, name):
        """
        
        Return a human readable string labelling the axes of the
        named data item.
        
        :Parameters:
        
        name: str
            The name of the data array.
        
        :Returns:
        
        label: str
            A string describing the structure of the data axes.
        
        """
        axes = self.get_data_axes(name)
        labels = None
        if hasattr(self, name):
            data = getattr(self, name)
            if data is not None and hasattr(data,'shape'):
                labels = _shape_to_string(data.shape, axes)
        return labels

    def get_data_units(self, name):
        """
        
        Return the units of the named data item (data array or data table)
        
        :Parameters:
        
        name: str
            The name of the data array or data table
        
        :Returns:
        
        units: str
            The data item units. Returns None if the data name does not 
            match, or if it has no units.
            Returns a single string for a data array.
            Returns a string containing a list of column units for a data table.
        
        """
        # The data units should be defined in the metadata.
        units = None
        if hasattr( self.meta, name):
            metaname = getattr(self.meta, name)
            if hasattr(metaname, 'units'):
                # A data array with data units defined in meta.units
                units = getattr(metaname, 'units')
            elif hasattr(metaname, 'tunit1'):
                # A data table with column units defined in
                # meta.tunit1, meta.tunit2, etc...
                col = 1
                unitname = 'tunit%s' % col
                colunits = []
                while (hasattr(metaname, unitname)):
                    colunits.append('\"' + str(getattr(metaname, unitname)) + '\"')
                    col += 1
                    unitname = 'tunit%s' % col
                # Return a string containing the column units formatted
                # as a comma-separated list of strings.
                units = "(" +  ", ".join(colunits) + ")"
                # Convert undefined units to null strings.
                units = units.replace('\"None\"', '\"\"')
        
        # If no units are defined in the metadata, see if there are
        # default units in the schema.
        if not units:
            unitlist = self.find_schema_components('units')
            if name in unitlist:
                units = unitlist[name]
        return units

    def set_data_units(self, name, units=None):
        """
        
        Define the units of the named data array.
        
        :Parameters:
        
        name: str
            The name of the data array
        units: str (optional)
            The required data units to be defined.
            If not specified (or None) the data units will be set
            to the default value defined in the schema. If there
            is no default, no unit will be defined.
            A warning will be issued if the data model already
            defines units that are different from those expected
            by the schema.
        
        :Returns:
        
        units: str
            The actual data units as defined.
            
        :Raises:
        
        AttributeError:
            Raised if there are no data units defined in the schema
            or if the data item does not exist.
        
        """
        # If the units have not been explicitly specified, get a default
        # from the schema. Search for a .units components of the named item.
        defaulted = False
        if units is None:
            unitlist = self.find_schema_components('units')
            for keyword in unitlist:
                if name in keyword:
                    units = str(unitlist[str(keyword)])
                    # Convert undefined units to a null string.
                    if units == 'None':
                        units = ''
                    defaulted = True
                    break
#                 print("+++Default of %s found in schema" % units)

        # If there are valid units, define them by setting the metadata
        # (if available).
        if units is not None:      
            if hasattr( self.meta, name):
                metaname = getattr(self.meta, name)
                if hasattr(metaname, 'units'):
                    metaunits = str(getattr(metaname, 'units'))
                    # Convert undefined units to a null string.
                    if metaunits == 'None':
                        metaunits = ''
                    if defaulted and metaunits:
                        # The metadata already defines the units.
                        # Leave the units as defined, but give a warning if
                        # there is a discrepancy with the default units
                        # defined in data model schema.
                        if metaunits.strip() != units.strip():
                            strg = "\n***Data units for \'%s\' " % name
                            strg += " defined in the metadata "
                            strg += "(%s) do not match the default units " % \
                                metaunits
                            strg += "defined in the schema (%s)." % \
                                units
                            warnings.warn(strg)
                    else:
                        # Data units have been explicitly provided or have
                        # only been found in the data model schema.
                        setattr(metaname, 'units', units)
                else:
                    # Cannot define units because the schema doesn't
                    # define the appropriate metadata.
                    strg = "Cannot define units. Schema does not define "
                    strg += "a self.meta.%s.units attribute." % name
                    raise AttributeError(strg)
            else:
                # Cannot define units because the schema doesn't
                # define the appropriate metadata.
                strg = "Cannot define data units. Schema does not define "
                strg += "a self.meta.%s attribute." % name
                raise AttributeError(strg)
        return units

#     def set_table_units(self, name, units=None):
#         """
#          
#         Define the column units of the named data table.
#          
#         :Parameters:
#          
#         name: str
#             The name of the data table
#         units: list of str (optional)
#             The required table column units to be defined.
#             If not specified (or None) the column units will be set
#             to the default value defined in the schema. If there
#             is no default, no unit will be defined.
#             A warning will be issued if the data model already
#             defines units that are different from those expected
#             by the schema.
#          
#         :Returns:
#          
#         units: list of str
#             The actual table column units as defined.
#              
#         :Raises:
#          
#         AttributeError:
#             Raised if there are no table units defined in the schema
#             or if the data item does not exist.
#          
#         """
#         if units is not None and len(units) > 0:
# #             print("+++Table column units for %s of %s provided explicitly" % (name,str(units)))
#             # Table column units have been explicitly provided.
#             # Set each column of the metadata to the defined units.
#             if hasattr( self.meta, name):
#                 metaname = getattr(self.meta, name)
#             for col in range(0,len(units)):
#                 # Convert null units into null strings
#                 if units[col] is None:
#                     units[col] = ''
#                 unitname = 'tunit%s' % (col + 1)
#                 if metaname and hasattr(metaname, unitname):
#                     metaunits = getattr(metaname, unitname)
# #                     if metaunits:
# #                         print("+++Changing %s.%s from %s to %s" % \
# #                             (name, unitname, str(metaunits), str(units[col])))
# #                     else:
# #                         print("+++Setting %s.%s to %s." % \
# #                             (name, unitname, str(units[col])))
#                     setattr(metaname, unitname, str(units[col]))
#                 else:
#                     strg = "\n***Table column unit %s.%s is not defined in the schema metadata" % (name,unitname)
#                     warnings.warn(strg)
#         else:
#             # No units provided explicitly.
#             # Leave the metadata as it is (but check that something
#             # is defined).
#             if hasattr( self.meta, name):
#                 metaname = getattr(self.meta, name)
#                 MAXCOLS = 24
#                 defined = False
#                 for col in range(0,MAXCOLS):
#                     unitname = 'tunit%s' % (col + 1)
#                     if metaname and hasattr(metaname, unitname):
#                         defined = True
# #                         metaunits = getattr(metaname, unitname)
# #                         print("+++Leaving table column unit %s.%s set to %s." % \
# #                                 (name, unitname, str(metaunits)))
#                     else:
#                         break
#                 if not defined:
#                     strg = "\n***No table units for table \'%s\' defined in metadata or in schema" % name
#                     warnings.warn(strg)
#         return units
    
    def get_data_stats(self, name):
        """
        
        Get a string summarising the contents of the named data item.
        The following statistics are included:
        
        * Number of zero elements
        * Number of non-zero elements
        * Minimum value
        * Maximum value
        * Mean value
        * Median value (if appropriate)
        * Standard devision (if appropriate)
        
        NOTE: If a masked version of a data array is available,
        the masked version is used to derive the statistics.
        
        :Parameters:
        
        name: str
            The name of the data item.
        
        :Returns:
        
        description: str
            A string giving statistics for the named data item.
        
        """
        strg = ''
        data_title = self.get_data_title(name)
        strg += data_title
        data_units = self.get_data_units(name)
        if data_units is not None and data_units:
            strg += " in \'%s\'" % data_units
        data_labels = self.get_data_labels(name)
        if data_labels is not None and data_labels:
            strg += " (%s)" % data_labels
        if hasattr(self, name):
            # NOTE: The "_masked" variation of a data array is a feature
            # of the HasDataErrAndDq class in operations.py.
            if self.maskable() and hasattr(self, name + "_masked"):
                data = getattr(self, name + "_masked")
                # Fill a masked array to prevent a numpy warning.
                data = data.filled()
                strg += " - masked array"
            else:
                data = getattr(self, name)
            if data is not None and hasattr(data, 'size') and (data.size > 0):
                nnzero = np.count_nonzero(data)
                nzero = data.size - nnzero
                strg += "\n  zero=%d, non-zero=%d; " % (nzero, nnzero)
                strg += "min=%f, max=%f, mean=%f" % \
                    (np.min(data), np.max(data), np.mean(data))
                if data.size > 2:
                    strg += ", median=%f" % \
                        (np.median(data))
                if data.size > 3:
                    strg += ", std=%f" % \
                        (np.std(data))
                strg += "\n"
            else:
                strg += "\n  EMPTY DATA STRUCTURE\n"
        else:
            strg += "\n  NO DATA STRUCTURE FOUND\n"
        return strg

    def get_data_str(self, name, underline=False, underchar="-"):
        """
        
        Get a string describing the named data item.
        
        :Parameters:
        
        name: str
            The name of the data item.
        underline: bool, optional
            Set to True if the title should be underlined. The default
            is False.
        underchar: str, optional
            If underline is True, the underline character to be used
            (default '-').
         
        :Returns:
        
        description: str
            A string describing the named data item.
        
        """
        strg = '\n'
        data_title = self.get_data_title(name)
        strg += data_title
        data_units = self.get_data_units(name)
        if data_units is not None and data_units:
            strg += " in \'%s\'" % data_units
        data_labels = self.get_data_labels(name)
        if data_labels is not None and data_labels:
            strg += " (%s)" % data_labels
        len1 = len(strg)
        fieldnames = self.get_field_names(name)
        strg2 = ""
        if fieldnames:
            strg2 += "\nColumns: \'"
            strg2 += "\', \'".join(fieldnames)
            strg2 += "\'"
            strg += strg2
        len2 = len(strg2)
        if underline:
            len3 = max(len1, len2)
            strg += "\n"
            strg += underchar * len3
        if hasattr(self, name):
            data = getattr(self, name)
            strg += "\n%s\n" % str(data)
        else:
            strg += "\nNO DATA STRUCTURE FOUND\n"
        return strg

    def stats(self):
        """
        
        Return statistics about a generic MIRI data object as a
        readable string.
        
        The string contains a list of data arrays together with
        statistics for the range of values contained in those
        data arrays.
        
        :Returns:
        
        description: str
            A string containing statistics of the contents.
        
        """
        # Start with the data object title
        strg = self.get_title(underline=True, underchar="=") + "\n"

        # Display the data arrays.
        list_of_arrays = self.list_data_arrays()
        narrays = len(list_of_arrays)
        for dataname in list_of_arrays:
            strg += self.get_data_stats(dataname)
        
        # Count the data tables.
        list_of_tables = self.list_data_tables()
        ntables = len(list_of_tables)
        if ntables > 0:
            if ntables == 1:
                ss = ''
            else:
                ss = 's'
            if narrays > 0:
                strg += "Plus "
            strg += "%d data table%.1s: \n  " % (ntables, ss)
            strg += ', '.join(str(s) for s in list_of_tables)
            if ntables == 1:
                strg += "\n"
                strg += self.get_data_str(list_of_tables[0], underline=True,
                                          underchar="-")
        return strg

    def get_title_and_metadata(self):
        """
        
        Return the title and metadata of a MIRI data object as a readable
        string.
        
        :Returns:
        
        description: str
            A string containing a description.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        strg += self.get_history_str()
        return strg

    def __str__(self):
        """
        
        Return the contents of a generic MIRI data object as a readable
        string.
        
        :Returns:
        
        description: str
            A string containing a summary of the contents.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()
        
        # Display the data arrays.
        list_of_arrays = self.list_data_arrays()
        for dataname in list_of_arrays:
            strg += self.get_data_str(dataname, underline=True, underchar="-")
        
        # Display the data tables.
        list_of_tables = self.list_data_tables()
        for tablename in list_of_tables:
            strg += self.get_data_str(tablename, underline=True, underchar="-")
        return strg

    def plot(self, description='', visitor=None, **kwargs):
        """
        
        Plot the data product using the algorithm contained in the
        specified visitor class.

        NOTE: This plot method can be considered a quick look method.
        
        :Parameters:
        
        description: str, optional
            Description to be added to the plot title, if required.
            This is only used if a visitor object is not provided.
        visitor: object, optional
            A visitor object which implements the algorithm for plotting the
            contents of the data structure. If not specified, a temporary
            object will be created from the matplotlib-based visitor class
            and used.
        \*\*kwargs:
            All other keyword arguments will be passed to the visitor's visit
            method (for data arrays). These can be keywords recognised by
            the plotting package.

        """
        # If a visitor object is not given, create one using the default class.
        if visitor is None:
            visitor = DataModelPlotVisitor(description=description)
            
        # Invoke the visitor plot method on this data object.
        visitor.visit(self, **kwargs)

    def dump_schema(self):
        """
        
        Dump the schema to a heirarchical string
        
        :Returns:
        
        strg: str
            String version of the schema layout.
        
        :Requires:
        
        Makes use of the search facilities in jwst.datamodels.schema
        
        """
        # Define a function to be applied at each level of the schema.
        def schema_to_string(subschema, path, combiner, ctx, recurse):
            keys = list(subschema.keys())
            #level = len(path.split('.'))
            for key in keys:
                value = subschema.get(key)
                if isinstance(value, dict):
                    results.append('.'.join(path) + ' | ' + key + \
                                   "= vv (expanded below) vv")
                else:
                    results.append('.'.join(path) + ' | ' + key + "= " + \
                                   str(value))

        # Walk through the schema and apply the search function.
        results = []
        mschema.walk_schema(self.schema, schema_to_string, results)
        
        strg = self.__class__.__name__ + " schema dump:\n"
        for line in results:
            strg += line + "\n"
        return strg

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriDataModel module.")

    PLOTTING = False
    SAVE_FILES = False

    # Create an empty MIRI data model and add some metadata.
    print("Empty MIRI data model")
    with MiriDataModel() as testdata:
        # Dump the schema to a string
        dump_strg = testdata.dump_schema()
        # Check the string is valid. It can be printed out, but generates
        # a lot of output.
        assert(dump_strg is not None)
        assert(len(dump_strg) > 0)
#         print("\n" + dump_strg + "\n")
        # Define some metadata.
        testdata.set_housekeeping_metadata('UKATC', 'Joe Bloggs', 'GROUND', 'V1.0')
        #testdata.set_observation_metadata() # Use current data/time
        testdata.set_observation_metadata(obsid='001',
                                 obsnumber='001', programnumber='1235',
                                 visitid='001', visitnumber='001', visitgroup='01',
                                 activityid='01', exposurenumber='0001')
        testdata.set_pointing_metadata(ra_v1=12.20, dec_v1=-7.15, pa_v3=2.0)
        testdata.set_target_metadata(ra=12.0, dec=-7.0)
        testdata.set_instrument_metadata(detector='MIRIMAGE', filt='F560W',
                                channel='', ccc_pos='OPEN', 
                                deck_temperature=10.0,
                                detector_temperature=7.0)
        testdata.set_wcs_metadata(wcsaxes=2, waverange_start=5.0,
                                  waverange_end=25.0, ra_ref=12.10,
                                  dec_ref=-7.03)
        testdata.set_exposure_metadata(readpatt='SLOW', nints=1, ngroups=3)
        testdata.set_exposure_type()
        testdata.set_subarray_metadata('FULL')
        testdata.set_subarray_metadata('GENERIC')
        testdata.set_subarray_metadata( (1,1,16,16) )
        testdata.set_subarray_metadata('BRIGHTSKY')
        testdata.add_comment( \
            "This object has been created to try out the MiriDataModel class.")
        testdata.add_history( \
            "This line has been added to the history as a test.")
        
        # Try reading back the metadata
        (detector, modelnam, detsetng) = testdata.get_instrument_detector()
        print("Read back detector metadata:", detector, modelnam, detsetng)
        filt = testdata.get_instrument_filter()
        (channel, band) = testdata.get_instrument_channel_band()
        print("Read back instrument metadata:", filt, channel, band)
        (readpatt, nints, ngroups, nframes) = testdata.get_exposure_readout()
        print("Read back exposure metadata:", readpatt, nints, ngroups, nframes)
        (name, xstart, ystart, xsize, ysize) = testdata.get_subarray_metadata()
        print("Read back subarray metadata:", name, xstart, ystart, xsize, ysize)
        
        # Display the metadata as a formatted string.
        # Also tests the fits_metadata_dict() function.
        print("Metadata string:\n" + testdata.get_meta_str())
        
        # Search for specific items in the metadata.
        print("All occurences of BUNIT: " + \
              str(testdata.find_fits_values('BUNIT')))
        print("All occurences of DATE: " + \
              str(testdata.find_fits_values('DATE')))
        print("All comments:\n" + testdata.get_comments_str())
        print("All history:\n" + testdata.get_history_str())
        
        print("Get INSTRUME FITS keyword: " + \
              testdata.get_fits_keyword('INSTRUME'))
        print("Get TELESCOP FITS keyword: " + \
              testdata.get_fits_keyword('TELESCOP'))
        print("Get SUBARRAY FITS keyword: " + \
              testdata.get_fits_keyword('SUBARRAY'))
        
        # The following will be empty lists
        dataarrays = testdata.list_data_arrays()
        assert len(dataarrays) == 0
        datatables = testdata.list_data_tables()
        assert len(datatables) == 0
        
        print(testdata)
        if PLOTTING:
            testdata.plot(description="Empty data model")
        if SAVE_FILES:
            testdata.save("test_base_model.fits", overwrite=True)

        del testdata
        
    print("Test finished.")
