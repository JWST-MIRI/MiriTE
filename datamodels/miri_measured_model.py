#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which adds the ability to 
mask the DATA and ERR arrays based on the contents of the DQ array. It 
also adds the ability to add, subtract, multiply and divide the data 
objects. Used for MIRI measured data.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

20 Nov 2012: Created
19 Dec 2012: Extended to cover image and ramp variants of the measured
             data model. Typo corrected in _generate_masks function.
             MIRI schema names prefixed with "miri.".
             Implemented scalar operations, with placeholders for other
             arithmetic operator functions.
11 Jan 2013: Array masking reworked and arithmetical operators almost
             completed.
23 Jan 2013: Added plotting.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
28 Feb 2013: Extended the class with the STScI world coordinates
             methods, HasFitsWcs.
13 May 2013: Added MiriSlopeModel to describe MIRI slope data
             (which is different from "ImageModel" data because it
             preserves integrations). N.B. FINAL MODEL IS TBD.
17 May 2013: Create a default table object explicitly to stop
             the underlying data model filling the table with
             junk by default.
21 May 2013: Added rough plot_ramp function to MiriRampModel.
04 Jun 2013: Shortened the names of the ramp, slope and image models.
02 Jul 2013: MiriCubeModel added. Don't display masked data when it
             is the same as unmasked data.
05 Aug 2013: Removed unnecessary setting of field_def to None.
14 Aug 2013: Modified MiriRampModel to replace the DQ array with PIXELDQ
             and GROUPDQ arrays (as done in jwst.datamodels). A maskwith parameter
             is used to determine which array is used to mask the data.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
04 Oct 2013: Data quality array changed to 32-bit. Changed default field_def
             table to use MIRI reserved flags.
07 Oct 2013: GROUP_DEF table added to MIRI ramp data. 'both' option added
             to maskwith parameter for ramp data. The exception raised
             when attempting to set the DQ array in ramp data replaced by
             a warning, so arithmetic operations still work. (This needs
             further work.)
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
             Added some even valued data quality flags to the module tests.
23 Jan 2014: Modified for jsonschema draft 4: Functions made more
             independent of schema structure. Additional data arrays added
             to MiriSlopeData string conversion.
11 Apr 2014: Modified to define data units using the set_data_units method.
20 Jun 2014: Issue a warning if ramp data could not be stored in a uint16
             format file.
25 Jun 2014: Added dq_def and groupdq_def alongside field_def and group_def.
             field_def and group_def are deprecated.
             Automatically create DQ_n header keywords from the dq_def,
             pixeldq_def and groupdq_def binary tables. 
23 Jul 2014: DQ definitions removed from metadata.
25 Sep 2014: reserved_flags replaced by master_flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             DQ_DEF table changed from 3 columns to 4 columns. Do not attempt
             to convert a FIELD_DEF table into a DQ_DEF table automatically.
29 Sep 2014: Corrected a mistake in the logic of the creation of default
             DQ_DEF, PIXELDQ_DEF and GROUPDQ_DEF tables.
30 Sep 2014: Added the convert_flags method, to aid CDP-3 delivery.
09 Dec 2014: Obsolete field_def and group_def tables removed.
11 Jun 2015: Obsolete methods commented out at CDP-3 delivery removed.
09 Jul 2015: Reference output array (refout) added to MiriRampModel schema.
18 Aug 2015: Obsolete dq_metadata schema and associated test code removed.
             Duplicated dq_def schema also removed.
             Removed MiriImageModel and MiriCubeModel.
20 Aug 2015: Data array metadata moved into a separate schema.
08 Sep 2015: Ramp and slope data given a more descriptive data type.
11 Sep 2015: Removed duplicated plot method.
25 Jan 2016: Changes for compatibility with ASDF version of jwst.datamodels:
             Removed reference to HasFitsWcs.
28 Jan 2016: Modified to meet the new requirements of the ASDF release of
             the jwst.datamodels package. Made more robust against missing
             attributes. MiriSimpleModel added.
05 Apr 2016: Changed ramp model schema from miri_ramp to miri_ramp_withmeta
             to avoid namespace clash with jwst.datamodels.
04 May 2016: Zero frame array (zeroframe) added to MiriRampModel.
             ERR array removed from ramp data model.
15 Jul 2016: Reworded the warning message issued when ramp data are
             outside the range to be stored in a uint16 data array.
06 Sep 2016: Operations work again. Test code restored.
07 Feb 2017: MiriRampModel made more resilient against missing data quality
             attributes.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
12 Sep 2017: Added 4-axis WCS metadata to the MIRI ramp data model.
04 Oct 2017: Define the meta.filetype metadata.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import warnings, logging
import numpy as np
import numpy.ma as ma

# Import the MIRI reserved data quality flags and flags table class
from miri.datamodels.dqflags import master_flags, \
    pixeldq_flags, groupdq_flags, FlagsTable, insert_value_column, convert_dq

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.operations import HasData, HasDataErrAndDq

# List all classes and global functions here.
__all__ = ['MiriSimpleModel', 'MiriMeasuredModel', 'MiriRampModel',
           'MiriSlopeModel']


class MiriSimpleModel(MiriDataModel, HasData):
    """
    
    A simple data model for MIRI, with a data array only, based on
    the STScI base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the data parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    data: numpy array (optional)
        An array containing the science data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_simple.schema.yaml"
    
    def __init__(self, init=None, data=None, **kwargs):
        """
        
        Initialises the MiriSimpleModel class.
        
        Parameters: See class doc string.

        """
        super(MiriSimpleModel, self).__init__(init=init, **kwargs)
        
        # Update the data array if it has been explicitly provided.
        # Otherwise use the array defined by the init or by the default
        # data object.
        HasData.__init__(self, data)

        # Copy the units of the data array from the schema, if defined.
        dataunits = self.set_data_units('data')

    def __str__(self):
        """
        
        Return the contents of the simple data object as a readable
        string. The function will display the contents of the data
        array.
        
        :Parameters:
        
        None
            
        :Returns:
        
        strg: str
            A string displaying the contents of the simple data object.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        
        # Display the data array.
        strg += self.get_data_str('data', underline=True, underchar="-")
        return strg

# EXTRA FUNCTIONS TEMPORARILY COPIED FROM HasData abstract class
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


class MiriMeasuredModel(MiriDataModel, HasDataErrAndDq):
    """
    
    A data model for MIRI data with error handling and masking, based on
    the STScI base model, DataModel.

    After a data model has been created, the data, err and dq arrays are
    available as attributes .data, .err and .dq. Metadata items are available
    within a ".meta" attribute tree.
    
    See http://ssb.stsci.edu/doc/jwst_dev/jwst.datamodels.doc/html/index.html.
        
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the data parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    data: numpy array (optional)
        An array containing the science data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    err: numpy array (optional)
        An array containing the error data.
        Must be broadcastable onto the data array.
    dq: numpy array (optional)
        An array containing the quality data.
        Must be broadcastable onto the data array.
    dq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'bad','Bad data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_measured.schema.yaml"
    dq_def_names = ('BIT', 'VALUE', 'NAME', 'TITLE')

    # Set the default dq_def table to the JWST master flags
    # TODO: Can this default be defined in the schema?
    _default_dq_def = master_flags
    
    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 rampdata=False, **kwargs):
        """
        
        Initialises the MiriMeasuredModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMeasuredModel, self).__init__(init=init, **kwargs)
        
        # Update the data arrays if they have been explicitly provided.
        # Otherwise use the arrays defined by the init or by the default
        # data object.
        self.rampdata = rampdata
        if self.rampdata:
            HasDataErrAndDq.__init__(self, data, None, dq, noerr=True)
        else:
            HasDataErrAndDq.__init__(self, data, err, dq, noerr=False)

        # Copy the units of the data array from the schema, if defined.
        # The error array should have the same units as the data array,
        # unless explicitly defined differently.
        dataunits = self.set_data_units('data')
        if not self.rampdata:
            errunits = self.set_data_units('err')
            if not errunits:
                self.set_data_units('err', dataunits)

        # The ramp data model has its own way of defining data quality.
        if not self.rampdata:
            # Set the data quality bit field definitions table, if provided.
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
                # TODO: Can the default declared in the schema be used?
                self.dq_def = self._default_dq_def

    def on_save(self, path):
        """
        
        Override the on_save method of DataModel to ensure the metadata
        and DQ flags table are in step before saving the model.

        :Parameters:
        
        path : str
            The path to the file that we're about to save to.

        """
#         if hasattr( self.meta, "dq"):
#             flags_table_to_metadata( self.flags_table, self.meta.dq )
        super(MiriMeasuredModel, self).on_save(path)

    def __str__(self, extra_objects=True):
        """
        
        Return the contents of the measured data object as a readable
        string. The string will include the contents of the data, err
        and dq arrays and the dq_def table. It will also include any
        extra data arrays found in the data structure, unless the
        extra_objects flags is set False.
        
        :Parameters:
        
        extra_objects: boolean (optional)
            A hidden extra parameter, which can be set to False to turn
            off the display of extra data objects.
            
        :Returns:
        
        strg: str
            A string displaying the contents of the measured data object.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        
        # Display the data arrays, including their masked aliases
        # (unless the data quality array is full of good values).
        strg += self.get_data_str('data', underline=True, underchar="-")
        if hasattr(self, 'data_masked'):
            data_masked = self.data_masked
            if self.maskable() and (data_masked is not None):
                title = self.get_data_title('data') + " (masked)"
                len2 = len(title)
                title += "\n" + (len2 * '~')
                strg += title + "\n"  + str(data_masked) + "\n"
        else:
            warnings.warn("%s object does not have a \'data_masked\' attribute!" % \
                          self.__class__.__name__)

        if not self.rampdata:
            strg += self.get_data_str('err', underline=True, underchar="-")
            if hasattr(self, 'err_masked'):
                err_masked = self.err_masked
                if self.maskable() and (err_masked is not None):
                    title = self.get_data_title('err') + " (masked)"
                    len2 = len(title)
                    title += "\n" + (len2 * '~')
                    strg += title + "\n" + str(err_masked) + "\n"
            else:
                warnings.warn("%s object does not have an \'err_masked\' attribute!" % \
                          self.__class__.__name__)
            
        strg += self.get_data_str('dq', underline=True, underchar="-")
        strg += self.get_data_str('dq_def', underline=True, underchar="-")
        
        if extra_objects:
            extras = self.extra_objects_str()
            if extras:
                strg += "Extra Data Objects\n=-=-=-=-=-=-=-=-=-\n"
                strg += extras
        return strg
    
    def extra_objects_str(self):
        """
        
        Return the contents of extra objects (in addition to those
        expected in a measured model) as a readable string.
        
        """
        strg = ''
        list_of_arrays = self.list_data_arrays()
        if 'data' in list_of_arrays:
            list_of_arrays.remove('data')
        if 'err' in list_of_arrays:
            list_of_arrays.remove('err')
        if 'dq' in list_of_arrays:
            list_of_arrays.remove('dq')
        if list_of_arrays:
            for dataname in list_of_arrays:
                strg += self.get_data_str(dataname, underline=True,
                                          underchar="-")

        list_of_tables = self.list_data_tables()
        if 'dq_def' in list_of_tables:
            list_of_tables.remove('dq_def')
        if list_of_tables:
            for tablename in list_of_tables:
                strg += self.get_data_str(tablename, underline=True,
                                          underchar="-")
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


class MiriRampModel(MiriMeasuredModel):
    """
    
    A data model for MIRI ramp data with error handling and masking, like
    MiriMeasuredModel, but with additional restrictions to ensure the
    underlying data model is compatible with the STScI RampModel.
    
    NOTE: Raw, uncalibrated JWST data comes in 3 forms:
    
    * level 0: Compressed science telemetry data.
    
    * level 1a: Uncompressed original FITS files.
    
    * level 1b: Uncalibrated FITS files, with additional metadata.
    
    See JWST-STScI-002111, SM-12, "DMS Level 1 and 2 Data Product Design
    Technical Report", Revision - A, 20 Dec 2012.
    
    This data model is designed to represent level 1b data, which is
    stored in floating point format. The original science telemetry
    data is communicated in unsigned 16-bit integer format (uint16).
    A warning will be issued if the values stored in the level 1b
    data arrays could not have originated from uint16 data. 
   
    :Parameters:
    
    The same as MiriMeasuredModel, except without the err array and
    with the addition of
    
    refout: numpy array (optional)
        A 4-D array containing reference output information.
    zeroframe: numpy array (optional)
        A 3-D array containing zero frame information.
    
    The dq and dq_def parameters are replaced by:
    
    pixeldq: numpy array (optional)
        A 2-D array containing the pixel quality data (for all groups).
    groupdq: numpy array (optional)
        A 4-D array containing the specific quality data for particular
        groups.
    maskwith: str (optional)
        A string declaring which of the two data quality arrays
        ('pixeldq', 'groupdq' or 'both') should be used to mask the data arrays
        during arithmetic operations. The default is 'both', or whichever
        of the two arrays is not None if only one array is available.
    pixeldq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the pixeldq array. For
        example: [(0, 'bad','Bad data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    groupdq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the groupdq array. For
        example: [(0, 'bad','Bad group'), (1, 'saturated', 'Saturated group')]
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the default MIRI ramp group
        flags.
    
    """
    # Floating point data is assumed because the STScI RampModel also
    # declares floating point data.
    schema_url = "miri_ramp_withmeta.schema.yaml"

    _default_dq_def = pixeldq_flags
    _default_groupdq_def = groupdq_flags

    def __init__(self, init=None, data=None, pixeldq=None, groupdq=None,
                 maskwith=None, err=None, refout=None, zeroframe=None,
                 dq_def=None, pixeldq_def=None, groupdq_def=None, **kwargs):
        super(MiriRampModel, self).__init__(init=init, data=data, dq=None,
                                            dq_def=None, rampdata=True,
                                            **kwargs)
        # Data type is MIRI ramp data (level 1b).
        self.meta.model_type = 'Ramp (level 1b)'
        self.meta.filetype = 'Ramp (level 1b)'
        
        # Define the REFOUT and ZEROFRAME data arrays which are unique to
        # ramp data.
        if refout is not None:
            self.refout = refout
        if zeroframe is not None:
            self.zeroframe = refout

        # Define the pixeldq and groupdq arrays if they are explicitly given.
        # Define the maskwith parameter if explcitly given, or if not given,
        # to 'both' or to whichever of the two dq arrays are valid.
        if pixeldq is not None:
            self.pixeldq = pixeldq
        if groupdq is not None:
            self.groupdq = groupdq
        if maskwith is not None:
            self.maskwith = maskwith
        elif self._isvalid(self.pixeldq) and self._isvalid(groupdq):
            self.maskwith = 'both'
        elif self._isvalid(self.pixeldq):
            self.maskwith = 'pixeldq'
        elif self._isvalid(groupdq):
            self.maskwith = 'groupdq'
        else:
            self.maskwith = 'none'

        # If 4-D data is given, the 3rd dimension gives the number of groups
        # and the 4th dimension gives the number of integrations.
        if data is not None and hasattr(data, "shape"):
            if len(data.shape) == 4 :
                nints = data.shape[0]
                ngroups = data.shape[1]
                if self.meta.exposure.nints is None:
                    self.meta.exposure.nints = nints
                if self.meta.exposure.ngroups is None:
                    self.meta.exposure.ngroups = ngroups

        # Set the pixel data quality bit field definitions table, if provided.
        if pixeldq_def is not None:
            try:
                self.pixeldq_def = pixeldq_def
            except (ValueError, TypeError) as e:
                strg = "pixeldq_def must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        elif self.pixeldq_def is None or len(self.pixeldq_def) < 1:
            # No pixeldq_def is provided.
            # Explicitly create a PIXELDQ_DEF table with default values.
            # TODO: Can the default declared in the schema be used?
            self.pixeldq_def = self._default_dq_def

#         # Initialise the pixel DQ metadata keywords
#         if init is not None and self.flags_table is not None:
#             flags_table_to_metadata( self.flags_table, self.meta.pixeldq)

        # Set the group bit field definitions table, if provided.
        if groupdq_def is not None:
            try:
                self.groupdq_def = groupdq_def
            except (ValueError, TypeError) as e:
                strg = "groupdq_def must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        elif self.groupdq_def is None or len(self.groupdq_def) < 1:
            # No groupdq_def is provided.
            # Explicitly create a GROUPDQ_DEF table with default values.
            # TODO: Can the default declared in the schema be used?
            self.groupdq_def = self._default_groupdq_def

#         # Initialise the GROUPDQ metadata keywords
#         if init is not None and self.groupdq_flags_table is not None:
#             flags_table_to_metadata( self.groupdq_flags_table, self.meta.groupdq)

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
            NOT IMPLEMENTED.
        show_ints: boolean, optional, default=False
            Set to True to distinguish the ramps belonging to different
            integrations. Only works when a single row and column is
            specified.
        description: string, optional
            Additional description to be shown on the plot, if required. 
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # TODO: TO BE COMPLETED PROPERLY USING VISITOR CLASS
        
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt

        # Generate a data cube view of the ramp data, where integrations
        # and groups are combined.
        if len(self.data.shape) == 4 :
            nints = self.data.shape[0]
            ngroups = self.data.shape[1]
            npoints = nints*ngroups
            datacube = self.data.reshape(npoints, self.data.shape[2], self.data.shape[3])
        else:
            datacube = self.data
            
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
            if description:
                tstrg += "\n" + description

            if show_ints:
                # Plot integrations separately.
                # Create a matplotlib figure and axis.
                fig = mplt.new_figure(1, stitle=tstrg)
                ax = mplt.add_subplot(fig, 1, 1, 1)

                linefmtlist = ['bo', 'rx', 'g+', 'ko', \
                               'bx', 'r+', 'go', 'kx', \
                               'b+', 'ro', 'gx', 'k+' ]
                itmin = 0.0
                for integ in range(0, nints):
                    itmax = itmin + ngroups * stime
                    itime = np.linspace(itmin, itmax, ngroups)
                    itmin = itmax + stime

                    # Extract the ramp for this integration
                    # from the data cube.
                    ramp = self.data[integ,:,rows,columns]

                    # The plots will cycle around the defined line formats,
                    ifmt = integ % len(linefmtlist)
                    lfmt = linefmtlist[ifmt] 

                    # Plot each ramp as an XY plot overlaid on the same axis.
                    mplt.plot_xy(itime, ramp, plotfig=fig, plotaxis=ax,
                                 linefmt=lfmt, xlabel=xlabel,
                                 ylabel=ylabel, title='')
            else:
                # Plot all integrations together.
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
            if description:
                tstrg += "\n" + description
            
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
        mplt.show_plot()

    def on_save(self, path):
        """
        
        Override the on_save method of MiriMeasuredModel to issue a warning
        if the data are outside the range to have originated from uint16
        telemetry data.

        :Parameters:
        
        path : str
            The path to the file that we're about to save to.

        """
        if self.data is not None and len(self.data) > 0:
            _MAXUINT16 = 65535
            data_max = self.data.max()
            if data_max > _MAXUINT16:
                strg = "\nNOTE: The maximum value stored in the data array (%g) " % data_max
                strg += "is outside the range of unsigned 16-bit telemetry data."
                warnings.warn(strg)

        super(MiriRampModel, self).on_save(path)

    def __str__(self):
        """
        
        Return the contents of the ramp data object as a readable
        string.
        
        """
        # Start with the data object title and metadata
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        
        # Display the data arrays, including their masked aliases
        # (unless the data quality array is full of good values).
        strg += self.get_data_str('data', underline=True, underchar="-")
        data_masked = self.data_masked
        if self.maskwith == 'both':
            maskedwith = 'both pixeldq and groupdq'
        else:
            maskedwith = self.maskwith
        if self.maskable() and (data_masked is not None):
            title = self.get_data_title('data') + " (masked with %s)" % \
                maskedwith
            len2 = len(title)
            title += "\n" + (len2 * '~')
            strg += title + "\n"  + str(data_masked) + "\n"

        strg += self.get_data_str('refout', underline=True, underchar="-")
        strg += self.get_data_str('zerodata', underline=True, underchar="-")
           
        if hasattr(self, 'pixeldq'):
            strg += self.get_data_str('pixeldq', underline=True, underchar="-")
        if hasattr(self, 'groupdq'):
            strg += self.get_data_str('groupdq', underline=True, underchar="-")
        if hasattr(self, 'pixeldq_def'):
            strg += self.get_data_str('pixeldq_def', underline=True, underchar="-")
        if hasattr(self, 'groupdq_def'):
            strg += self.get_data_str('groupdq_def', underline=True, underchar="-")
        return strg

    # "flags_table" is a FlagsTable object created on the fly
    # from the contents of the pixeldq_def table
    @property
    def flags_table(self):
        if hasattr(self, 'pixeldq_def') and self.pixeldq_def is not None:
            # Convert the pixeldq_def table into a FlagsTable object
            # and return it.
            return FlagsTable( self.pixeldq_def )
        else:
            return None

    @flags_table.setter
    def flags_table(self, data):
        raise AttributeError("The flags_table object is read-only")

    # "groupdq_flags_table" is a FlagsTable object created on the fly
    # from the contents of the groupdq_def table
    @property
    def groupdq_flags_table(self):
        if hasattr(self, 'groupdq_def') and self.groupdq_def is not None:
            # Convert the groupdq_def table into a FlagsTable object
            # and return it.
            return FlagsTable( self.groupdq_def )
        else:
            return None

    @groupdq_flags_table.setter
    def groupdq_flags_table(self, data):
        raise AttributeError("The groupdq_flags_table object is read-only")

    @property
    def dq(self):
        # Alias dq for pixeldq, groupdq or both combined as specified.
        if self.maskwith == 'groupdq':
            if hasattr(self, 'groupdq'):
                return self.groupdq
            else:
                return None
        elif self.maskwith == 'pixeldq':
            if hasattr(self, 'pixeldq'):
                return self.pixeldq
            else:
                return None
        else:
            # Combine both sets of flags together, if they both exist.
            if hasattr(self, 'groupdq') and self._isvalid(self.groupdq) and \
               hasattr(self, 'pixeldq') and self._isvalid(self.pixeldq):
                return self.groupdq | self.pixeldq
            elif hasattr(self, 'groupdq') and self._isvalid(self.groupdq):
                return self.groupdq
            elif hasattr(self, 'pixeldq') and self._isvalid(self.pixeldq):
                return self.pixeldq
            else:
                return None

    # TODO: The following function needs further work.
    # It's used during arithmetic operations.
    @dq.setter
    def dq(self, data):
        # Alias dq for pixeldq, groupdq or both combined as specified.
        if self.maskwith == 'groupdq':
            logging.info("Mask results written to GROUPDQ array.")
            self.groupdq = data
        elif self.maskwith == 'pixeldq':
            logging.info("Mask results written to PIXELDQ array.")
            self.pixeldq = data
        else:
            # One set of data can't be used to update both the PIXELDQ
            # and GROUPDQ arrays.
            strg = "\n***DQ array is a combined view of PIXEL DQ and GROUPDQ. "
            strg += "Mask results written to GROUPDQ array only."
            warnings.warn(strg)
            # TODO: Create the pixeldq by shrinking the given array?
            self.groupdq = data


class MiriSlopeModel(MiriMeasuredModel):
    """
    
    A data model for MIRI slope data with error handling and masking, like
    MiriMeasuredModel, but with additional restrictions to ensure the
    underlying data model is compatible with the STScI RampModel.
    
    :Parameters:
    
    The same as MiriMeasuredModel, plus
    
    zeropt: numpy array
        An array containing the zero point of the fit.
        Must be broadcastable onto the data array.
    nreads: numpy array
        An array containing the number of good frames used in the fits
        Must be broadcastable onto the data array.
    readsat: numpy array
        An array containing the frame number of the first saturated frame.
        Must be broadcastable onto the data array.
    ngoodseg: numpy array
        An array containing the number of good segments used in the slope fit.
        Must be broadcastable onto the data array.
    fiterr: numpy array
        An array containing the RMS error in the slope fit.
        Must be broadcastable onto the data array.
    meandata: numpy array (optional)
        An array containing the mean slope data.
        Must be broadcastable onto the data array.
    meanerr: numpy array (optional)
        An array containing the uncertainty in the mean slope data.
        Must be broadcastable onto the data array.
    meandq: numpy array (optional)
        An array containing the quality of the mean slope data.
        Must be broadcastable onto the data array.
    
    """
    schema_url = "miri_slope.schema.yaml"

    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 zeropt=None, nreads=None, readsat=None, ngoodseg=None,
                 fiterr=None, fitinfo=None, **kwargs):
        super(MiriSlopeModel, self).__init__(init=init, data=data, dq=dq,
                                             err=err, dq_def=dq_def,
                                             **kwargs)
        # Data type is MIRI slope data.
        self.meta.model_type = 'Slope (level 2)'
        self.meta.filetype = 'Slope (level 2)'
        
        # If 3-D data is given, the 3rd dimension gives the number of
        # integrations.
        if data is not None and hasattr(data, "shape"):
            if len(data.shape) == 3:
                nints = data.shape[0]
                if self.meta.exposure.nints is None:
                    self.meta.exposure.nints = nints

        # NOTE: The additional metadata and housekeeping arrays
        # added below are TO BE DECIDED. Some of the extra data is
        # not going to be needed by the STScI pipeline.
                    
        # If fit information is provided, use it to update the metadata.
        if fitinfo is not None and isinstance(fitinfo, dict):
            if 'NPINT' in fitinfo:
                self.meta.fit.npint = fitinfo['NPINT']
            if 'NSFITS' in fitinfo:
                self.meta.fit.nsfits = fitinfo['NSFITS']
            if 'NSFITE' in fitinfo:
                self.meta.fit.nsfite = fitinfo['NSFITE']
            if 'HIGHSAT' in fitinfo:
                self.meta.fit.highsat = fitinfo['HIGHSAT']

        if nreads is not None:
            self.nreads = nreads
#         self._nreads_mask = None
#         self._nreads_fill = 0
#         self._nreads_fill_value = None

        if readsat is not None:
            self.readsat = readsat
#         self._readsat_mask = None
#         self._readsat_fill = -1
#         self._readsat_fill_value = None

        if ngoodseg is not None:
            self.ngoodseg = ngoodseg
#         self._ngoodseg_mask = None
#         self._ngoodseg_fill = 0
#         self._ngoodseg_fill_value = None
        if zeropt is not None:
            self.zeropt = zeropt
        self._zeropt_mask = None
        self._zeropt_fill = 'min'
        self._zeropt_fill_value = None

        if fiterr is not None:
            self.fiterr = fiterr
        self._fiterr_mask = None
        self._fiterr_fill = 'max'
        self._fiterr_fill_value = None

        # Copy the units of the these arrays, if defined.
        self.set_data_units('zeropt')
        self.set_data_units('fiterr')
        
    def __str__(self):
        """
        
        Return the contents of the slope object as a readable
        string.
        
        """
        # First obtain a string describing the underlying measured
        # model.
        strg = super(MiriSlopeModel, self).__str__(extra_objects=False)
        
        # Add the extras
        maxdq = self.dq.max()

        strg += self.get_data_str('nreads', underline=True, underchar="-")
        strg += self.get_data_str('readsat', underline=True, underchar="-")
        strg += self.get_data_str('ngoodseg', underline=True, underchar="-")

        if self.zeropt is not None:
            strg += self.get_data_str('zeropt', underline=True, underchar="-")
            if maxdq > 0:
                if self.maskable():
                    zeropt_masked = self.zeropt_masked
                    if zeropt_masked is not None:
                        title = self.get_data_title('zeropt') + " (masked)"
                        len2 = len(title)
                        title += "\n" + (len2 * '~')
                        strg += title + "\n" + str(zeropt_masked) + "\n"
        if self.fiterr is not None:
            strg += self.get_data_str('fiterr', underline=True, underchar="-")
            if maxdq > 0:
                if self.maskable():
                    fiterr_masked = self.fiterr_masked
                    if fiterr_masked is not None:
                        title = self.get_data_title('fiterr') + " (masked)"
                        len2 = len(title)
                        title += "\n" + (len2 * '~')
                        strg += title + "\n" + str(fiterr_masked) + "\n"
        return strg

    @property
    def zeropt_masked(self):
        # Generate the masked data on the fly. This ensures the
        # masking is always up to date with the latest dq array.
        # TODO: Can this result be cached and the cache invalidated
        # when either the zeropt or dq arrays change?
        if self.zeropt is not None and self.dq is not None:
            self._zeropt_mask = self._generate_mask(self.zeropt, self.dq)
            self._zeropt_fill_value = self._generate_fill(self.zeropt,
                                                        self._zeropt_fill)
            return ma.array(self.zeropt, mask=self._zeropt_mask,
                            fill_value=self._zeropt_fill_value)
        else:
            return self.zeropt

    @property
    def zeropt_filled(self):
        masked = self.zeropt_masked
        if masked is not None:
            return masked.filled(self._zeropt_fill_value)
        else:
            return self.zeropt

    @property
    def fiterr_masked(self):
        # Generate the masked data on the fly. This ensures the
        # masking is always up to date with the latest dq array.
        # TODO: Can this result be cached and the cache invalidated
        # when either the fiterr or dq arrays change?
        if self.fiterr is not None and self.dq is not None:
            self._fiterr_mask = self._generate_mask(self.fiterr, self.dq)
            self._fiterr_fill_value = self._generate_fill(self.fiterr,
                                                          self._fiterr_fill)
            return ma.array(self.fiterr, mask=self._fiterr_mask,
                            fill_value=self._fiterr_fill_value)
        else:
            return self.fiterr

    @property
    def fiterr_filled(self):
        masked = self.fiterr_masked
        if masked is not None:
            return masked.filled(self._fiterr_fill_value)
        else:
            return self.fiterr


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriMeasuredModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,2,1],[0,1,0]])
    dqa3x3 = np.array([[0,2,0],[1,0,1],[0,0,0]])
    dqb3x3 = np.array([[0,1,0],[0,1,0],[2,1,0]])

    dqdef = master_flags
    grpdef = groupdq_flags
    
    tenarray = np.ones_like(data3x3) * 10.0

    print("Completely null measured data")
    with MiriMeasuredModel(title='null data') as nulldata:
        print(nulldata)
        print(nulldata.stats())
 
    print("Scalar measured data")
    with MiriMeasuredModel( data=42, title='scalar data' ) as scalardata:
        print(scalardata)
        print(scalardata.stats())

    print("Measured data with data + err + dq:")
    with MiriMeasuredModel(data=data3x3, err=err3x3, dq=dq3x3, dq_def=dqdef,
                           title='data+err+dq') \
            as testdata:
        print("Data arrays: ", testdata.list_data_arrays())
        print("Data tables: ", testdata.list_data_tables())
        print(testdata)
        print(testdata.stats())
        if PLOTTING:
            testdata.plot("Test data 1 - with data + err + dq")
        if SAVE_FILES:
            testdata.save("test_measured_model.fits", overwrite=True)
 
        print("Add 1 to testdata")
        newdata = testdata + 1
        print(newdata)
        del newdata
 
        print("Subtract 1 from testdata")
        newdata = testdata - 1
        print(newdata)
        del newdata
 
        print("Multiply testdata by 10.")
        newdata = testdata * 10
        print(newdata)
        del newdata
 
        print("Divide testdata by 10.")
        newdata = testdata / 10
        print(newdata)
        del newdata
     
        print("Add an array full of 10s to testdata.")
        newdata = testdata + tenarray
        print(newdata)
        del newdata
 
        print("Subtract an array full of 10s from testdata.")
        newdata = testdata - tenarray
        print(newdata)
        del newdata
 
        print("Multiply testdata by an array full of 10s.")
        newdata = testdata * tenarray
        print(newdata)
        del newdata
 
        print("Divide testdata by an array full of 10s.")
        newdata = testdata / tenarray
        print(newdata)
        del newdata

        print("Image data with data + err + dq:")
        with MiriMeasuredModel(data=data3x3, err=err3x3, dq=dq3x3,
                                    dq_def=dqdef) \
                as testdata2:
            print(testdata2)
            print(testdata2.stats())
            if PLOTTING:
                testdata2.plot("Test data 2 - image with data + err + dq")
            if SAVE_FILES:
                testdata2.save("test__image_model.fits",
                                  overwrite=True)
        
            print("Add two image objects together")
            newdata = testdata + testdata2
            print(newdata)
            del newdata
 
            print("Subtract one image object from another")
            newdata = testdata - testdata2
            print(newdata)
            del newdata
 
            print("Multiply two image objects together")
            newdata = testdata * testdata2
            print(newdata)
            del newdata
 
            print("Divide one image object by another")
            newdata = testdata / testdata2
            print(newdata)
            del newdata, testdata2
        del testdata

    print("\nRamp data with data + refout + pixeldq:")
    ramp3x3x2x2 = np.array([[data3x3,1.2*data3x3],[1.4*data3x3,1.6*data3x3]])
    refout = np.ones_like(ramp3x3x2x2)
    with MiriRampModel(data=ramp3x3x2x2, refout=refout,
                       pixeldq=dq3x3, pixeldq_def=dqdef) \
            as testdata3:
        print(testdata3)
        print(testdata3.stats())
        if PLOTTING:
            testdata3.plot("Test data 3 - ramp with data + refout + pixeldq")
            testdata3.plot_ramp(1,1)
        if SAVE_FILES:
            testdata3.save("test__ramp_model.fits", overwrite=True)
        del testdata3

    print("\nRamp data with data + refout + groupdq:")
    ramp3x3x2x2 = np.array([[data3x3,1.2*data3x3],[1.4*data3x3,1.6*data3x3]])
    dqa3x3x2 = np.array([dqa3x3,dqb3x3])
    dqb3x3x2 = np.array([dqb3x3,dqa3x3])
    dq3x3x2x2 = np.array([dqa3x3x2,dqa3x3x2])
    with MiriRampModel(data=ramp3x3x2x2, refout=refout,
                       pixeldq=None, groupdq=dq3x3x2x2, groupdq_def=grpdef) \
            as testdata3a:
        print(testdata3a)
        print(testdata3a.stats())
        if PLOTTING:
            testdata3a.plot("Test data 3a - ramp with data + err + refout + groupdq")
            testdata3a.plot_ramp(1,1)
        if SAVE_FILES:
            testdata3a.save("test__ramp_modela.fits", overwrite=True)
        del testdata3a

    print("\nRamp data with data + refout + pixeldq + groupdq:")
    with MiriRampModel(data=ramp3x3x2x2, refout=refout, pixeldq=dq3x3,
                       groupdq=dq3x3x2x2, pixeldq_def=dqdef, groupdq_def=grpdef) \
            as testdata3b:
        print(testdata3b)
        print(testdata3b.stats())
        if PLOTTING:
            testdata3b.plot("Test data 3b - ramp with data + refout " \
                            + "+ pixeldq + groupdq")
            testdata3b.plot_ramp(1,1)
        if SAVE_FILES:
            testdata3b.save("test__ramp_modelb.fits", overwrite=True)
        del testdata3b

    print("\nRamp data with data only")
    with MiriRampModel( data=ramp3x3x2x2 ) \
            as testdata3c:
        print(testdata3c)
        print(testdata3c.stats())
        if PLOTTING:
            testdata3c.plot("Test data 3c - ramp with data only")
            testdata3c.plot_ramp(1,1)
        if SAVE_FILES:
            testdata3c.save("test__ramp_modelc.fits", overwrite=True)
        del testdata3c

    print("\nEmpty ramp data")
    with MiriRampModel( ) \
            as testdata3d:
        print(testdata3d)
        if SAVE_FILES:
            testdata3d.save("test__ramp_modeld.fits", overwrite=True)
        del testdata3d

    print("\nSlope data with data + err + dq only:")
    slope3x3x4 = np.array([data3x3,data3x3,data3x3,data3x3])
    err3x3x4 = np.array([err3x3,err3x3,err3x3,err3x3])
    dq3x3x4 = np.array([dq3x3,dq3x3,dq3x3,dq3x3])
    with MiriSlopeModel( data=slope3x3x4, err=err3x3x4, dq=dq3x3x4,
                         dq_def=dqdef) \
            as testdata4:
        print(testdata4)
        print(testdata4.stats())
        if PLOTTING:
            testdata4.plot("Test data 4 - slope with data + err + dq only")
        if SAVE_FILES:
            testdata4.save("test__slope_model1.fits", overwrite=True)
        del testdata4

    print("\nSlope data with data + err + dq + zeropt + nreads, etc...:")
    zeropt3x3 = np.array([[0.,4.,3.],[2.,0.,3.5],[3.2,0.5,0.7]])
    zeropt3x3x4 = np.array([zeropt3x3,zeropt3x3+1,zeropt3x3,zeropt3x3*2])
    nreads3x3 = np.array([[12,11,10],[10,12,9],[10,9,11]])
    nreads3x3x4 = np.array([nreads3x3,nreads3x3,nreads3x3,nreads3x3])
    ngoodseg3x3x4 = np.ones_like(nreads3x3x4)
    fiterr3x3x4 = np.array([err3x3,err3x3,err3x3,err3x3])
    fitinfo = {'NPINT':4, 'NSFITS':2, 'NSFITE':5, 'HIGHSAT':5800}
    with MiriSlopeModel(data=slope3x3x4, err=err3x3x4, dq=dq3x3x4,
                               dq_def=dqdef, zeropt=zeropt3x3x4,
                               nreads=nreads3x3x4, ngoodseg=ngoodseg3x3x4,
                               fiterr=fiterr3x3x4, fitinfo=fitinfo) \
            as testdata5:
        print(testdata5)
        print(testdata5.stats())
        if PLOTTING:
            testdata5.plot("Test data 5 - slope with data + err + dq + ...")
        if SAVE_FILES:
            testdata5.save("test__slope_model2.fits", overwrite=True)
        del testdata5

    print("\nData cube with data + err + dq.")
    cube3x3x4 = np.array([data3x3,data3x3,data3x3,data3x3])
    err3x3x4 = np.array([err3x3,err3x3,err3x3,err3x3])
    dq3x3x4 = np.array([dq3x3,dq3x3,dq3x3,dq3x3])
    with MiriMeasuredModel(data=cube3x3x4, err=err3x3x4, dq=dq3x3x4,
                               dq_def=dqdef) \
            as testdata6:
        # Test the setting of some metadata.
        testdata6.set_instrument_metadata(detector='MIRIMAGE', filt='F560W',
                                          ccc_pos='OPEN',
                                          deck_temperature=10.0,
                                          detector_temperature=7.0)
        print(testdata6)
        print(testdata6.stats())
        print("Data arrays=", testdata6.list_data_arrays())
        print("Data tables=", testdata6.list_data_tables())
        print("FlagsTable:\n" + str(testdata6.flags_table))       
        if PLOTTING:
            testdata6.plot("Test data 6 - data cube with data + err + dq")
        if SAVE_FILES:
            testdata6.save("test__slope_model1.fits", overwrite=True)
        del testdata6
        
    print("Test finished.")
