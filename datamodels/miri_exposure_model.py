#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the MIRI data model for "level 1b" exposure data. Similar to the
MIRI ramp data model, except the data may be built up a piece
at a time instead of being created all at once. Groups and/or integrations
can also be averaged together before being written to disk.

This module is intended to be used by a MIRI simulator for managing
exposure data.

:Reference:

The STScI jwst.datamodels documentation:

https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

and the documentation on JWST file formats:

https://jwst-docs.stsci.edu/display/JDAT/JWST+File+Names%2C+Formats%2C+and+Data+Structures

:History:

06 Feb 2013: Created (based on the ExposureDataProduct in the old data model).
07 Mar 2013: Corrected some typos.
04 Jun 2013: Shortened the names of the ramp, slope and image models.
10 Jun 2013: Resolved problem with rows and columns being declared in the
             wrong order. Use non-square test data. Remodelled in parallel
             with changes to the SCASim exposure data model.
13 Jun 2013: Modified to use the miri.exposure.schema.yaml schema instead
             of miri.ramp.schema.yaml
14 Aug 2013: Modified for new MIRI ramp model.
07 Oct 2013: GROUP_DEF table added to MIRI ramp data.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
09 Jun 2014: Added ability to generate slope data (inefficiently).
20 Jun 2014: Added set_simulation_metadata, get_group and get_integration
             methods.
02 Mar 2015: SUBBURST keyword added to simulator metadata schema.
04 Mar 2015: Corrected typo in metadata schema.
06 Aug 2015: Added set_exposure_times method. Added mean amplifier property
             keywords to schema metadata. The groupdq array is now intialised
             to zero and updated as the exposure data arrays are constructed.
31 Aug 2015: Exposure start, middle and end time are now floating point
             values rather than date-time strings.
03 Sep 2015: Rearranged the data model schema to put simulator metadata
             into a separate file. Removed unwanted print statements in
             cube_data method.
22 Sep 2015: get_exposure_times method added.
07 Oct 2015: Made exception catching Python 3 compatible.
27 Oct 2015: Raise an exception if expected metadata attributes are missing.
30 Nov 2015: Tightened up the data typing so that all type conversions inherit
             the data types declared in the constructor.
26 Jan 2016: Date objects converted to string before writing to metadata.
04 Apr 2016: Ensure no class attributes are defined before calling the
             superclass constructor. Convert parameters to integer just once.
06 Apr 2016: Undid the renaming of some attributes.
03 May 2016: Added the REFOUT extension needed by the STScI pipeline.
             Added refrows and refcolumns to the constructor.
             Added _extract_refout method, which is called from the save method.
             Added default ZEROFRAM keyword to metadata.
04 May 2016: Added get_exp_type function and set_exposure_type method.
             Made the PIXELDQ and GROUPDQ arrays optional.
             ERR array removed from ramp data model.
09 May 2016: Moved get_exposure_type and get_exp_type to miri_model_base.py.
08 Jun 2016: Corrected a typo in the definition of __all__.
03 Aug 2016: Included a check on the integrity of the size of the reference
             output data.
20 Oct 2016: Explicitly truncate reference rows to an integer before applying
             a numpy slice.
07 Feb 2017: pixeldq and groupdq arrays removed from exposure data when the
             include_pixeldq and include_groupdq flags are not enabled. The
             pixeldq_def and groupdq_def data tables will be removed as well.
17 Mar 2017: Corrected some documentation typos.
23 Jun 2017: Added the ability to add DARK data in situ (for simulation).
07 Jul 2017: Check the size of the dark before attempting to add it.
             Skip the reference rows when adding a DARK.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
18 Jul 2017: More work on the add_dark function. Clip the output at 1 DN.
20 Jul 2017: Also clip the output above 65535 DN in add_dark.
27 Jul 2017: Added cosmic ray count to metadata.
11 Sep 2017: Ensure EXPSTART, EXPMID and EXPEND keywords are MJD.
13 Dec 2017: Corrected the default frame time for the slope_data function.
19 Feb 2018: Added the ability to make nonlinearity correction in situ
             (for simulation).
27 Feb 2018: Much more efficient nonlinearity correction.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.exposure_model") # Get a default parent logger

# import warnings
import numpy as np
#import numpy.ma as ma
import scipy.stats
import copy

# Import the MIRI ramp data model utilities.
from miri.datamodels.miri_measured_model import MiriRampModel
from miri.datamodels.plotting import DataModelPlotVisitor

# List all classes and global functions here.
__all__ = ['MiriExposureModel']

#
# Global helper function
#
def linear_regression(x, y):
    """
    
    Linear regression function. Fit an intercept and a slope
    to a set of x,y coordinates
    
    """
    #print("Fitting", y, "\nagainst", x)
    matrix = np.vstack( [x, np.ones_like(x)] ).T
    slope, intercept = np.linalg.lstsq(matrix,y)[0]
    #print("gives slope=", slope, "intercept=", intercept)
    return (slope, intercept)


class MiriExposureModel(MiriRampModel):
    """
    
    Class MiriExposureModel - science exposure data product.

    A variation of MiriRampModel used to create science exposure data.
    Unlike the a MiriRampModel object, which is created in one go, the
    MiriExposureModel object is initialized full of zeros and can be created
    all at once or built up one group at a time.
    This class is intended to be used by MIRI data simulators.
    
    :Parameters:
    
    rows: int
        The number of rows making up the data.
    columns: int
        The number of columns making up the data.
    ngroups: int
        The number of readout groups making up each integration.
    nints: int
        The number of integrations making up each exposure.
    readpatt: str
        Name of detector readout pattern.
    refrows: int
        The number of reference output rows included with the data.
        Set to zero if there is no reference output.
    refcolumns: int
        The number of reference output columns included with the data.
        Set to zero if there is no reference output.
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
    pixeldq: numpy array (optional)
        A 2-D array containing the pixel quality data (for all groups).
        Typically the bad pixel mask which applies to an exposure.
        NOTE: The groupdq array is expected to be constructed at the
        same time as the exposure data.
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
        If not specified, it will default to the JWST PIXELDQ flags.
    groupdq_def: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (value:int, name:str, title:str),
        giving the meaning of values stored in the groupdq array. For
        example: [(0, 'bad','Bad group'), (1, 'saturated', 'Saturated group')]
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the default JWST GROUPDQ flags
        flags.
    include_pixeldq: bool, optional
        Set to False to remove the pixeldq data array and pixeldq_def table
        from the exposure data. The default is True.
    include_groupdq: bool, optional
        Set to False to remove the groupdq data array and the groupdq_def table
        from the exposure data. The default is True.
    
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    # Exposure data has, by definition, the same structure as ramp data.
    schema_url = "miri_exposure.schema.yaml"

    def __init__(self, rows, columns, ngroups, nints, readpatt, refrows,
                 refcolumns, grpavg=1, intavg=1, nframes=1,
                 groupgap=0, pixeldq=None, maskwith=None, pixeldq_def=None,
                 groupdq_def=None, include_pixeldq=True, include_groupdq=True,
                 **kwargs):
        """
        
        Initialises the MiriExposureModel class.
        
        Parameters: See class doc string.

        """
        # Check the integrity of data array dimensions.
        int_rows = int(rows)
        int_columns = int(columns)
        if int_rows <= 0 or int_columns <= 0:
            strg = "The exposure data must have a non-zero size "
            strg += "(rows=%d, columns=%d given)." % (rows, columns)
            raise ValueError(strg)

        int_ngroups = int(ngroups)
        if int_ngroups <= 0:
            strg = "There must be at least one group "
            strg += "(%d given)." % ngroups
            raise ValueError(strg)
        int_nframes = int(nframes)
        int_groupgap = int(groupgap)
        int_grpavg = int(grpavg)
        
        int_nints = int(nints)
        if int_nints <= 0:
            strg = "There must be at least one integration "
            strg += "(%d given)." % nints
            raise ValueError(strg)
        int_intavg = int(intavg)
        
        # Check the integrity of the requested data averaging. It is
        # assumed the caller will provide nicely rounded values.
        if int_grpavg < 1 or (int_ngroups % int_grpavg != 0):
            strg = "%d groups does not divide exactly by %d. " % \
                (int_ngroups, int_grpavg)
            strg += "Please adjust the number of groups."
            raise ValueError(strg)    
        if int_intavg < 1 or (int_nints % int_intavg != 0):
            strg = "%d integrations does not divide exactly by %d. " % \
                (int_nints, int_intavg)
            strg += "Please adjust the number of integrations."
            raise ValueError(strg)    
        
        if (pixeldq is not None) and (not include_pixeldq):
            strg = "Incompatible arguments. A pixeldq array is provided "
            strg += "when include_pixeldq=False. The array is ignored."
            LOGGER.error(strg)
            pixeldq = None
            
        if not include_pixeldq:
            pixeldq_def = None

        # Check the integrity of the reference output dimensions (if given)
        int_refrows = int(refrows)
        int_refcolumns = int(refcolumns)
        if int_refrows > 0 and int_refcolumns > 0:
            refarea = int_refrows * int_refcolumns
            if (refarea % int_columns) > 0:
                strg = "Reference output data area (%d rows x %d columns) " % \
                    (int_refrows, int_refcolumns)
                strg += "will not divide by the data columns (%d)." % int_columns
                raise ValueError(strg)
       
        # Create zero-filled arrays of the correct shape.
        # Floating point is used here because the STScI RampModel also
        # declares floating point data.
        datashape = [int_nints, int_ngroups, int_rows, int_columns]
        zerodata = np.zeros(datashape, dtype=np.float32)
        if include_groupdq:
            zerogdq = np.zeros(datashape, dtype=np.uint8)
        else:
            zerogdq = None
            groupdq_def = None
        
        if int_refrows > 0 and int_refcolumns > 0:
            refoutshape = [int_nints, int_ngroups, int_refrows, int_refcolumns]
            zerorefout = np.zeros(refoutshape, dtype=np.float32)
        else:
            # Create an empty reference output with a single value (otherwise
            # it ends up the same dimensions as the main data array).
            zerorefout = np.zeros([1,1,1,1], dtype=np.float32)
            #zerorefout = None
        
        super(MiriExposureModel, self).__init__(data=zerodata, err=None,
                                                refout=zerorefout,
                                                pixeldq=pixeldq,
                                                groupdq=zerogdq,
                                                maskwith=maskwith,
                                                pixeldq_def=pixeldq_def,
                                                groupdq_def=groupdq_def,
                                                title="MIRI exposure data",
                                                **kwargs)

        # Initially, no reference output frame has been extracted.
        self._refout_extracted = False
        self.refrows = int_refrows
        self.refcolumns = int_refcolumns
        # By default, there is no zero frame.
        self.meta.exposure.zero_frame = False
        self.include_pixeldq = include_pixeldq
        self.include_groupdq = include_groupdq

        # Add exposure parameter attributes
        self.rows = int_rows
        self.columns = int_columns
        self.ngroups = int_ngroups
        self.grpavg = int_grpavg
        self.nints = int_nints
        self.intavg = int_intavg
        self.nframes = int_nframes
        self.groupgap = int_groupgap
        self.readpatt = readpatt
        self.ngroups_file = self.ngroups // self.grpavg
        self.nints_file = self.nints // self.intavg

        # Start with no averaged data
        self._data_averaged = None

        # Ensure the metadata has the correct exposure and subarray
        # parameters.
        self.set_exposure_metadata(readpatt=readpatt, nints=int_nints,
                                   ngroups=int_ngroups, nframes=int_nframes,
                                   grpavg=int_grpavg, intavg=int_intavg,
                                   groupgap=int_groupgap)
        subtuple = (1, 1, int_columns, int_rows)
        self.set_subarray_metadata(subtuple)

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

    def set_simulation_metadata(self, crmode, burst_mode=None):
        """
        
        Convenience function to define simulation metadata.
        Used when the exposure data is created by a simulator.
        
        :Parameters:
        
        crmode: str
            Cosmic ray mode.
        burst_mode: boolean, optional
            Set True if subarrays are read out in burst mode.
            The parameter is only written to the metadata if
            specified.
            
        """
        if hasattr(self, 'meta') and hasattr(self.meta, 'simulator'):
            self.meta.simulator.cosmic_ray_mode = crmode
            if burst_mode is not None:
                self.meta.subarray_burst_mode = burst_mode
        else:
            strg = "***Simulation metadata attributes missing from data model"
            raise AttributeError(strg)

    def set_exposure_times(self, exposure_time=None, duration=None,
                           start_time=None, mid_time=None, end_time=None):
        """
        
        Convenience function to define exposure times in the metadata.
        If called with no arguments, the function will attempt to set
        the exposure end time from the existing start time and duration.
        
        Exposure times are stored as MJD in days.
        
        :Parameters:
        
        exposure_time: number, optional
            Exposure time (seconds).
        duration: number, optional
            The exposure duration (seconds). Defaults to the same as the exposure time.
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
        if hasattr(self, 'meta') and hasattr(self.meta, 'exposure'):
            if exposure_time is not None:
                self.meta.exposure.exposure_time =  exposure_time
            if duration is not None:
                self.meta.exposure.duration = duration
            elif exposure_time is not None:
                self.meta.exposure.duration = exposure_time
                
            if start_time == 'NOW':
                start_time = MJD_ZEROPOINT + (time.time()/SECONDS_PER_DAY)
            if start_time is not None:
                self.meta.exposure.start_time = float(start_time)
                
            if mid_time == 'NOW':
                mid_time = MJD_ZEROPOINT + (time.time()/SECONDS_PER_DAY)
            if mid_time is not None:
                self.meta.exposure.mid_time = float(mid_time)
                
            if end_time == 'NOW':
                end_time = time.time()
            elif self.meta.exposure.start_time is not None and \
                 self.meta.exposure.duration is not None and end_time is None:
                # Set the end time to start_time + duration
                end_time = self.meta.exposure.start_time + \
                    (self.meta.exposure.duration/SECONDS_PER_DAY)
            if end_time is not None:
                self.meta.exposure.end_time = float(end_time)
        else:
            strg = "Exposure metadata attributes missing from data model"
            raise AttributeError(strg)

    def get_exposure_times(self):
        """
        
        Return the exposure time metadata to the caller
        
        :Returns:
        
        (exposure_time, duration, start_time, mid_time, end_time): tuple of 5 floats
            The exposure time in seconds, clock duration in seconds, exposure
            start time, mid time and end time in MJD days. Undefined values are
            returned None.
        
        """
        exposure_time = self.meta.exposure.exposure_time
        duration = self.meta.exposure.duration
        start_time = self.meta.exposure.start_time
        mid_time = self.meta.exposure.mid_time
        end_time = self.meta.exposure.end_time
        return (exposure_time, duration, start_time, mid_time, end_time)
        
    def __str__(self):
        """
        
        Return the contents of the exposure model object as a readable
        string.
        
        """
        # First obtain a string describing the underlying ramp
        # model.
        strg = super(MiriExposureModel, self).__str__()
        
        # Also display the averaged data.
        if self.grpavg > 1 or self.intavg > 1:
            title = self.get_data_title('data')
            title += " (%d integrations and %d groups averaged)" % \
                (self.intavg, self.grpavg)
            len2 = len(title)
            title += "\n" + (len2 * '~')
            strg += title + "\n"  + str(self.data_averaged) + "\n"
            
        
        return strg

    def statistics(self, integration=None, group=-1,
                   rowstart=None, rowstop=None, colstart=None, colstop=None,
                   clip=3.0):
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
            
        :Returned:
        
        stats: str
            String describing the exposure data statistics.
            
        """
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
        return strg

    def plot_averaged(self, description='', visitor=None, **kwargs):
        """
        
        Plot the averaged data product using the algorithm contained in the
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

        # Invoke the visitor plot method on the averaged data.
        title = self.get_title()
        if description:
            title += " - " + description
        title += "\n"
        title += self.get_data_title('data')
        title += " (%d integrations and %d groups averaged)" % \
            (self.intavg, self.grpavg)
        visitor.visit(self.data_averaged, title=title, **kwargs)

    def set_exposure(self, data, dq=None):
        """
        
        Define all the exposure data in one go.
        
        :Parameters:
        
        data: numpy array
            An array containing the exposure data. It's size must match the
            number of rows, columns, groups and integrations defined
            when the ExposureDataProduct object was created.
            A 4-D array of shape (nints, ngroups, columns, rows) is
            expected, but a 3-D array of shape (nints+groups, columns,
            rows) is also acceptable and will be reshaped to match the
            4-D array.
        dq: numpy array, optional.
            An array containing the GROUPDQ data relevant for this exposure.
            It must be the same 3-D or 4-D shape as the data array.
            If include_groupdq=False was specified when the data model
            was created and this array is provided, an error is reported
            and the array is ignored.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        # The given data must be the same size as the existing SCI_data array.
        data = np.asarray(data, dtype=self.data.dtype)
        if data.size == self.data.size:
            if data.ndim == 4:
                # 4-D data is expected
                self.data = data
                if dq is not None:
                    if self.include_groupdq:
                        dq = np.asarray(dq, dtype=self.groupdq.dtype)  # Convert to same data type.
                        self.groupdq |= dq
                    else:
                        strg = "Incompatible arguments. A groupdq array is "
                        strg += "provided when include_groupdq=False. "
                        strg += "The array is ignored."
                        LOGGER.error(strg)

            else:
                # But 3-D data will be reshaped
                data.shape = (self.nints, self.ngroups, self.rows,
                              self.columns)
                self.data = data
                if dq is not None:
                    if self.include_groupdq:
                        dq = np.asarray(dq, dtype=self.groupdq.dtype)  # Make sure dq has shape attribute.
                        dq.shape = data.shape
                        self.groupdq |= dq
                    else:
                        strg = "Incompatible arguments. A groupdq array is "
                        strg += "provided when include_groupdq=False. "
                        strg += "The array is ignored."
                        LOGGER.error(strg)
            # Invalidate the averaged data
            self._data_averaged = None
        else:
            strg = "Exposure data array (%d-D) has the wrong size " % data.ndim
            strg += "(%d instead of %d)." % (data.size, self.data.size)
            raise TypeError(strg)

    def set_integration(self, data, intg, dq=None):
        """
        
        Add a new integration to the exposure data.
        
        :Parameters:
        
        data: numpy array
            An array containing the integration data to be added to the
            exposure data. It must have the same shape as one integration
            of the exposure data (i.e. a 3-D array matching the last 3
            dimensions of the exposure data).
        intg: int
            The integration to which the data belongs,
            starting from 0.
        dq: numpy array, optional.
            An array containing the GROUPDQ data relevant for this integration.
            It must be the same 3-D shape as the data array.
            If include_groupdq=False was specified when the data model
            was created and this array is provided, an error is reported
            and the array is ignored.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        # Copy the input data to a 2-D plane for this intg.
        # NOTE: This only works if data array is broadcastable so the shape
        # of the data array is checked.
        #
        data = np.asarray(data, dtype=self.data.dtype)
        if data.shape == self.data.shape[-3:]:
            self.data[intg, :, :, :] = data
            # Invalidate the averaged data
            self._data_averaged = None
            # Update the group data quality array if necessary.
            if dq is not None:
                if self.include_groupdq:
                    dq = np.asarray(dq, dtype=self.groupdq.dtype)  # Convert to same data type.
                    self.groupdq[intg, :, :, :] |= dq
                else:
                    strg = "Incompatible arguments. A groupdq array is "
                    strg += "provided when include_groupdq=False. "
                    strg += "The array is ignored."
                    LOGGER.error(strg)
        else:
            strg = "Integration data array has the wrong shape "
            strg += "(%s compared with %s)." % (str(data.shape),
                                                str(self.data.shape))
            raise TypeError(strg)

    def get_integration(self, intg):
        """
        
        Return the 3-D data belonging to the specified integration.
        
        :Parameters:

        intg: int
            The integration required, starting from 0.
        
        """
        return self.data[intg, :, :, :]

    def set_group(self, data, group, intg, dq=None):
        """
        
        Add a new group to the exposure data.
        
        :Parameters:
        
        data: numpy array
            An array containing the readout data to be added to the
            exposure data. It must have the same shape as one frame of
            the exposure data (i.e. a 2-D array with the same number of
            rows and columns specified in the MiriExposureData
            constructor).
        group: int
            The group to which the data readout belongs,
            starting from 0.
        intg: int
            The integration to which the data readout belongs,
            starting from 0.
        dq: numpy array, optional.
            An array containing the GROUPDQ data relevant for this group.
            It must be the same 2-D shape as the data array.
            If include_groupdq=False was specified when the data model
            was created and this array is provided, an error is reported
            and the array is ignored.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        # TODO: Include a 2-D DQ array to be combined with the GROUPDQ array
        #
        # Copy the input data to a 2-D plane for this group/intg combination.
        # NOTE: This only works if data array is broadcastable so the shape
        # of the data array is checked.
        #
        data = np.asarray(data, dtype=self.data.dtype)
        detector_shape = (self.rows, self.columns)
        if data.shape == detector_shape:
            self.data[intg, group, :, :] = data  
            # Invalidate the averaged data
            self._data_averaged = None
            # Update the group data quality array if necessary.
            if dq is not None:
                if self.include_groupdq:
                    dq = np.asarray(dq, dtype=self.groupdq.dtype)  # Convert to same data type.
                    self.groupdq[intg, group, :, :] |= dq
                else:
                    strg = "Incompatible arguments. A groupdq array is "
                    strg += "provided when include_groupdq=False. "
                    strg += "The array is ignored."
                    LOGGER.error(strg)
        else:
            strg = "Group data array has the wrong shape "
            strg += "(%s instead of %s)." % (str(data.shape),
                                             str(detector_shape))
            raise TypeError(strg)

    def get_group(self, group, intg):
        """
        
        Return the 2-D data belonging to the specified group and
        integration.
        
        :Parameters:

        group: int
            The group required, starting from 0.
        intg: int
            The integration required, starting from 0.
        
        """
        return self.data[intg, group, :, :]
    
    def add_dark(self, darkarray, clipvalue=65535.0 ):
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
        clipvalue: float (optional)
            An upper limit to which data are clipped after adding the DARK.
            The default value is 65535.0, which keeps the exposure data within
            the range of 16-bit telemetry. Set to None to turn off the
            clipping.
            
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
            
    def apply_translation(self, translation_table, fromcolumn=None,
                          tocolumn=None, clipvalue=65535.0 ):
        """
        
        Apply a translation table to the exposure data in-situ (for simulation
        purposes). This function may be used to apply a linearity correction
        table derived from a MIRI LINEARITY CDP file.
        The exposure data and translation table are assumed both to be in
        DN units.

        The translation table is contained in an array of integers so that
        
           output_dn = translation_table[ input_dn ]
            
        :Parameters:
        
        translation_table: array of int
            An array of integers containing the translation table.
        fromcolumn: int (optional)
            If given, the starting row to be converted
        tocolumn: int (optional)
            If given, the finishing row to be converted
        clipvalue: float (optional)
            An upper limit to which data are clipped after adding the DARK.
            The default value is 65535.0, which keeps the exposure data within
            the range of 16-bit telemetry. Set to None to turn off the
            clipping.
            
        :Raises:
    
        TypeError
            Raised if the data array has the wrong type, size or shape.
            
        """
        translation_table = np.asarray(translation_table)
        if translation_table.ndim == 1:
            maxtable = len(translation_table)
            if fromcolumn is None:
                fromcolumn = 0
            if tocolumn is None:
                tocolumn = self.data.shape[3]
            LOGGER.debug("Applying linearity translation table of length %s" % maxtable + \
                  " to columns %d-%d" % (fromcolumn,tocolumn))
            selection = np.clip(self.data[:, :, :, fromcolumn:tocolumn].astype(int), 0, maxtable)
            self.data[:, :, :, fromcolumn:tocolumn] = translation_table[selection]
        else:
            strg = "Translation table array must be 1-D"
            raise TypeError(strg)
    
        # The translation table must not allow the resulting data to go
        # negative or contain zeros. It must also not allow the data to go
        # above the maximum value for 16-bit telemetry data.
        whereneg = np.where( self.data < 1.0 )
        if whereneg and (len(whereneg[0]) > 0):
            strg = "%d negative pixels after translation." % len(whereneg[0])
            LOGGER.debug(strg)
            self.data[whereneg] = 1.0
        if clipvalue is not None:
            whereclipped = np.where( self.data > clipvalue )
            if whereclipped and (len(whereclipped[0]) > 0):
                strg = "%d saturated pixels after translation." % len(whereclipped[0])
                LOGGER.debug(strg)
                self.data[whereclipped] = clipvalue 

    def save(self, path, *args, **kwargs):
        """
        
        Override the save method of DataModel to allow the data
        to be averaged before saving.

        :Parameters:
        
        path : str
            The path to the file that we're about to save to.
        \*args, \*\*kwargs:
            Any additional arguments are passed along to
            the save function (e.g. overwrite=True).
        
        """
        # TODO: Is it still correct to hold back on the averaging
        # until the data are saved to disk? What about when plotting
        # the data? Should the data model hold the averaged data and
        # more sophisticated accumulation algorithms be inserted into
        # the "set_" methods?

        # Decide whether the SCI array contents need to be
        # averaged before saving.
        if self.grpavg <= 1 and self.intavg <= 1:
            # Ensure the reference output data has been extracted.
            self._extract_refout()
            
            # Explicitly delete the data quality arrays before saving
            # if they are not needed.
            if not self.include_pixeldq:
                if hasattr(self, 'pixeldq'):
                    del self.pixeldq
                if hasattr(self, 'pixeldq_def'):
                    del self.pixeldq_def
            if not self.include_groupdq:
                if hasattr(self, 'groupdq'):
                    del self.groupdq
                if hasattr(self, 'groupdq_def'):
                    del self.groupdq_def
            
            # No averaging - skip straight to the parent save method.
            super(MiriExposureModel, self).save(path, *args, **kwargs)
        else:
            # Make a copy of the exposure data which is reduced in size by
            # averaging the groups and/or integrations and save that
            # product instead.
            with MiriExposureModel(self.rows, self.columns,
                                   self.ngroups_file, self.nints_file,
                                   self.readpatt, self.refrows, self.refcolumns,
                                   nframes=self.nframes, groupgap=self.groupgap,
                                   include_pixeldq=self.include_pixeldq,
                                   include_groupdq=self.include_groupdq,
                                   pixeldq=self.pixeldq, maskwith=self.maskwith) \
                    as newproduct:
                # Copy the metadata and update the number of groups and
                # integrations.
                newproduct.copy_metadata( self, ignore=['NINTS', 'NGROUPS'] )
                newproduct.meta.exposure.nints = self.nints
                newproduct.meta.exposure.ngroups = self.ngroups
                newproduct.meta.exposure.groups_averaged = self.grpavg
                newproduct.meta.exposure.integrations_averaged = self.intavg
                newproduct.set_exposure(self.data_averaged)
                # Explicitly delete the data quality arrays before saving
                # if they are not needed.
                if not self.include_pixeldq:
                    if hasattr(newproduct, 'pixeldq'):
                        del newproduct.pixeldq
                    if hasattr(newproduct, 'pixeldq_def'):
                        del newproduct.pixeldq_def
                if not self.include_groupdq:
                    if hasattr(newproduct, 'groupdq'):
                        del newproduct.groupdq
                    if hasattr(newproduct, 'groupdq_def'):
                        del newproduct.groupdq_def
                super(MiriExposureModel, newproduct).save(path, *args, **kwargs)
                del newproduct

    def slope_data(self, grptime=2.77504, diff_only=False):
        """
        
        Return a copy of the SCI data array where all group planes
        have been fitted with a straight line to make slope data.

        NOTE: This function is very inefficient!
        
        :Parameters:
        
        grptime: float, optional, default=2.77504s
            The time interval between groups (which determines the
            time axis for the slope calculation).
        diff_only: bool, optional, default=False
            If True, implement a quick and dirty estimate of the
            slope by subtracting the last frame from the first.
            
        :Returned:
        
        slope_data: array_like float32
            The SCIdata array with groups converted to slopes.
            The is a slope image for each integration.
            It can be 2-D or 3-D.
            
        """
        assert self.ngroups > 1
        # TODO: Very inefficient and does not calculate ERR or DQ arrays.
        output = np.zeros([self.nints, self.rows, self.columns],
                            dtype=self.data.dtype)
        
        if diff_only:
            # Quick and dirty estimate which subtracts the last
            # ramp from the first. timediff should never be zero because
            # self.ngroups is forced to be > 1.
            timediff = grptime * (self.ngroups - 1)
            output = (self.data[:,-1,:,:] - self.data[:,0,:,:]) / float(timediff)
        else:
            # Full straight line fit
            timearray = grptime * np.array( list(range(0, self.ngroups)) )
            for intg in range(0, self.nints):
                for row in range(0, self.rows):
                    for column in range(0, self.columns):
                        # Ouch! Slope calculated one (row,column) at a time.
                        # Can the efficiency be improved?
                        (slope, ic) = linear_regression( timearray,
                                            self.data[intg,:,row,column] )
                        output[intg,row,column] = slope
        return output

    def cube_data(self):
        """
        
        Generate a copy of the SCI data array which has been changed in
        shape from a hypercube to a cube, by combining all groups and
        integrations together in the third dimension.
        
        :Parameters:
        
        None.
            
        :Returned:
        
        cube_data: array-like float32
            The exposure data array changed into a cube
            
        """
        cube_data = copy.deepcopy(self.data)
        cube_data.shape = [self.nints * self.ngroups, self.rows, self.columns]
        return cube_data
    
    def _extract_refout(self):
        """
        
        Extract the reference output data from the SCI data array
        and write it to the REFOUT data array. Also, remove the
        reference output from the PIXELDQ array.
        
        This function can only be executed once.
        
        """
        if not self._refout_extracted and \
           self.refcolumns > 0 and self.refrows > 0:
            top_rows = int(self.refrows * self.refcolumns / self.columns)           
            self.refout = copy.deepcopy(self.data[:, :, -top_rows:, :])
            self.refout.shape = [self.nints, self.ngroups, self.refrows,
                                 self.refcolumns]
            
            new_data = self.data[:, :, :-top_rows, :]
            self.data = new_data
            
            if self.include_pixeldq and self.pixeldq is not None:
                new_pixeldq = self.pixeldq[:-top_rows, :]
                self.pixeldq = new_pixeldq

            if self.include_groupdq and self.groupdq is not None:
                new_groupdq = self.groupdq[:, :, :-top_rows, :]
                self.groupdq = new_groupdq
            
            self._refout_extracted = True

    def _average_data(self):
        """
        
        Generate a copy of the SCI data array which has been reduced
        in size by averaging integrations and groups, as prescribed
        by the intavg and grpavg parameters. The averaged data is only
        generated once and stored internally, unless a recalculation
        is forced with recalculate=True.
        
        :Parameters:
        
        None.
            
        :Returned:
        
        average_data: array-like float32
            The averaged exposure data array.
            
        """
        # TODO: This function looks inefficient - can it be improved?
        output = np.zeros([self.nints_file, self.ngroups_file,
                           self.rows, self.columns], dtype=self.data.dtype)
        count = np.zeros([self.nints_file, self.ngroups_file, 1, 1],
                         dtype=self.data.dtype)

        for intg in range(0, self.nints):
            intg_file = intg // self.intavg
            for grp in range(0, self.ngroups):
                grp_file = grp // self.grpavg
                output[intg_file, grp_file, :, :] += self.data[intg, grp, :, :]
                count[intg_file, grp_file, 0, 0] += 1.0

            # Avoid divide by zero.
            iszero = np.where(count < 1.0)
            if iszero:
                count[iszero] = 1.0
            output = output / count

        del count            
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

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriExposureModel module.")

    PLOTTING = False
    SAVE_FILES = False
    INCLUDE_DQ = True

    data4x5 = np.array([[1.,2.,3.,4.,5.],
                        [5.,6.,7.,8.,9.],
                        [8.,7.,6.,5.,6.],
                        [4.,3.,2.,1.,2.]])
    dq4x5 = np.array([[0,0,0,1,0],
                      [0,1,0,0,1],
                      [0,0,1,0,0],
                      [1,0,0,0,0]])

    print("Exposure data with 2 integrations and 3 groups of 4x4 pixels:")
    with MiriExposureModel(readpatt='FAST', nints=2, ngroups=3, rows=4,
                           columns=5, refrows=0, refcolumns=0, pixeldq=dq4x5,
                           include_pixeldq=INCLUDE_DQ,
                           include_groupdq=INCLUDE_DQ) \
            as testdata1:
        testdata1.set_instrument_metadata(detector='MIRIMAGE', filt='F560W',
                                channel='', ccc_pos='OPEN', 
                                deck_temperature=10.0,
                                detector_temperature=7.0)
        testdata1.set_simulation_metadata( 'SOLAR_MIN' )
        testdata1.set_exposure_times(exposure_time=3.0, start_time='NOW')
        testdata1.set_exposure_type()
        # Build the exposure one group at a time
        for intg in range(0, 2):
            for group in range(0, 3):
                data = data4x5 + group + intg*group
                testdata1.set_group(data, group, intg)
        
        print(testdata1)
        print(testdata1.statistics())
        if PLOTTING:
            testdata1.plot(description="testdata1")
            testdata1.plot_averaged(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_exposure_model1.fits", overwrite=True)
            
        print("Testing the translation table function")
        table_l = np.array([13,12,11,10,9,8,7,6,5,4,3,2,1,0])
        table_r = np.array([130,120,110,100,90,80,70,60,50,40,30,20,10,0])
        testdata1.apply_translation( table_l, fromcolumn=0, tocolumn=2 )
        testdata1.apply_translation( table_r, fromcolumn=2, tocolumn=5)
        print(testdata1)
        print(testdata1.statistics())
        del testdata1

    print("Exposure data with 1 integration and 12 groups of 4x4 pixels" + \
          " (and every 4 groups averaged):")
    with MiriExposureModel(readpatt='FAST', nints=2, ngroups=12, rows=4,
                           grpavg=4, columns=5, refrows=0, refcolumns=0,
                           pixeldq=dq4x5, include_pixeldq=INCLUDE_DQ,
                           include_groupdq=INCLUDE_DQ) \
            as testdata2:
        testdata2.set_simulation_metadata( 'NONE' )
        testdata2.set_exposure_type()
        # Build the exposure one group at a time
        for intg in range(0, 2):
            for group in range(0, 12):
                data = data4x5 + group + intg*group
                testdata2.set_group(data, group, intg)
         
        print(testdata2)
        print(testdata2.statistics())
        if PLOTTING:
            testdata2.plot(description="testdata2")
            testdata2.plot_averaged(description="testdata2")
        if SAVE_FILES:
            testdata2.save("test_exposure_model2.fits", overwrite=True)
        del testdata2

    print("Test finished.")
