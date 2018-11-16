#!/usr/bin/env python

"""

Module miri_measurement - Contains the MiriMeasurement class and 
associated functions, which are used to store a set of measurements
which vary with a variable.

:History:

14 Jul 2010: Original MeasuredVariable product created (in measured_variable.py)
10 Jul 2013: Modified to work with jsonschema draft 4 schemas.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
27 May 2015: Replaced pyfits with astropy.io.fits
01 Apr 2016: Changed to YAML schemas.
31 Aug 2016: Removed obsolete code from table setting. Read the test tables
             from FITS files, not ASCII files. ASCII support withdrawn.
15 Jun 2017: TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema. Changed 'ANY' to 'N/A' in test data.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


#import sys
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    pyfits = None

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel
#from miri.datamodels.plotting import DataModelPlotVisitor

# List all classes and global functions here.
__all__ = ['MiriMeasurement']


class MiriMeasurement(MiriDataModel):
    """
    
    A generic data model for a MIRI measurement, based on the MIRI
    base model, MiriDataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the flux parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    measurement_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (variable:number,
        measurement:number), giving the transmission of the filter
        as a function of wavelength.
        Or: A numpy record array containing the same information as above.
        A transmission table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
        NOTE: The table MUST contain 1 variable column plus 5 measurement
        columns.
    title: str (optional)
        A title describing the measurement.
    name: str (optional)
        The name of the property being measured.
        If specified, the name should be applicable to all the columns in
        the values array.
    unit: str (optional)
        The units of the property being measured.
    vname: str (optional)
        The name of the variable the property is measured against.
    vunit: str  (optional)
        The units of the variable.
    mcolumns: int (optional), default=5
        The number of valid measurement columns contained in the table
        (not including the column containing the variable parameter).
        NOTE: The yaml schema defines 5 columns, so 5 columns will
        always be written to the data structure. This parameter declares
        how many of those columns are valid. The rest are padding.
    interptype: str (optional), default='LINEAR'
        The type of interpolation to be used between the data points.
        Possible values are:
        
        * 'LINEAR' - Apply linear interpolation
        * 'LOGLIN' - Interpolate the logarithm of the variables
                     against the unchanged adjustable parameter.
        * 'LINLOG' - Interpolate the variables linearly against the
                     logarithm of the adjustable parameter.
        * 'LOGLOG' - Interpolate the logarithm of the variables
                     against the logarithm of the adjustable parameter.
                     
        If not specified will default to 'LINEAR'.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_measurement.schema.yaml"
    fieldnames = ('VARIABLE', 'MEASUREMENT1', 'MEASUREMENT2', 'MEASUREMENT3',
                  'MEASUREMENT4', 'MEASUREMENT5')
    
    def __init__(self, init=None, measurement_table=None, title=None,
                 name=None, unit=None, vname=None, vunit=None, mcolumns=5,
                 interptype=None, **kwargs):
        """
        
        Initialises the MiriFilter class.
        
        Parameters: See class doc string.

        """
        super(MiriMeasurement, self).__init__(init=init, title=title, **kwargs)

        # Data type is measurement.
        self.meta.model_type = 'Measurement'
        
        if name is not None:
            self.meta.measurement.name = name
        if unit is not None:
            self.meta.measurement.unit = unit
        if vname is not None:
            self.meta.measurement.variable = vname
        if vunit is not None:
            self.meta.measurement.vunit = vunit
            self.meta.measurement_table.tunit1 = vunit
        if unit is not None:
            self.meta.measurement_table.tunit2 = unit
            self.meta.measurement_table.tunit3 = unit
            self.meta.measurement_table.tunit4 = unit
            self.meta.measurement_table.tunit5 = unit
        self.meta.measurement_table.mcolumns = mcolumns
        if interptype is None or not interptype:
            interptype = 'LINEAR'
        self.meta.measurement.interp = interptype

        if measurement_table is not None:
            try:
                #phot_table = np.recarray(phot_table)
                self.measurement_table = measurement_table
            except (ValueError, TypeError) as e:
                strg = "measurement_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # Copy the table column units from the schema, if defined.
        measurement_units = self.set_table_units('measurement_table')

        # Cached arrays
        self._parameters = None
        self._values = None

    def get_table_arrays(self):
        """
        
        Helper function to extract the data arrays from the
        measurements table.
        
        :Returns:
        
        (parameters, val1, val2, val3, val4, val5): tuple of 6 numpy ndarrays
        
        """
        #ftable = np.asarray(self.measurement_table)
        ftable = self.measurement_table
        parameters = []
        val1 = []
        val2 = []
        val3 = []
        val4 = []
        val5 = []
        for item in ftable:
            parameters.append(item[0])
            val1.append(item[1])
            val2.append(item[2])
            val3.append(item[3])
            val4.append(item[4])
            val5.append(item[5])
        parameters = np.asarray(parameters)
        val1 = np.asarray(val1)
        val2 = np.asarray(val2)
        val3 = np.asarray(val3)
        val4 = np.asarray(val4)
        val5 = np.asarray(val5)
        return (parameters, val1, val2, val3, val4, val5)

    def lookup(self, pvalue, column=0):
        """
        
        Lookup the value of a particular instance of this variable.
        
        :Parameters:
        
        pvalue: float (or array_like)
            The value(s) of the adjustable parameter.
        column: int
            The index of the measurement column to be looked up
            (starting from 0)
        
        :Returned:
        
        value: float 
            The interpolated variable value.
        
        """
        if column < 0 or (self.meta.measurement_table.mcolumns is not None and column >= self.meta.measurement_table.mcolumns):
            strg = "Specified value column, %d, is out of range (0-%d)" % \
                (column, self.meta.measurement_table.mcolumns)
            raise ValueError(strg)
        parameters = self.parameters
        values = self.values[column]
        if self.meta.measurement.interp == 'LOGLIN':
            # The measured values are interpolated logarithmically.
            logvals = np.log10(values)
            logres = np.interp(float(pvalue), parameters, logvals)
            result = 10.0 ** logres
            del logvals
        elif self.meta.measurement.interp == 'LINLOG':
            # The varied parameter is interpolated logarithmically.
            logp = np.log10(float(pvalue))
            logparams = np.log10(parameters)
            result = np.interp(logp, logparams, values)
            del logparams
        elif self.meta.measurement.interp == 'LOGLOG':
            # Both values are interpolated logarithmically.
            logp = np.log10(float(pvalue))
            logvals = np.log10(values)
            logparams = np.log10(parameters)
            logres = np.interp(logp, logparams, logvals)
            result = 10.0 ** logres
            del logvals, logparams
        else:
            # Both values are interpolated linearly.
            result = np.interp(float(pvalue), parameters, values) 
        return result

    def get_comments(self):
        """
        
        Return a comment string describing the data object.
        This function exists for backwards compatibility.
        
        """
        strg = "Measurement of %s against %s in %d columns" % \
            (self.meta.measurement.name, self.meta.measurement.variable,
             self.meta.measurement_table.mcolumns)
        return strg

#     def get_comments(self):
#         """
#         
#         Return a list of comments describing the measured data,
#         formatted so they may be inserted in a FITS header.
#             
#         :Returns:
#         
#         comments: list of str
#             A list of comment strings describing this MeasuredVariable
#             object.
#             
#         """
#         _MAX_COMMENT_LEN = 72
#         comments = []
#         strg = self.get_title(underline=False)
#         if 'VERSION' in self.metadata:
#             strg += ": %s" % self.metadata['VERSION']
#         comments.append(strg[:_MAX_COMMENT_LEN])
#         
#         if self.filename is not None:
#             if self.filename:
#                 if self.name:
#                     MAX_LEN = _MAX_COMMENT_LEN - len(self.name) - 8
#                     if len(self.filename) > MAX_LEN:
#                          The filename is truncated if it is too long
#                          for one comment string.
#                         MAX_LEN -= 3
#                         strg = "%s file: ...%s" % \
#                             (self.name, self.filename[-MAX_LEN:])
#                     else:
#                         strg = "%s file: %s" % \
#                             (self.name, self.filename)
#                 elif isinstance(self.vnames, str):
#                     MAX_LEN = _MAX_COMMENT_LEN - len(self.vnames) - 8
#                     if len(self.filename) > MAX_LEN:
#                          The filename is truncated if it is too long
#                          for one comment string.
#                         MAX_LEN -= 3
#                         strg = "%s file: ...%s" % \
#                             (self.vnames, self.filename[-MAX_LEN:])
#                     else:
#                         strg = "%s file: %s" % \
#                             (self.vnames, self.filename)
#                 else:
#                     MAX_LEN = _MAX_COMMENT_LEN - 8
#                     if len(self.filename) > MAX_LEN:
#                          The filename is truncated if it is too long
#                          for one comment string.
#                         MAX_LEN -= 3
#                         strg = "File: ...%s" % self.filename[-MAX_LEN:]
#                     else:
#                         strg = "File: %s" % self.filename
#                     
#                 comments.append(strg[:_MAX_COMMENT_LEN])
#         
#         if self.comment:
#             comments.append(self.comment[:_MAX_COMMENT_LEN])
#       
#         return comments
    
    def plot(self, plotfig=None, plotaxis=None, columns=None, overlay=False,
             xlabel='', ylabel='', title=None, description='', **kwargs):
        """

        Plot one or more columns of measured variables against a
        parameter.

        A plot of one column is added to a single matplotlib axis.
        If more than one column is plotted, the function creates
        a new matplotlib figure and plots exach column in a separate
        axis.
        
        :Parameters:
        
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
            It is only sensible to supply an axis when one column
            is being plotted or when overlay=True.
        columns: int or list of int (optional)
            Either: An integer containing the column number of be plotted
            Or: A list of integers giving column numbers to be plotted.
            If not specified, all columns will be plotted.
        overlay: bool (optional)
            When True, multi-column plots are overlayed together onto a
            single plot. The default is False.
        xlabel: str (optional)
            An X label for the plot. The default is the name of the
            parameter.
        ylabel: str (optional)
            A Y label for the plot. The default is the name of
            each variable.
        title: string (optional)
            Optional title for the plot. The default is a string
            describing the MeasuredVariable object.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
        description: str (optional)
            Optional description to be appended to the plot title, if
            required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
        \*\*kwargs:
            All other keyword arguments will be passed to matplotlib.plot
            For example label, linewidth, linestyle, color, marker, etc...
            See the matplotlib documentation for a full list.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot

        """
        # Begin by constructing the overall plot title and X label,
        # which no not change from column to column.
        if title is None:
            if title is None:
                tstrg = self.get_title(underline=False)
            else:
                tstrg = self.title
        else:
            tstrg = title
        if description:
            tstrg += " - " + description

        if xlabel:
            xstrg = xlabel
        else:
            xstrg = self.meta.measurement.variable
            if self.meta.measurement.vunit:
                xstrg += " (%s)" % self.meta.measurement.vunit

        # Match the scale of the axes to the interpolation type.
        xscale = 'linear'
        yscale = 'linear'
        if self.meta.measurement.interp == 'LOGLIN':
            yscale = 'log'
        elif self.meta.measurement.interp == 'LINLOG':
            xscale = 'log'
        elif self.meta.measurement.interp == 'LOGLOG':
            xscale = 'log'
            yscale = 'log'
        
        # Get the parameters array, which doesn't change from column
        # to column.    
        parameters = self.parameters
        
        if columns is None:
            columns = list(range(0,self.meta.measurement_table.mcolumns))
        elif isinstance(columns, (int,float)):
            columns = [int(columns)]
        # If there are 2 or more columns, produce multiple xy plots in
        # a column. A single column will generate 1 plot.
        if len(columns) > 1:
            # For each column, build up a list of values arrays and
            # a list of Y labels.
            valuelist = []
            ylabels = []
            for column in columns:
                if ylabel:
                    ystrg = ylabel
                else:
                    if self.meta.measurement.name:
                        ystrg = "%s %d" % (self.meta.measurement.name, column)
                    else:
                        ystrg = "column %d" % column
                    if self.meta.measurement.unit:
                        if overlay:
                            ystrg += " (%s)" % self.meta.measurement.unit
                        else:
                            ystrg += "\n(%s)" % self.meta.measurement.unit
                ylabels.append(ystrg)
        
                values = self.values[column]
                valuelist.append(values)

            if overlay:
                # Turn the list of Y values into a 2-D array and
                # generate a single plot showing all the columns
                # overlayed on top of oneanother. The Y labels are
                # joined into a single string.
                varray = np.transpose(valuelist)
                ystrg = "\n".join(ylabels)
                mplt.plot_xy(parameters, varray, plotfig=plotfig, plotaxis=plotaxis,
                             xscale=xscale, yscale=yscale, xlabel=xstrg,
                             ylabel=ystrg, title=tstrg, **kwargs)
                
            else:
                # Generate an XY column plot.
                mplt.plot_xycolumn(parameters, valuelist, xscale=xscale,
                               yscale=yscale, xlabel=xstrg, ylabels=ylabels,
                               title=tstrg, **kwargs)
        else:
            # For a single column there is only 1 values array and 1 Y label.
            column = columns[0]
            values = self.values[column]
            if ylabel:
                ystrg = ylabel
            else:
                if self.meta.measurement.name:
                    ystrg = self.meta.measurement.name
                else:
                    ystrg = "Column %d" % column
                if self.meta.measurement.unit:
                    ystrg += " (%s)" % self.meta.measurement.unit
            mplt.plot_xy(parameters, values, plotfig=plotfig, plotaxis=plotaxis,
                         xscale=xscale, yscale=yscale, xlabel=xstrg,
                         ylabel=ystrg, title=tstrg, **kwargs)
#         mplt.close()

    def __str__(self):
        """
        
        Returns a string representation of the object MiriMeasurement
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()

        # Describe the filter transmission table
        strg += "\nColumn names: " + str(self.fieldnames) + "\n"
        if self.measurement_table is not None:
            strg += self.get_data_str('measurement_table', underline=True,
                                      underchar="-")
        return strg

    @property
    def parameters(self):
        if self._parameters is None:
            (parameters,val1,val2,val3,val4,val5) = self.get_table_arrays()
            self._parameters = parameters
            self._values = [val1,val2,val3,val4,val5]
            #self.meta.measurement_table.mcolumns = len(self._values)
            return self._parameters
        else:
            return self._parameters

    @parameters.setter
    def parameters(self, data):
        raise AttributeError("parameters attribute is read-only")

    @property
    def values(self):
        if self._values is None:
            (parameters,val1,val2,val3,val4,val5) = self.get_table_arrays()
            self._parameters = parameters
            self._values = [val1,val2,val3,val4,val5]
            #self.meta.measurement_table.mcolumns = len(self._values)
            return self._values
        else:
            return self._values

    @values.setter
    def values(self, data):
        raise AttributeError("values attribute is read-only")

#     # The parameter arrays is accessible through the parameters attribute,
#     # or through the get_parameters and set_parameters methods.
#     @property
#     def parameters(self):
#         return getattr(self.measurements.tabledata, self.meta.measurement.variable)
# 
#     @parameters.setter
#     def parameters(self, data):
#         setattr(self.measurements.tabledata, self.meta.measurement.variable, data )
#         self.modified = True
# 
#     def get_parameters(self):
#         """
#         
#         Return the parameter values as a floating point array.
#         
#         :Parameters:
#         
#         None.
#         
#         :Returns:
#         
#         parameters: numpy array (float)
#             The array of parameter values.
#         
#         """
#         params = getattr(self.measurements.tabledata, self.meta.measurement.variable)
#         return np.asarray(params, dtype=np.float)
# 
#     def set_parameters(self, data):
#         """
#         
#         Write a new set of parameter values.
#         
#         :Parameters:
#         
#         parameters: numpy array (float)
#             The array of parameter values.
#         
#         """
#         setattr(self.measurements.tabledata, self.meta.measurement.variable, data )
#         self.modified = True
# 
#     def get_columns(self):
#         """
#         
#         Return the number of columns of values
#         
#         """
#         return self.meta.measurement_table.mcolumns
# 
#     # The values array can be accessed one column at a time.
#     def get_values(self, column):
#         """
#         
#         Return a set of variable values as a floating point array
#         
#         :Parameters:
#         
#         column: int
#             The column whose values are to be returned.
#         
#         :Returns:
#         
#         values: numpy array (float)
#             The array of variable values for the given column.
#         
#         """
#         if column < 0 or vcolumn >= self.meta.measurement_table.mcolumns:
#             strg = "Specified value column, %d, is out of range (0-%d)" % \
#                 (vcolumn, self.meta.measurement_table.mcolumns)
#             raise ValueError(strg)
#         columndata = self.measurements.get_column( vcolumn+1 )
#         return np.asarray(columndata, dtype=np.float)
#     
#     def set_values(self, vcolumn, data):
#         """
#         
#         Write a new column of variable values.
#         
#         :Parameters:
#         
#         vcolumn: int
#             The column whose values are to be updated.
#         values: numpy array (float)
#             The array of variable values for the given column.
#         
#         """
#         if vcolumn < 0 or vcolumn >= self.meta.measurement_table.mcolumns:
#             strg = "Specified value column, %d, is out of range (0-%d)" % \
#                 (vcolumn, self.meta.measurement_table.mcolumns)
#             raise ValueError(strg)
#         self.measurements.set_column( vcolumn+1, data )
#         self.modified = True

def ascii_to_measurement(filename, title=None, name=None, unit=None, vname=None,
                         vunit=None, interptype=None, **kwargs):
    """
        
    Create a MiriMeasurement data product from an ASCII file.
                
    :Parameters:
    
    filename: str
        The name of an ASCII file from which to create the Filter Data Product.

    """
    strg = "Reading a MiriMeasurement model from an ASCII file "
    strg += "is not longer supported."
    raise NotImplementedError(strg)
#     # This function no longer works
#     try:
#         data = np.loadtxt(filename)
# 
#     except Exception as e:
#         # If the file could not be read re-raise the exception
#         # with a more meaningful error message.
#         strg = "%s: Could not read Measurement Data Product ASCII file.\n   %s" % \
#             (e.__class__.__name__, e)
#         raise IOError(strg)
# 
#     mcolumns = len(data[0, :]) - 1
#     padding = 5 - mcolumns
#  
#     variable = data[:, 0]
#     meas = []
#     for col in range(1,mcolumns+1):
#         measurement = data[:, col]
#         meas.append(measurement)
#   
#     measurement_table = []
#     for ii in range(0, len(variable)):
#         record = [variable[ii]]
#         for jj in range(0,mcolumns):
#             record += [meas[jj][ii]]
#         # If there are too few columns, extra empty ones need to be added
#         for kk in range(0,padding):
#             record += [1.0]
#         measurement_table.append(record)
# 
#     # Initialise the data product from the table just created.
#     miri_measurement = MiriMeasurement( measurement_table=data,
#                                         title=title, name=name, unit=unit,
#                                         vname=vname, vunit=vunit,
#                                         mcolumns=mcolumns, interptype=interptype,
#                                         **kwargs)
#     
#     return miri_measurement

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see
# miri.tools/tests/test_measured_variable.py.
#
if __name__ == '__main__':
    print("Testing the MiriMeasurement class")
    import os
    import miri.simulators.data as simdata
    sim_datapath = simdata.__path__[0]

    PLOTTING = False       # Set to False to turn off plotting.
    SAVE_FILES = False     # Set to False to avoid saving files.
    
    import miri.datamodels.data
    datapath = miri.datamodels.data.__path__[0]
    test_file_name = os.path.join(datapath, "example_measurement.fits")
    
    title = "Read noise for channels 1 to 5"
    name= 'Read noise'
    unit = 'electrons'
    vname = 'Temperature'
    vunit = 'K'
    mv = MiriMeasurement(test_file_name, title=title, name=name, unit=unit,
                              vname=vname, vunit=vunit)
#     mv = ascii_to_measurement(test_file_name, title=title, name=name, unit=unit,
#                               vname=vname, vunit=vunit)
    print(mv)
    print(mv.get_comments())
    if PLOTTING:
        mv.plot(description='\nplot of mv in columns')
        mv.plot(overlay=True, description='\nplot of mv overlayed')
        
    print("Interpolated value for column 0 at 6.7 is", mv.lookup(6.7, 0))
    print("Interpolated value for column 3 at 6.7 is", mv.lookup(6.7, 3))
    if SAVE_FILES:
        mv.save('example_measurement.fits', overwrite=True)
        readback = MiriMeasurement('example_measurement.fits')
        print("Read back from FITS file")
        print(readback)

    mv2 = MiriMeasurement(test_file_name, title=title, name=name, unit=unit,
                vname=vname, vunit=vunit, interptype='LOGLIN')
#     mv2 = ascii_to_measurement(test_file_name, title=title, name=name, unit=unit,
#                 vname=vname, vunit=vunit, interptype='LOGLIN')
    if PLOTTING:
        mv2.plot(description='\nplot of mv2 in LOGLIN space')
    print("Interpolated value for column 0 at 6.7 is", mv2.lookup(6.7, 0))
    print("Interpolated value for column 3 at 6.7 is", mv2.lookup(6.7, 3))

    dark_file_name = os.path.join(datapath, "example_single_measurement.fits")
    mv3 = MiriMeasurement(dark_file_name, title="Single measurement column",
                               name='Dark Current', unit=unit,
                               vname=vname, vunit=vunit, interptype='LOGLIN')
#     mv3 = ascii_to_measurement(dark_file_name, title="Single measurement column",
#                                name='Dark Current', unit=unit,
#                                vname=vname, vunit=vunit, interptype='LOGLIN')
    print(mv3)
    print(mv3.get_comments())
    if PLOTTING:
        mv3.plot(description='\nplot of mv3 dark current')
    print("Interpolated value at 6.7 is", mv.lookup(6.7))
    if SAVE_FILES:
        mv3.save('example_single_measurement_out.fits', overwrite=True)

    del mv, mv2, mv3
    
# -----------------------------------------------------------------------

    # NOTE: The following tests assume the miri.simulators package has been
    # built and the ASCII files in the following list of ASCII files are
    # available.
    # The following code (with SAVE_FILES=True) can be used to recreate the
    # FITS versions of these simulator files.
    print("\nMeasurements read from simulator ASCII files:")
    fits_file_names = [os.path.join(sim_datapath, "amplifiers/read_noiseIM"),
                        os.path.join(sim_datapath, "amplifiers/read_noiseLW"),
                        os.path.join(sim_datapath, "amplifiers/read_noiseSW")]
    title = "Amplifier read noise for channels 1 to 5"
    name = "Read noise"
    unit = "electrons"
    vname = "Temperature"
    vunit = "K"
    detectors = ['MIRIMAGE', 'MIRIFULONG', 'MIRIFUSHORT']
    ii = 0
    for filnam in fits_file_names:
        detector = detectors[ii]
        full_filename = filnam + ".fits"
        print("Reading read noise for", detector, "from", full_filename)
        meas1 = MiriMeasurement(full_filename, title=title, name=name,
                                    unit=unit, vname=vname, vunit=vunit,
                                    interptype='LINEAR')
#         meas1 = ascii_to_measurement(full_filename, title=title, name=name,
#                                     unit=unit, vname=vname, vunit=vunit,
#                                     interptype='LINEAR')
        meas1.set_instrument_metadata(detector=detector, modelnam='FM',
                                filt='N/A', channel='N/A', band='N/A',
                                ccc_pos='OPEN', deck_temperature=None,
                                detector_temperature=6.7)
        print(meas1)
        print(meas1.get_comments())
        if PLOTTING:
            meas1.plot()
        if SAVE_FILES:
            output_name = os.path.basename(filnam) + "_out.fits"
            print("Saving to", output_name)
            meas1.save(output_name, overwrite=True)
        del meas1
        ii += 1

    fits_file_names = [os.path.join(sim_datapath, "detector/dark_currentIM"),
                        os.path.join(sim_datapath, "detector/dark_currentLW"),
                        os.path.join(sim_datapath, "detector/dark_currentSW")]
    title = "Detector dark current"
    name = "Dark current"
    unit = "electrons"
    vname = "Temperature"
    vunit = "K"
    detectors = ['MIRIMAGE', 'MIRIFULONG', 'MIRIFUSHORT']
    ii = 0
    for filnam in fits_file_names:
        detector = detectors[ii]
        full_filename = filnam + ".fits"
        fulltitle = title + " for %s" % detector
        print("Reading", fulltitle, "from", full_filename)
        meas2 = MiriMeasurement(full_filename, title=fulltitle, name=name,
                                    unit=unit, vname=vname, vunit=vunit,
                                    interptype='LOGLIN')
#         meas2 = ascii_to_measurement(full_filename, title=fulltitle, name=name,
#                                     unit=unit, vname=vname, vunit=vunit,
#                                     interptype='LOGLIN')
        meas2.set_instrument_metadata(detector=detector, modelnam='FM',
                                filt='N/A', channel='N/A', band='N/A',
                                ccc_pos='OPEN', deck_temperature=None,
                                detector_temperature=6.7)
        print(meas2)
        print(meas2.get_comments())
        if PLOTTING:
            meas2.plot()
        if SAVE_FILES:
            output_name = os.path.basename(filnam) + "_out.fits"
            print("Saving to", output_name)
            meas2.save(output_name, overwrite=True)
        del meas2
        ii += 1

    print("Test finished.")
