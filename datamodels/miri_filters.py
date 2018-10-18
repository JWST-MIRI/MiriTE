#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

A package of data models describing the properties of spectral filters
and other devices whose transmission or efficiency varies with wavelength.

:Reference:

The STScI jwst.dataproduct documentation.

:History:

03 Jun 2010: Original Filter product created (in filters.py) (jmorin)
13 Apr 2012: Reworking of Filter product to the old MIRI data model.
22 Aug 2013: Converted to new data model.
10 Jul 2013: Modified to work with jsonschema draft 4 schemas.
             Extra code for regenerating the FITS files.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
09 Sep 2015: Corrected file names.
01 Apr 2016: Changed to YAML schemas
18 Aug 2016: Import from new JWST libraries.
31 Aug 2016: Trap and log an interpolation error. Removed obsolete code
             from table setting. Read the test tables from FITS files,
             not ASCII files. ASCII support withdrawn.
04 Jan 2017: Removed unnecessary import of DataModelPlotVisitor and
             duplication of the parent class' "plot" method.
15 Jun 2017: TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Oct 2018: The string 'ANY' is no longer recommended within CDP metadata.
             'N/A' should be used instead.


@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


# Python logging facility
import logging
logging.basicConfig(level=logging.WARNING) # Default level is warning output 
LOGGER = logging.getLogger("miri.datamodels") # Get a logger

import sys
import numpy as np
#import numpy.ma as ma
from scipy import interpolate

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Import the MIRI base data model and utilities.
from miri.datamodels.miri_model_base import MiriDataModel

# List all classes and global functions here.
__all__ = ['MiriFilter', 'MiriBandPassFilter', 'MiriQuantumEfficiency']

# Define a limit within which two floating point wavelengths are
# considered the same.
#EPS = 1.0e-12
EPS = 10 * sys.float_info.epsilon


class MiriFilter(MiriDataModel):
    """
    
    A generic data model for a MIRI filter, based on the MIRI
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
        
    filter_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (wavelength:number,
        transmission:number), giving the transmission of the filter
        as a function of wavelength.
        Or: A numpy record array containing the same information as above.
        A transmission table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    filter_name: str (optional)
        Standard name of the filter Example: 'F560W'
    filter_type: str (optional)
        Filter type. Example: 'BandPass' for band-pass filter.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_filter.schema.yaml"
    fieldnames = ('WAVELENGTH', 'TRANSMISSION')
    
    def __init__(self, init=None, filter_table=None, filter_name=None,
                 filter_type=None, **kwargs):
        """
        
        Initialises the MiriFilter class.
        
        Parameters: See class doc string.

        """
        super(MiriFilter, self).__init__(init=init, **kwargs)

        # Data type is filter.
        self.meta.model_type = 'FILTER'
        
        # Define the filter name and type, if given
        if filter_name is not None:
            self.meta.instrument.filter = filter_name
        if filter_type is not None:
            self.meta.instrument.filter_type = filter_type

        if filter_table is not None:
            try:
                self.filter_table = filter_table
            except (ValueError, TypeError) as e:
                strg = "filter_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Define the wavelength units.
#         units = self.get_data_units('filter_table')
        
        # Cached arrays
        self._wavelength = None
        self._transmission = None
        self._interptransmission = None

    def remove_negative(self, filter_table):
        """

        Replace any negative values found in the transmission measurements
        with 0.0.

        """
        new_table = []
        for record in filter_table:
            if record[1] < 0.0:
                record[1] = 0.0
            new_table.append(record)
        return new_table
       
    def get_table_arrays(self):
        """
        
        Helper function to extract the data arrays from the
        filter table.
        
        :Returns:
        
        (wavelength, transmission): tuple of 2 numpy ndarrays
        
        """
#         ftable = np.asarray(self.filter_table)
        ftable = self.filter_table
        wavelength = []
        transmission = []
        for item in ftable:
            wavelength.append(item[0])
            transmission.append(item[1])
        wavelength = np.asarray(wavelength)
        transmission = np.asarray(transmission)
        return (wavelength, transmission)

    def apply_filter(self, flux, wavelength=None, ignoreneg=False ):
        """

        Apply the filter to a given input flux.

        :Parameters:
        
        flux: array-like
            Input flux of the SED to filter.
        wavelength: array_like (optional)
            Corresponding wavelength array. If None, it is assumed that the
            filter and the input flux are defined on the same wavelength array.
        ignorenegneg: bool (optional)
            If True, negative transmission values are treated as zero when
            applying the filter.

        :Returns:
        
        filteredflux: array-like
            Flux array resulting of the filtering.

        :Note:
       
        If both wavelength arrays are identical (to within EPS), the filter
        is directly applied. Otherwise the filter's transmission function is
        interpolated using scipy.interpolate.interp1d.

        """
        # If no wavelength array provided assume wavelength == self._wavelength
        if wavelength is None: 
            wavelength = self.wavelength
            hassamewavelength = True
        else: 
            if wavelength.shape == self.wavelength.shape:
                hassamewavelength = \
                    ((wavelength-self.wavelength).max() \
                     < EPS)
            else:
                hassamewavelength = False

        if hassamewavelength:
            filteredflux = self.transmission * flux
        else:
            # Compute linear interpolation function if it doesn't exist yet.
            if self._interptransmission is None:       
                self._interptransmission = \
                    interpolate.interp1d(self.wavelength, \
                                         self.transmission)
            try:
                filteredflux = self._interptransmission(wavelength) * flux
            except ValueError as e:
                # Trap a mathematical error during the interpolation.
                strg = "apply_filter: Interpolation failure!"
                strg += " Flux will be returned unmodified by the filter."
                strg += "\n  %s" % str(e)
                LOGGER.error(strg)
                filteredflux = flux

        if ignoreneg:
            areneg = np.where(filteredflux < 0.0)
            filteredflux[areneg] = 0.0

        return filteredflux

    def get_comments(self):
        """
        
        Return a comment string describing the data object.
        This function exists for backwards compatibility.
        
        """
        strg = "Filter measurement for %s filter %s" % \
            (self.meta.instrument.filter_type, self.meta.instrument.filter)
        return strg
        
    # TODO: Is this function needed?
    def __str__(self):
        """
        
        Return the contents of the filter transmission object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()

        # Describe the filter transmission table
        strg += "\nColumn names: " + str(self.fieldnames) + "\n"
        if self.filter_table is not None:
            strg += self.get_data_str('filter_table', underline=True,
                                      underchar="-")
        return strg

    def plot_new(self, description='', visitor=None, **kwargs):
        """
        
        Plot the data product using the algorithm contained in the
        specified visitor class.

        NOTE: This plot method can be considered a quick look method.
        
        :Parameters:
        
        description: str (optional)
            Description to be added to the plot title, if required.
            This is only used if a visitor object is not provided.
        visitor: object (optional)
            A visitor object which implements the algorithm for plotting the
            contents of the data structure. If not specified, a temporary
            object will be created from the matplotlib-based visitor class
            and used.
        \*\*kwargs:
            All other keyword arguments will be passed to the visitor's visit
            method (for data arrays). These can be keywords recognised by
            the plotting package.

        """
        # Invoke the standard plot method for the parent object.
        super(MiriFilter, self).plot(description=description, visitor=visitor,
                                     **kwargs)

    def plot(self, plotaxis=None, xlabel='', ylabel='', description='',
             **kwargs):
        """

        Plot transmission as a function of wavelength.
        
        :Parameters:
        
        plotaxis: matplotlib axis object (optional)
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xlabel: str (optional)
            An X label for the plot. The default is the name of the
            wavelength table array.
        ylabel: str (optional)
            A Y label for the plot. The default is the name of
            the transmission table array.
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
        units = self.get_data_units('filter_table')
        if not xlabel:
            if units[0]:
                xlabel = "%s (%s)" % (self.fieldnames[0], units[0])
            else:
                xlabel = self.fieldnames[0]
        if not ylabel:
            if units[1]:
                ylabel = "%s (%s)" % (self.fieldnames[1], units[1])
            else:
                ylabel = self.fieldnames[1]
            
        tstrg = self.get_title(underline=False)
        if description:
            tstrg += " - " + description
        else:
            tstrg += " - %s filter %s" % (self.meta.instrument.filter_type,
                                       self.meta.instrument.filter)

        mplt.plot_xy(self.wavelength, self.transmission,
                     plotaxis=plotaxis, xlabel=xlabel, ylabel=ylabel,
                     title=tstrg, grid=True, **kwargs)
        mplt.close()

    # The wavelength and transmission arrays are accessible through
    # wavelength and transmission attributes.
    @property
    def wavelength(self):
        if self._wavelength is None:
            (wavelength,transmission) = self.get_table_arrays()
            self._wavelength = wavelength
            self._transmission = transmission
            return wavelength
        else:
            return self._wavelength

    @wavelength.setter
    def wavelength(self, data):
        raise AttributeError("wavelength attribute is read-only")

    @property
    def transmission(self):
        if self._transmission is None:
            (wavelength,transmission) = self.get_table_arrays()
            self._wavelength = wavelength
            self._transmission = transmission
            return transmission
        else:
            return self._transmission

    @transmission.setter
    def transmission(self, data):
        raise AttributeError("transmission attribute is read-only")


class MiriBandPassFilter(MiriFilter):
    """
    
    A data model for a MIRI band pass filter.
    
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
        
    filter_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (wavelength:number,
        transmission:number), giving the transmission of the filter
        as a function of wavelength.
        Or: A numpy record array containing the same information as above.
        A transmission table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    filter_name: str (optional)
        Standard name of the filter Example: 'F560W'
    wavcentre: float (optional)
        Central wavelength of the filter.
    fwhm: float (optional)
        Full-width at half maximum of the bandpass.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
        
    """
    
    def __init__(self, init=None, filter_table=None, filter_name=None,
                 wavecent=None, fwhm=None, **kwargs):
        """
        
        Initialises the MiriQuantumEfficiency class.
        
        Parameters: See class doc string.

        """
        super(MiriBandPassFilter, self).__init__(init=init,
                                             filter_table=filter_table,
                                             filter_name=filter_name,
                                             filter_type='BandPass',
                                             **kwargs)

        # Define the central wavelength and fwhm, if defined.
        if wavecent:
            self.meta.instrument.filter_wavecent = wavecent
        if fwhm:
            self.meta.instrument.filter_fwhm = fwhm

    def remove_out_of_band(self, filter_table, nfwhm=5, wavecentre=None,
                           fwhm=None):
        """

        Remove out-of-bandpass wiggles, transmission outside the domain
        [wavcentre-removeout*fwhm; wavcentre+removeout*fwhm] is set to zero.

        """
        if wavecentre is None:
            wavecentre = self.meta.instrument.filter_wavecent
        if fwhm is None:
            fwhm = self.meta.instrument.filter_fwhm = fwhm
        new_table = []
        lower_edge = wavecentre - fwhm * nfwhm
        upper_edge = wavecentre + fwhm * nfwhm
        for record in filter_table:
            if record[0] < lower_edge:
                record[1] = 0.0
            if record[0] > upper_edge:
                record[1] = 0.0
            new_table.append(record)
        return new_table


class MiriQuantumEfficiency(MiriFilter):
    """
    
    A data model describing the quantum efficiency of a MIRI detector.
    This uses an almost identical data model to MiriFilter, except the
    values stored are an efficiency rather than a transmission, and the
    table is identified by detector name rather than filter name.
    
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
        
    qe_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (wavelength:number,
        efficiency:number), giving the quantum efficiency of the
        detector as a function of wavelength.
        Or: A numpy record array containing the same information as above.
        A QE table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    detector: str (optional)
        Name of the detector to which this QE measurement belongs
    temperature: float (optional)
        Detector temperature at the time of the QE measurement.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
        
    """
    
    def __init__(self, init=None, qe_table=None, detector=None,
                 temperature=None, **kwargs):
        """
        
        Initialises the MiriQuantumEfficiency class.
        
        Parameters: See class doc string.

        """
        super(MiriQuantumEfficiency, self).__init__(init=init,
                                                    filter_table=qe_table,
                                                    filter_name='N/A',
                                                    filter_type='QE', **kwargs)

        # Data type is QE.
        self.meta.model_type = 'QE'
        
        # Define the detector name and detector temperature, if given
        if detector is not None:
            self.meta.instrument.detector = detector
        if temperature is not None:
            self.meta.instrument.detector_temperature = temperature

    def get_comments(self):
        """
        
        Return a comment string describing the data object.
        This function exists for backwards compatibility.
        
        """
        strg = "QE measurement for detector %s" % self.meta.instrument.detector
        return strg

    @property
    def efficiency(self):
        # efficiency is a synonym for filter transmission
        return self.transmission

    @efficiency.setter
    def efficiency(self, data):
        # efficiency is a synonym for filter transmission
        self.transmission = data

    def plot(self, plotaxis=None, xlabel='', ylabel='', description='', **kwargs):
        """
 
        Plot quantum efficiency as a function of wavelength.
         
        :Parameters:
         
        plotaxis: matplotlib axis object (optional)
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xlabel: str (optional)
            An X label for the plot. The default is the name of the
            wavelength table array.
        ylabel: str (optional)
            A Y label for the plot. The default is 'Quantum Efficiency'.
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
        if not ylabel:
            ylabel = 'Quantum EFficiency'
        if not description:
            description = "detector %s" % self.meta.instrument.detector
        super(MiriQuantumEfficiency, self).plot(plotaxis=plotaxis, xlabel=xlabel,
                                                ylabel=ylabel, description=description,
                                                **kwargs)

def ascii_to_filter(filename, filter_name=None, detector=None, temperature=None, 
                    filter_type=None, wcol=0, tcol=None, **kwargs):
    """
        
    Create a MiriFilter data product from an ASCII file.
                
    :Parameters:
    
    filename: str
        The name of an ASCII file from which to create the Filter Data Product.
    filter_name: str (optional)
        Standard name of a filter Example: 'F560W'
    detector: str (optional)
        Name of the detector to which a QE measurement belongs
    temperature: float (optional)
        Detector temperature at the time of a QE measurement.    
    filter_type: str (optional)
        Filter type. Example: 'BandPass' for band-pass filter.
    wcol: int (optional)
        The column number containing the wavelength.
        Defaults to 0.
    tcol: int (optional)
        The column number containing the transmission or efficiency.
        Defaults to 1 for a QE measurement or 2 for any other
        kind of filter.

    """
    strg = "Reading a MiriFilter model from an ASCII file "
    strg += "is not longer supported."
    raise NotImplementedError(strg)
#     # This function no longer works
#     if filter_type is not None and filter_type == 'QE':
#         if tcol is None:
#             tcol = 1
#     else:
#         if tcol is None:
#             tcol = 2
#     
#     try:
#         data = np.loadtxt(filename)
#     except Exception as e:
#         # If the file could not be read re-raise the exception
#         # with a more meaningful error message.
#         strg = "%s: Could not read Filter Data Product ASCII file.\n   %s" % \
#             (e.__class__.__name__, e)
#         raise IOError(strg)
# 
#     wavelength = data[:, wcol]
#     transmission = data[:, tcol]
# 
#     filter_table = []
#     for wval, tval in zip(wavelength, transmission):
#         filter_table.append([wval, tval])
# 
#     # Initialise the data product from the table just created.
#     if filter_type == 'QE':
#         miri_filter = MiriQuantumEfficiency( qe_table=filter_table,
#                                 detector=detector, temperature=temperature,
#                                 **kwargs)
#     else:
#         miri_filter = MiriFilter( filter_table=filter_table,
#                                   filter_name=filter_name,
#                                   filter_type=filter_type, **kwargs)
#     
#     return miri_filter

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriFilter module.")
    
    import os
    import miri.simulators.data as simdata
    sim_datapath = simdata.__path__[0]

    PLOTTING = False       # Set to False to turn off plotting.
    SAVE_FILES = False     # Set to False to avoid saving files.

    print("\nTesting MiriFilter class." + \
          "\n=========================")
    print("Filter with transmissions derived from list of tuples:")
    transmissions = [
                    (5.80, 0.8698470),
                    (5.81, 0.8759494),
                    (5.82, 0.8944225),
                    (5.83, 0.8899569),
                    (5.84, 0.8760563),
                    (5.85, 0.8726164),
                    (5.86, 0.8782486),
                    (5.87, 0.8753881),
                    (5.88, 0.8844002),
                    (5.89, 0.8682995),
                    (5.90, 0.8495247),
                    (5.91, 0.8289118),
                    (5.92, 0.8211463),
                    (5.93, 0.8199366),
                    (5.94, 0.8202344),
                    (5.95, 0.7952070),
                    (5.96, 0.7884885),
                    (5.97, 0.7938501),
                    (5.98, 0.7938051),
                    (5.99, 0.8033671),
                    (6.00, 0.7985086)
                    ]
 
    with MiriFilter( filter_table=transmissions, filter_name='N/A',
                     filter_type='LowPass' ) as testfilter1:
        print(testfilter1)
        if PLOTTING:
            testfilter1.plot(description="testfilter1")
            testfilter1.plot_new(description="testfilter1")
        if SAVE_FILES:
            testfilter1.save("test_filter_model1.fits", overwrite=True)
        del testfilter1
    del transmissions

# NOTE: This code (with SAVE_FILES=True) recreates the example_filter.fits file.
#     print("\nExample filter read from ASCII file")
#     import miri.datamodels.data
#     datapath = miri.datamodels.data.__path__[0]
#     test_file_name = os.path.join(datapath, "example_filter")
#     example = ascii_to_filter(test_file_name + ".txt", filter_name='F560W',
#                               filter_type='BandPass')
#     example.meta.instrument.filter_wavecent = 5.6
#     example.meta.instrument.filter_fwhm = 1.2
#     if PLOTTING:
#         example.plot()
#         example.plot_new()
#     if SAVE_FILES:
#         example.save("example_filter.fits", overwrite=True)
#     del example
    
# -----------------------------------------------------------------------

    # NOTE: The following tests assume the miri.simulators package has been
    # built and following list of FITS files are available.
    print("\nFilters with transmissions read from FITS files:")
    fits_file_names = [os.path.join(sim_datapath, "filters/IM-01_F560W"),
                       os.path.join(sim_datapath, "filters/IM-02_F770W"),
                       os.path.join(sim_datapath, "filters/IM-03_F1000W"),
                       os.path.join(sim_datapath, "filters/IM-04_F1130W"),
                       os.path.join(sim_datapath, "filters/IM-05_F1280W"),
                       os.path.join(sim_datapath, "filters/IM-06_F1500W"),
                       os.path.join(sim_datapath, "filters/IM-07_F1800W"),
                       os.path.join(sim_datapath, "filters/IM-08_F2100W"),
                       os.path.join(sim_datapath, "filters/IM-09_F2550W"),
                       os.path.join(sim_datapath, "filters/IM-10_F2550WR")]
    wavecents = [5.6, 7.7, 10.0, 11.3, 12.8, 15.0, 18.0, 21.0, 25.0, 25.0]
    fwhms =     [1.2, 2.2,  2.0,  2.2,  2.4,  3.0,  3.0,  5.0,  5.0,  5.0]

    ii = 0
    for filnam in fits_file_names:
        full_filename = filnam + ".fits"
        filter_name = filnam.split('_')[-1]
        print("Reading filter", filter_name, "from", full_filename)
        flt2 = MiriBandPassFilter(full_filename, filter_name=filter_name)
#         flt2 = ascii_to_filter(full_filename, filter_name=filter_name,
#                                filter_type='BandPass')
        flt2.meta.instrument.filter_wavecent = wavecents[ii]
        flt2.meta.instrument.filter_fwhm = fwhms[ii]
        print(flt2)
        print(flt2.get_comments())
        if PLOTTING:
            flt2.plot()
            flt2.plot_new()
        if SAVE_FILES:
            output_name = os.path.basename(filnam) + "_out.fits"
            print("Saving to", output_name)
            flt2.save(output_name, overwrite=True)
        del flt2
        
        if SAVE_FILES:
            with MiriFilter(init=output_name) as new_model:
                print("Read back:", new_model)
            del new_model
        ii += 1
    del fits_file_names, wavecents, fwhms

# -----------------------------------------------------------------------

    print("\nTesting MiriQuantumEfficiency class." + \
          "\n====================================")
    
    print("QE measurements read from FITS files:")
    # The following code (with SAVE_FILES=True) can be used to recreate the
    # FITS versions of these simulator files.
    fits_file_names = [os.path.join(sim_datapath, "detector/qe_measurementIM"),
                       os.path.join(sim_datapath, "detector/qe_measurementLW"),
                       os.path.join(sim_datapath, "detector/qe_measurementSW")]
    detectors = ['MIRIMAGE', 'MIRIFULONG', 'MIRIFUSHORT']

    ii = 0
    for filnam in fits_file_names:
        detector = detectors[ii]
        full_filename = filnam + ".fits"
        print("Reading QE for", detector, "from", full_filename)
        flt3 = MiriQuantumEfficiency( full_filename, detector=detector,
                               temperature=6.7)
#         flt3 = ascii_to_filter(full_filename, detector=detector,
#                                temperature=6.7, filter_type='QE')
        print(flt3)
        print(flt3.get_comments())
        if PLOTTING:
            flt3.plot()
            flt3.plot_new()
        if SAVE_FILES:
            output_name = os.path.basename(filnam) + "_out.fits"
            print("Saving to", output_name)
            flt3.save(output_name, overwrite=True)
        del flt3
          
        if SAVE_FILES:
            with MiriQuantumEfficiency(init=output_name) as new_model:
                print("Read back:", new_model)
            del new_model
        ii += 1
    del fits_file_names, detectors

    print("Tests finished.")
