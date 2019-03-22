#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

A package of data models describing the photon conversion efficiency
of the MIRI optics, which includes the reflectance of the optics, the
transmission of the spectral filters and the QE of the detectors as a
function of wavelength.

:Reference:

The STScI jwst.datamodels documentation.
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

22 Aug 2013: Original version as filter model, MiriFilter.
22 Jun 2016: Converted to new PCE model, MiriPceModel. This model has 3 table
             columns instead of 2, but the overall API is the same.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
             Do not set observation or target metadata. Neither are
             appropriate for a reference file.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema. Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).

@author: Steven Beard (UKATC)

"""

import sys
import numpy as np
#import numpy.ma as ma
from scipy import interpolate

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.plotting import DataModelPlotVisitor

# List all classes and global functions here.
__all__ = ['MiriPceModel', 'ascii_to_pce']

# Define a limit within which two floating point wavelengths are
# considered the same.
EPS = 10 * sys.float_info.epsilon


def ascii_to_pce(filename, component=None, detector=None, ncols=3,
                 clip_values=True, **kwargs):
    """
        
    Create a MiriPceModel data product from an ASCII file.
                
    :Parameters:
    
    filename: str
        The name of an ASCII file from which to create the Filter Data Product.
    component: str (optional)
        Standard ETC name of a component Example: 'F560W'
    detector: str (optional)
        Name of the detector to which the QE part of a PCE measurement belongs.
    ncols: int (optional)
        Number of columns to read from file. Defaults to 3.
        Can be set to 2 to ignore a conversion column and assume conversion
        is always 1.0.
    clip_values: boolean (optional)
        Set to True (the default) to clip efficiencies to the range 0-1.

    :Returns:
    
    miri_pce: MiriPceModel data model.

    """
    # Read the table columns from the specified ASCII file into a numpy array.
    try:
        data = np.loadtxt(filename)
    except Exception as e:
        # If the file could not be read re-raise the exception
        # with a more meaningful error message.
        strg = "%s: Could not read PCE data product ASCII file.\n   %s" % \
            (e.__class__.__name__, e)
        raise IOError(strg)

    wavelength = data[:, 0]
    efficiency = data[:, 1]
    if ncols > 2:
        conversion = data[:, 2]
    else:
        conversion = np.ones_like( wavelength )
        
    if clip_values:
        efficiency = np.clip(efficiency, 0.0, 1.0)

    # Construct a record structure from the columns.
    pce_table = []
    for wval, eval, cval in zip(wavelength, efficiency, conversion):
        pce_table.append((wval, eval, cval))

    # Initialise the data product from the table just created.
    miri_pce = MiriPceModel( pce_table=pce_table, component=component,
                        detector=detector, **kwargs)
    
    return miri_pce


class MiriPceModel(MiriDataModel):
    """
    
    A generic data model for a MIRI photon conversion efficiency,
    based on the MIRI base model, MiriDataModel.
    
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
        
    pce_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (wavelength:number,
        efficiency:number, conversion:number), giving the efficiency
        of the component as a function of wavelength.
        Or: A numpy record array containing the same information as above.
        A efficiency table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    component: str (optional)
        Name of instrument component, as regonised by the ETC.
        Examples: 'F560W' or 'CHANNEL3SHORT'.
        If not given, a name will be constructed from a combination of
        filter name, channel and band.
    detector: str (optional)
        The name of the detector associated with this PCE data.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_pce.schema.yaml"
    fieldnames = ('WAVELENGTH', 'EFFICIENCY', 'CONVERSION')
    
    def __init__(self, init=None, pce_table=None, component=None,
                 detector=None, **kwargs):
        """
        
        Initialises the MiriPceModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPceModel, self).__init__(init=init, **kwargs)

        # Data type is PCE.
        self.meta.reftype = 'PCE'
        model_type = get_my_model_type( self.__class__.__name__ )
        if model_type is not None:
            self.meta.model_type = model_type        

        # This is a reference data model.
        self._reference_model()
        
        # Define the component name, if given
        if component is not None:
            self.meta.etc.component = component
        # If the component name is not defined, define it from the other metadata.
        if self.meta.etc.component is None:
            if self.meta.instrument.filter is not None:
                self.meta.etc.component = self.meta.instrument.filter
            elif self.meta.instrument.channel is not None and \
                 self.meta.instrument.band is not None:
                self.meta.etc.component = "CHANNEL" + \
                    str(self.meta.instrument.channel) + \
                    str(self.meta.instrument.band)
                    
        # Define the detector identifier, if specified.
        if detector is not None:
            self.meta.instrument.detector = detector

        if pce_table is not None:
            try:
                #pce_table = np.recarray(phot_table)
                self.pce_table = pce_table
            except (ValueError, TypeError) as e:
                strg = "pce_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # Copy the table column units from the schema, if defined.
        pce_units = self.set_table_units('pce_table')
                
        # Cached arrays
        self._wavelength = None
        self._efficiency = None
        self._interpefficiency = None

    def remove_negative(self, pce_table):
        """

        Replace any negative values found in the efficiency measurements
        with 0.0.
        
        :Parameters:
        
        pce_table: list of tuples or numpy record array
            Input table
            
        :Returns:
        
        new_table: list of tuples or numpy record array
            Output table

        """
        new_table = []
        for record in pce_table:
            if record[1] < 0.0:
                record[1] = 0.0
            new_table.append(record)
        return new_table
       
    def get_table_arrays(self):
        """
        
        Helper function to extract the data arrays from the
        PCE table.
        
        :Returns:
        
        (wavelength, efficiency, conversion): tuple of 3 numpy ndarrays
        
        """
        ftable = np.asarray(self.pce_table)
        wavelength = []
        efficiency = []
        conversion = []
        for item in ftable:
            wavelength.append(item[0])
            efficiency.append(item[1])
            conversion.append(item[2])
        wavelength = np.asarray(wavelength)
        efficiency = np.asarray(efficiency)
        conversion = np.asarray(conversion)
        return (wavelength, efficiency, conversion)

    def apply_filter(self, flux, wavelength=None, ignoreneg=False ):
        """

        Use the efficiency table to filter a given input flux.

        :Parameters:
        
        flux: array-like
            Input flux of the SED to filter.
        wavelength: array_like (optional)
            Corresponding wavelength array. If None, it is assumed that the
            filter and the input flux are defined on the same wavelength array.
        ignorenegneg: bool (optional)
            If True, negative efficiency values are treated as zero when
            applying the filter.

        :Returns:
        
        filteredflux: array-like
            Flux array resulting of the filtering.

        :Note:
       
        If both wavelength arrays are identical (to within EPS), the filter
        is directly applied. Otherwise the filter's efficiency function is
        interpolated using scipy.interpolate.interp1d.

        """
        # If no wavelength array provided assume wavelength == self.wavelength
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
            filteredflux = self.efficiency * flux
        else:
            # Compute linear interpolation function if it doesn't exist yet.
            if self._interpefficiency is None:       
                self._interpefficiency = \
                    interpolate.interp1d(self.wavelength, \
                                         self.efficiency)
            filteredflux = self._interpefficiency(wavelength) * flux

        if ignoreneg:
            areneg = np.where(filteredflux < 0.0)
            filteredflux[areneg] = 0.0

        return filteredflux

    def plot(self, plotaxis=None, xlabel='', ylabel='', description='',
             **kwargs):
        """
 
        Plot efficiency as a function of wavelength.
         
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
            the efficiency table array.
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
        units = self.get_data_units('pce_table').split(", ")
        if not xlabel:
            units0 = units[0]
            units0 = units0.lstrip('(')
            units0 = units0.lstrip('"')
            units0 = units0.rstrip(')')
            units0 = units0.rstrip('"')
            if units0:
                xlabel = "%s (%s)" % (self.fieldnames[0], units0)
            else:
                xlabel = self.fieldnames[0]
        if not ylabel:
            units1 = units[1]
            units1 = units1.lstrip('(')
            units1 = units1.lstrip('"')
            units1 = units1.rstrip(')')
            units1 = units1.rstrip('"')
            if units1:
                ylabel = "%s (%s)" % (self.fieldnames[1], units1)
            else:
                ylabel = self.fieldnames[1]
             
        tstrg = self.get_title(underline=False)
        if description:
            tstrg += " - " + description
        else:
            tstrg += " - %s" % self.meta.etc.component
 
        mplt.plot_xy(self.wavelength, self.efficiency,
                     plotaxis=plotaxis, xlabel=xlabel, ylabel=ylabel,
                     title=tstrg, grid=True, **kwargs)
        mplt.close()

    # The wavelength and efficiency arrays are accessible through
    # wavelength and efficiency attributes.
    @property
    def wavelength(self):
        if self._wavelength is None:
            (wavelength, efficiency, conversion) = self.get_table_arrays()
            self._wavelength = wavelength
            self._efficiency = efficiency
            self._conversion = conversion
            return wavelength
        else:
            return self._wavelength

    @wavelength.setter
    def wavelength(self, data):
        raise AttributeError("wavelength attribute is read-only")

    @property
    def efficiency(self):
        if self._efficiency is None:
            (wavelength,efficiency, conversion) = self.get_table_arrays()
            self._wavelength = wavelength
            self._efficiency = efficiency
            self._conversion = conversion
            return efficiency
        else:
            return self._efficiency

    @efficiency.setter
    def efficiency(self, data):
        raise AttributeError("efficiency attribute is read-only")

    @property
    def conversion(self):
        if self._efficiency is None:
            (wavelength,efficiency, conversion) = self.get_table_arrays()
            self._wavelength = wavelength
            self._efficiency = efficiency
            self._conversion = conversion
            return conversion
        else:
            return self._conversion

    @conversion.setter
    def conversion(self, data):
        raise AttributeError("conversion attribute is read-only")



#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriPceModel module.")
    
    import os
    import miri.simulators.data as simdata
    sim_datapath = simdata.__path__[0]

    PLOTTING = False       # Set to False to turn off plotting.
    SAVE_FILES = False     # Set to False to avoid saving files.

    print("\nTesting MiriPceModel class." + \
          "\n=========================")
    print("PCE data model created from alist of tuples:")
    efficiencies = [
                    (5.80, 0.8698470, 1.0),
                    (5.81, 0.8759494, 1.0),
                    (5.82, 0.8944225, 1.0),
                    (5.83, 0.8899569, 1.0),
                    (5.84, 0.8760563, 1.0),
                    (5.85, 0.8726164, 1.0),
                    (5.86, 0.8782486, 1.0),
                    (5.87, 0.8753881, 1.0),
                    (5.88, 0.8844002, 1.0),
                    (5.89, 0.8682995, 1.0),
                    (5.90, 0.8495247, 1.0),
                    (5.91, 0.8289118, 1.0),
                    (5.92, 0.8211463, 1.0),
                    (5.93, 0.8199366, 1.0),
                    (5.94, 0.8202344, 1.0),
                    (5.95, 0.7952070, 1.0),
                    (5.96, 0.7884885, 1.0),
                    (5.97, 0.7938501, 1.0),
                    (5.98, 0.7938051, 1.0),
                    (5.99, 0.8033671, 1.0),
                    (6.00, 0.7985086, 1.0)
                    ]
 
    with MiriPceModel( pce_table=efficiencies, component='F560W',
                  detector='MIRIMAGE' ) as testpce1:
        testpce1.set_instrument_metadata(detector='MIRIMAGE', filt='F560W',
                                         ccc_pos='OPEN', deck_temperature=11.0,
                                         detector_temperature=6.0)
        testpce1.set_subarray_metadata('GENERIC')
        testpce1.set_housekeeping_metadata('UK', author='MIRI team',
                                           version='1.0', date='TODAY',
                                           useafter='',
                                           description='PCE test data')
        testpce1.set_exposure_type()
        print(testpce1)
        if PLOTTING:
            testpce1.plot(description="test_pce_model1")
        if SAVE_FILES:
            testpce1.save("test_pce_model1.fits", overwrite=True)
            print("Saved to test_pce_model1.fits")
        del testpce1
    del efficiencies

    print("\nPCE data model created from an ASCII file")
    import miri.datamodels.data
    datapath = miri.datamodels.data.__path__[0]
    test_file_name = os.path.join(datapath, "example_filter")
    example1 = ascii_to_pce(test_file_name + ".txt", component='F560W',
                           detector='MIRIMAGE', ncols=2)
    print(example1)
    if PLOTTING:
        example1.plot(description="test_pce_ascii1")
    if SAVE_FILES:
        example1.save("test_pce_ascii1.fits", overwrite=True)
        print("Saved to test_pce_ascii1.fits")
    del example1

    test_file_name = os.path.join(datapath, "example_qe")
    example2 = ascii_to_pce(test_file_name + ".txt", component='MIRIMAGE',
                           detector='MIRIMAGE', ncols=2)
    print(example2)
    if PLOTTING:
        example2.plot(description="test_pce_ascii2")
    if SAVE_FILES:
        example2.save("test_pce_ascii2.fits", overwrite=True)
        print("Saved to test_pce_ascii2.fits")
    del example2

    print("Tests finished.")
