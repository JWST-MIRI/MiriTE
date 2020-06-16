#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module spectraldata - Contains data models common to all MIRI spectroscopy
packages.

:History:

03 Jul 2013: Created with Spectrum1D, Spectrum2D and Spectrum3D
19 Aug 2013: Corrected typo in the initialisation of alpha, beta
             and wavelength arrays from starting values and steps.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
20 May 2014: Modified for jsonschema draft 4. Define data units using
             the set_data_units method.
09 Dec 2014: Obsolete field_def table replaced by dq_def.
19 Aug 2015: Removed MiriImageModel.
21 Sep 2016: Try both methods for importing flattening_lib.
28 Sep 2016: Changed to new data models and YAML schemas.
29 Sep 2016: Removed all code relating to old data models.
30 Sep 2016: Corrected a typo in the YAML schema for Spectrum3D.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
17 Oct 2019: Moved from MiriTeam/spectroscopy to MiriTE/datamodels
             to work around a problem with AsdfExtension.

@author: Steven Beard (UKATC), Ruyman Azzollini (DIAS)

"""


import numpy as np
import numpy.ma as ma

#from pdb import set_trace as stop

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

from miri.datamodels.miri_measured_model import MiriMeasuredModel
from miri.tools.spec_tools import spec2d_to_spec1d, \
                                  spec3d_to_image, \
                                  spec3d_to_spec1d


class Spectrum1D(MiriMeasuredModel):
    """
    
    A data model for 1-D MIRI spectra.
    
    :Parameters:
    
    These parameters are passed to MiriMeasuredModel class.

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
        Either: A list of tuples containing (value:int, name:str, comment:str),
        giving the meaning of values stored in the data quality array. For
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to a trivial good/bad mask.

    These parameters are passed to the Spectrum1D class.

    wavelength: 1-D numpy array, optional
        An array containing the wavelength axis.
        If not given an array will be created full of pixel indices.
    wstart: float, optional
        Starting wavelength. An alternative to the wavelength array.
    wstep: float, optional
        Wavelength increment. An alternative to the wavelength array.
    wunit, str, optional
        Wavelength units (if different from the default specified in
        the schema). Examples: "microns", "angstroms", "nanometres".
    title, str, optional
        A title for the data object (if different from the default
        specified in the schema). Examples: "Spectrum of NGC1365".
    label, str, optional
        A label for the science data (if different from the default
        specified in the schema). Examples: "Intensity", "Photon count",
        "Calibrated Flux".
    unit, str, optional
        Science data units (if different from the default specified in
        the schema). Examples: "Jy", "Photons".
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst_lib documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_spectrum1D.schema"
    
    # Small errors can't weight a pixel by more than the average weighting
    # multiplied by this amount.
    MAX_WEIGHT_GAIN = 10000.0
    
    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 wavelength=None, wstart=0.0, wstep=1.0, wunit='', title='',
                 label='', unit='', **kwargs):
        super(Spectrum1D, self).__init__(init=init, data=data, dq=dq, err=err,
                                         dq_def=dq_def, title=title,
                                         **kwargs)
        # Data type is MIRI spectrum.
        self.meta.datatype = 'Spectrum'

        # The title is defined by the MiriDataModel class.
        # Set the science data label, if explicitly specified.
        # Otherwise leave it as defined in the schema.
        # COMMENTED OUT. DOES NOT WORK WITH NEW DATA MODEL SCHEMAS.
#         if label and ('title' in self.schema['properties']['data']):
#             self.schema['properties']['data']['title'] = label

        # Create an array of wavelengths. If not given explicitly,
        # the default is a linear wavelength scale. The wavelength array
        # must match the dimensions of the data array.
        # TODO: This array will eventually be obtained from the embedded WCS information.
        if wavelength is not None:
            self.wavelength = wavelength
        elif self.data is not None:
            if (not self._isvalid(self.wavelength)) or \
                (np.count_nonzero( self.wavelength ) == 0):
                # The wavelength array is empty or full of zeros (the default).
                # Fill the array using the start and increment supplied.
                wstop = wstart + (self.data.shape[0]-1) * wstep
                wavelength = np.linspace(wstart, wstop, self.data.shape[0])
                self.wavelength = wavelength
                
        # Check that the data and wavelength arrays are compatible.
        # (The schema already ensures the wavelength array is 1-D)
        if self.data is not None and self.wavelength is not None:
            if self.wavelength.size != self.data.shape[0]:
                strg = "Wavelength array has wrong size "
                strg += "(%d instead of %d)." % \
                    (self.wavelength.size, self.data.shape[0])
                raise ValueError(strg)
        
        # If the units of the wavelength array are explicitly specified,
        # copy the string to the metadata and the schema. Otherwise obtain
        # the default units from the schema.
        wavunits = self.set_data_units('wavelength', units=wunit)

    # In a Spectrum object the attribute "spectrum" is a synonym for the SCI array.
    @property
    def spectrum(self):
        return self.data

    @spectrum.setter
    def spectrum(self, data):
        self.data = data

    def _get_plot_labels(self, description):
        """
        
        Helper function to put together the spectrum plot labels
        
        """
        tstrg = self.get_title()
        if description:
            tstrg += " - %s" % description
            
        xlabel = "Wavelength"
        if self.meta.wavelength.units:
            xlabel += " (" + self.meta.wavelength.units + ")"

        ylabel = self.get_data_title('data')
        if self.meta.data.units:
            ylabel += " (" + self.meta.data.units + ")"
            
        return tstrg, xlabel, ylabel

    def plot_spectrum(self, plotaxis=None, errorbars=True, description=''):
        """
        
        Plot the spectrum as a 1-D plot (with optional error bars)
        
        """
        # Define fixed labels
        tstrg, xlabel, ylabel = self._get_plot_labels(description)

        # Plot the spectrum as an XY plot. If an ERR array exists it
        # will be plotted with error bars.
        # If an array is masked, the gaps are filled with
        # the fill value before the array is plotted.
        if errorbars and self.has_dataarray('ERR'):
            mplt.plot_xy(self.wavelength, self.data_filled,
                         yerr=self.err_filled, plotaxis=plotaxis,
                         linefmt='bo', xlabel=xlabel, ylabel=ylabel,
                         title=tstrg)
        else:
            mplt.plot_xy(self.wavelength, self.data_filled,
                         yerr=None, plotaxis=plotaxis,
                         linefmt='bo', xlabel=xlabel, ylabel=ylabel,
                         title=tstrg)


class Spectrum2D(Spectrum1D):
    """
    
    A data model for 2-D MIRI spectra.
    
    NOTE: In this data model the wavelength is the Y axis.
    
    :Parameters:
    
    The same as Spectrum1D, plus
    
    slit: 1-D numpy array, optional
        An array containing the along slit axis.
        If not given an array will be created full of pixel indices.
    sstart: float, optional
        Starting slit coordinate. An alternative to the slit array.
    sstep: float, optional
        Slit coordinate increment. An alternative to the slit array.
    sunit, str, optional
        Slit coordinate units (if different from the default specified in
        the schema).
    
    """
    schema_url = "miri_spectrum2D.schema"
    
    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 wavelength=None, wstart=0.0, wstep=1.0, wunit='', 
                 slit=None, sstart=0.0, sstep=1.0, sunit='',
                 title='', label='', unit='', **kwargs):
        super(Spectrum2D, self).__init__(init=init, data=data, dq=dq, err=err,
                                         dq_def=dq_def,
                                         wavelength=wavelength, wstart=wstart,
                                         wstep=wstep, wunit=wunit, title=title,
                                         label=label, unit=unit, **kwargs)
        # Define which of the axes represents wavelength.
        # TODO: This information should come from the schema.
        self.waveaxis = 0
        
        # The wavelength array has already been defined in the
        # above constructor. For 2-D data, the slit array needs to be added.
        # The slit array must match the dimensions of the data array.
        # TODO: This array will eventually be obtained from the embedded WCS information.
        if slit is not None:
            self.slit = slit
        elif self.data is not None:
            if (not self._isvalid(self.slit)) or \
                (np.count_nonzero( self.slit ) == 0):
                # The slit array is empty or full of zeros (the default).
                # Fill the array using the start and increment supplied.
                sstop = sstart + (self.data.shape[1]-1) * sstep
                slit = np.linspace(sstart, sstop, self.data.shape[1])
                self.slit = slit

        # Check that the data and slit arrays are compatible.
        # (The schema already ensures the slit array is 1-D)
        if self.data is not None and self.slit is not None:
            if self.slit.size != self.data.shape[1]:
                strg = "Slit array has wrong size "
                strg += "(%d instead of %d)." % \
                    (self.slit.size, self.data.shape[1])
                raise ValueError(strg)

        # If the units of the slit array are explicitly specified,
        # copy the string to the metadata and the schema. Otherwise obtain
        # the default units from the schema.
        slitunits = self.set_data_units('slit', units=sunit)

    def flatten_to_spectrum(self, weighted=True, perfect=False, masked=False):
        """
        
        Flatten a 2-D spectrum in the spatial direction to generate one
        combined spectrum.
        
        The output signal is the (weighted) average of the input signal.
        The output error is the RMS of the input error.
    
        :Parameters:
        
        weighted: bool, optional
            Set to True to use the error array (if present) to make a
            weighted average.
            Set to False for an unweighted average.
        perfect: bool, optional
            Set to True to regard as bad any pixel in the spectrum
            which has contributing bad pixel.
            Set to False to regard a pixel in the spectrum as bad only
            if all of the contributing pixels are bad (and use the
            error array to mark the reduction in quality).
        masked: bool, optional
            Set to True to use the masked versions of the data arrays.
            By default, the non-masked versions are used. This is more
            efficient but may generate runtime warning messages.

        :Returns:

        spectrum: Spectrum1D
            A 1-D spectrum object containing the extracted spectrum.        

        """
        if masked:
            (flat_data, flat_error, flat_quality) = \
                spec2d_to_spec1d(self.data_masked,
                                 self.err_masked, self.dq,
                                 waveaxis=self.waveaxis,
                                 weighted=weighted, perfect=perfect)
        else:
            (flat_data, flat_error, flat_quality) = \
                spec2d_to_spec1d(self.data, self.err, self.dq,
                                 waveaxis=self.waveaxis,
                                 weighted=weighted, perfect=perfect)
        title = "(flattened to 1-D spectrum)"

        # NOTE: The STScI data model doesn't not support masked arrays,
        # so they need to be converted back to contiguous arrays before
        # they are used to create a new data product.
        spectrum = Spectrum1D(data=np.ascontiguousarray(flat_data),
                              err=np.ascontiguousarray(flat_error),
                              dq=np.ascontiguousarray(flat_quality),
                              wavelength=self.wavelength, title=title)
        return spectrum

    def plot_spectrum(self, plotaxis=None, errorbars=True, description=''):
        """
        
        Plot the spectrum as a series of overlaid 1-D plots
        (with optional error bars)
        
        """
        # Define fixed labels
        tstrg, xlabel, ylabel = self._get_plot_labels(description)
            
        # Create a matplotlib figure and axis.
        tstrg += " - slit positions 1 to %d" % self.data.shape[1]
        fig = mplt.new_figure(1, stitle=tstrg)
        ax = mplt.add_subplot(fig, 1, 1, 1)

        # The plots will cycle around these defined line formats,
        linefmtlist = ['bo', 'rx', 'g+', 'ko', \
                       'bx', 'r+', 'go', 'kx', \
                       'b+', 'ro', 'gx', 'k+' ]
        
        # Cycle through the list of spectra, remembering that the
        # wavelength is the [0] axis
        # TODO: Can some of this duplicated code be packed up?
        for slc in range(0, self.data.shape[1]):
            spectrum = self.data_filled[:, slc]

            ifmt = slc % len(linefmtlist)
            lfmt = linefmtlist[ifmt] 

            # Plot each spectrum as an XY plot overlaid on the same axis.
            if errorbars and self.has_dataarray('ERR'):
                errors = self.err_filled[:, slc]
                mplt.plot_xy(self.wavelength, spectrum, yerr=errors,
                             plotfig=fig, plotaxis=ax, linefmt=lfmt,
                             xlabel=xlabel, ylabel=ylabel, title='')
            else:
                mplt.plot_xy(self.wavelength, spectrum, plotfig=fig,
                             plotaxis=ax, linefmt=lfmt, xlabel=xlabel,
                             ylabel=ylabel, title='')
        mplt.show_plot()


class Spectrum3D(Spectrum1D):
    """
    
    A data model for 3-D MIRI spectra.

    NOTE: In this data model the wavelength is the Z axis.
    
    :Parameters:
    
    The same as Spectrum1D, plus
    
    alpha: 1-D numpy array, optional
        An array containing the alpha axis.
        If not given an array will be created full of pixel indices.
    astart: float, optional
        Starting alpa coordinate. An alternative to the alpha array.
    astep: float, optional
        Alpha coordinate increment. An alternative to the alpha array.
    aunit, str, optional
        Alpha coordinate units (if different from the default specified in
        the schema).
    beta: 1-D numpy array, optional
        An array containing the beta (across slit) axis.
        If not given an array will be created full of pixel indices.
    bstart: float, optional
        Starting beta coordinate. An alternative to the beta array.
    bstep: float, optional
        Beta coordinate increment. An alternative to the beta array.
    bunit, str, optional
        Beta coordinate units (if different from the default specified in
        the schema).
    
    """
    schema_url = "miri_spectrum3D.schema"
    
    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 wavelength=None, wstart=0.0, wstep=1.0, wunit='',
                 alpha=None, astart=0.0, astep=1.0, aunit='',
                 beta=None, bstart=0.0, bstep=1.0, bunit='',
                 title='', label='', unit='', **kwargs):
        super(Spectrum3D, self).__init__(init=init, data=data, dq=dq, err=err,
                                         dq_def=dq_def,
                                         wavelength=wavelength, wstart=wstart,
                                         wstep=wstep, wunit=wunit, title=title,
                                         label=label, unit=unit, **kwargs)
        # Define which of the axes represents wavelength.
        # TODO: This information should come from the schema.
        self.waveaxis = 0
        
        # The wavelength and beta arrays have already been defined in the
        # above constructor. For 3-D data, the alpha and beta arrays needs
        # to be added. These arrays must match the dimensions of the data
        # array.
        # TODO: These arrays will eventually be obtained from the embedded WCS information.
        if alpha is not None:
            self.alpha = alpha
        elif self.data is not None:
            if (not self._isvalid(self.alpha)) or \
                (np.count_nonzero( self.alpha ) == 0):
                # The alpha array is empty or full of zeros (the default).
                # Fill the array using the start and increment supplied.
                astop = astart + (self.data.shape[2]-1) * astep
                alpha = np.linspace(astart, astop, self.data.shape[2])
                self.alpha = alpha

        # Check that the data and alpha arrays are compatible.
        # (The schema already ensures the alpha array is 1-D)
        if self.data is not None and self.alpha is not None:
            if self.alpha.size != self.data.shape[2]:
                strg = "Alpha array has wrong size "
                strg += "(%d instead of %d)." % \
                    (self.alpha.size, self.data.shape[2])
                raise ValueError(strg)

        if beta is not None:
            self.beta = beta
        elif self.data is not None:
            if (not self._isvalid(self.beta)) or \
                (np.count_nonzero( self.beta ) == 0):
                # The beta array is empty or full of zeros (the default).
                # Fill the array using the start and increment supplied.
                bstop = bstart + (self.data.shape[1]-1) * bstep
                beta = np.linspace(bstart, bstop, self.data.shape[1])
                self.beta = beta

        # Check that the data and beta arrays are compatible.
        # (The schema already ensures the beta array is 1-D)
        if self.data is not None and self.beta is not None:
            if self.beta.size != self.data.shape[1]:
                strg = "Beta array has wrong size "
                strg += "(%d instead of %d)." % \
                    (self.beta.size, self.data.shape[1])
                raise ValueError(strg)

        # If the units of the beta alpha and arrays are explicitly specified,
        # copy the strings to the metadata and the schema. Otherwise obtain
        # the default units from the schema.
        alphaunits = self.set_data_units('alpha', units=aunit)
        betaunits = self.set_data_units('beta', units=bunit)

    def flatten_to_image(self, weighted=True, perfect=False, masked=False):
        """
        
        Flatten a data cube in the wavelength direction to generate one
        combined image.
        
        The output signal is the (weighted) average of the input signal.
        The output error is the RMS of the input error.
        
        :Parameters:
        
        weighted: bool, optional
            Set to True to use the error array (if present) to make a
            weighted average.
            Set to False for an unweighted average.
        perfect: bool, optional
            Set to True to regard as bad any pixel in the image
            which has contributing bad pixel.
            Set to False to regard a pixel in the image as bad only
            if all of the contributing pixels are bad (and use the
            error array to mark the reduction in quality).
        masked: bool, optional
            Set to True to use the masked versions of the data arrays.
            By default, the non-masked versions are used. This is more
            efficient but may generate runtime warning messages.

        :Returns:

        image: MiriMeasuredModel
            A 2-D image object containing the averaged image.        
        
        """
        if masked:
            (flat_data, flat_error, flat_quality) = \
                spec3d_to_image(self.data_masked,
                                self.err_masked, self.dq,
                                waveaxis=self.waveaxis,
                                weighted=weighted, perfect=perfect)
        else:
            (flat_data, flat_error, flat_quality) = \
                spec3d_to_image(self.data, self.err, self.dq,
                                waveaxis=self.waveaxis,
                                weighted=weighted, perfect=perfect)
        title = "(flattened to image)"

        # NOTE: The STScI data model doesn't not support masked arrays,
        # so they need to be converted back to contiguous arrays before
        # they are used to create a new data product.
        image = MiriMeasuredModel(data=np.ascontiguousarray(flat_data),
                                  err=np.ascontiguousarray(flat_error),
                                  dq=np.ascontiguousarray(flat_quality),
                                  title=title)
        return image
    
    def flatten_to_spectrum(self, weighted=True, perfect=False, masked=False):
        """
        
        Flatten a data cube in the spatial directions to generate one
        combined spectrum.
        
        The output signal is the (weighted) average of the input signal.
        The output error is the RMS of the input error.
        
        Not a particularly useful function unless the entire data cube
        contains the spectrum of one object.
    
        :Parameters:
        
        weighted: bool, optional
            Set to True to use the error array (if present) to make a
            weighted average.
            Set to False for an unweighted average.
        perfect: bool, optional
            Set to True to regard as bad any pixel in the spectrum
            which has contributing bad pixel.
            Set to False to regard a pixel in the spectrum as bad only
            if all of the contributing pixels are bad (and use the
            error array to mark the reduction in quality).
        masked: bool, optional
            Set to True to use the masked versions of the data arrays.
            By default, the non-masked versions are used. This is more
            efficient but may generate runtime warning messages.

        :Returns:

        spectrum: Spectrum1D
            A 1-D spectrum object containing the extracted spectrum.        
        
        """
        if masked:
            (flat_data, flat_error, flat_quality) = \
                spec3d_to_spec1d(self.data_masked,
                                 self.err_masked, self.dq,
                                 waveaxis=self.waveaxis,
                                 weighted=weighted, perfect=perfect)
        else:
            (flat_data, flat_error, flat_quality) = \
                spec3d_to_spec1d(self.data, self.err, self.dq,
                                 waveaxis=self.waveaxis,
                                 weighted=weighted, perfect=perfect)
        title = "(flattened to 1-D spectrum)"

        # NOTE: The STScI data model doesn't not support masked arrays,
        # so they need to be converted back to contiguous arrays before
        # they are used to create a new data product.
        spectrum = Spectrum1D(data=np.ascontiguousarray(flat_data),
                              err=np.ascontiguousarray(flat_error),
                              dq=np.ascontiguousarray(flat_quality),
                              wavelength=self.wavelength, title=title)
        return spectrum

    def plot_spectrum(self, plotaxis=None, errorbars=True, description=''):
        """
        
        Plot the spectrum as several pages of a series of overlaid 1-D plots
        (with optional error bars)
        
        """
        # Define fixed labels
        tstrg, xlabel, ylabel = self._get_plot_labels(description)
            
        # Create a new plot for each beta plane.
        for plane in range(0, self.data.shape[1]):
            ptitle = tstrg
            ptitle += " - alpha positions 1 to %d" % self.data.shape[2]
            ptitle += "\n- beta position %d -" % (plane+1)
            
            # Create a matplotlib figure and axis.
            fig = mplt.new_figure(1, stitle=ptitle)
            ax = mplt.add_subplot(fig, 1, 1, 1)

            # The plots will cycle around these defined line formats,
            linefmtlist = ['bo', 'rx', 'g+', 'ko', \
                           'bx', 'r+', 'go', 'kx', \
                           'b+', 'ro', 'gx', 'k+' ]
        
            # Cycle through the list of spectra, remembering that the
            # wavelength is the [0] axis
            # TODO: Can some of this duplicated code be packed up?
            for slc in range(0, self.data.shape[2]):
                spectrum = self.data_filled[:, plane, slc]

                ifmt = slc % len(linefmtlist)
                lfmt = linefmtlist[ifmt] 

                # Plot each spectrum as an XY plot overlaid on the same axis.
                if errorbars and self.has_dataarray('ERR'):
                    errors = self.err_filled[:, plane, slc]
                    mplt.plot_xy(self.wavelength, spectrum, yerr=errors,
                                 plotfig=fig, plotaxis=ax, linefmt=lfmt,
                                 xlabel=xlabel, ylabel=ylabel, title='')
                else:
                    mplt.plot_xy(self.wavelength, spectrum, plotfig=fig,
                                 plotaxis=ax, linefmt=lfmt, xlabel=xlabel,
                                 ylabel=ylabel, title='')
            mplt.show_plot()

    def fov2spax(self, alpha, beta, lamb):
        """
        
        Converts between FOV coordinates (alpha-beta-lambda) and 
        cube spaxel indexes (x->alpha,y->beta,z->wavelength).
    
        :Parameters:
        
        alpha: numpy array
            alpha array
        
        beta: numpy array
            beta array
        
        lamb: numpy array
            wavelength array
        
        :Returns:
        
        (x, y, z): tuple of numpy arrays
            Arrays containing the x, y and z spaxel indices
        
        """
        # The input parameters must be arrays of same shape
        alpha = np.asarray(alpha)
        beta = np.asarray(beta)
        lamb = np.asarray(lamb)
        assert alpha.shape == beta.shape
        assert alpha.shape == lamb.shape
        
        def _trans1D_inv(axis, inaxis):
            assert inaxis.shape[0] > 1
            start = inaxis.min()
            step  = (inaxis.max()-inaxis.min())/(inaxis.shape[0]-1.)
            x = (axis - start)/step
            return x
        
        x = _trans1D_inv(alpha,self.alpha)
        y = _trans1D_inv(beta,self.beta)
        z = _trans1D_inv(lamb,self.wavelength)
    
        return x,y,z
        
    def spax2fov(self, x, y ,z):
        """
        
        Converts between cube spaxel indexes 
        (x->alpha,y->beta,z->wavelength) and FOV coordinates 
        (alpha-beta-lambda).

        :Parameters:
        
        x: numpy array
            x array
        
        y: numpy array
            y array
        
        z: numpy array
            z array
        
        :Returns:
        
        (alpha, beta, lamb)
        
        
        """
        
        from pdb import set_trace as stop
        # The input parameters must be arrays of same shape
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        assert x.shape == y.shape
        assert x.shape == z.shape
        
        # TO DO: start and step should perhaps be attributes of self?
        
        def _trans1D_dir(axis, inaxis):
            
            assert inaxis.shape[0] > 1
            start = inaxis.min()
            step  = (inaxis.max()-inaxis.min())/(inaxis.shape[0]-1.)
            
            y = axis * step + start
            return y
        
        alpha = _trans1D_dir(x,self.alpha)
        beta = _trans1D_dir(y,self.beta)
        lamb = _trans1D_dir(z,self.wavelength)
        
        return alpha,beta,lamb


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the spectral data models\n")

    PLOTTING = False
    SAVE_FILES = False

    a = np.array([10.0,15.0,10.0,20.0,10.0,30.0,90.0,30.0,20.0,50.0,10.0,12.0])
    b = np.array([5.0,4.0,3.0,2.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0])
    w = np.array([5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0,21.0,23.0,25.0,27.0])

    print("Testing Spectrum1D")
    spectrum = Spectrum1D(data=a, err=b, wavelength=w, unit='Jy',
                          title='My favourite 1-D spectrum')
    print(spectrum)
    if PLOTTING:
        spectrum.plot(description='Standard plot')
        spectrum.plot_spectrum(description='Spectrum plot', errorbars=True)
    if SAVE_FILES:
        spectrum.save("test_spectrum1D.fits", overwrite=True)

    # Note that in Spectrum2D, the wavelength array is the second axis.
    # The axes need to be swapped to make the wavelength the long axis.
    a2 = np.ascontiguousarray(np.swapaxes(np.array([a, a+10, a+2]), 0, 1))
    b2 = np.ascontiguousarray(np.swapaxes(np.array([b*2, b, b*2]), 0, 1))
    print("Testing Spectrum2D")
    spimage = Spectrum2D(data=a2, err=b2, wavelength=w, unit='Photons',
                         title='My favourite 2-D spectrum')
    print(spimage)
    if PLOTTING:
        spimage.plot(description='Standard plot')
        spimage.plot_spectrum(description='Spectrum plot', errorbars=True)
    if SAVE_FILES:
        spimage.save("test_spectrum2D.fits", overwrite=True)

    print("Flatten to 1-D spectrum")
    flat_spec = spimage.flatten_to_spectrum(weighted=True, perfect=False)
    print(flat_spec)
    if PLOTTING:
        flat_spec.plot_spectrum(description='Flattened spectrum')

    # Note that in Spectrum2D, the wavelength array is the third axis.
    # The axes need to be swapped to make the wavelength the long axis.
    a3 = np.ascontiguousarray(np.swapaxes(np.array([a2, a2+2, a2+4]), 0, 1))
    b3 = np.ascontiguousarray(np.swapaxes(np.array([b2, b2+1, b2+3]), 0, 1))
    print("Testing Spectrum3D")
    spcube = Spectrum3D(data=a3, err=b3, wavelength=w, unit='Flux units',
                        title='Casper the Spectral Cube')
    print(spcube)
    if PLOTTING:
        spcube.plot(description='Standard plot')
        spcube.plot_spectrum(description='Spectrum plot', errorbars=True)
    if SAVE_FILES:
        spcube.save("test_spectrum3D.fits", overwrite=True)
        
    print("Flatten to 2-D image")
    flat_image = spcube.flatten_to_image(weighted=True, perfect=False)
    print(flat_image)
    if PLOTTING:
        flat_image.plot(description='Image plot')

    print("Flatten to 1-D spectrum")
    flat_spec = spcube.flatten_to_spectrum(weighted=True, perfect=False)
    print(flat_spec)
    if PLOTTING:
        flat_spec.plot_spectrum(description='Spectrum plot')

    print("Test finished")
