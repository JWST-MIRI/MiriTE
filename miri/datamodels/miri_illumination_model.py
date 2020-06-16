#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines a means of 
describing the illumination at a plane within the MIRI instrument.

This module is intended to be used by a MIRI simulator for describing
the illumination seen by the instrument.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

08 Feb 2013: Created (based on the IlluminationDataProduct in the old
             data model).
23 Apr 2013: intensity element renamed data to keep up with changes
             in jwst.datamodel data model (where inclusion of a data
             element is now compulsory). "intensity" alias created.
21 May 2013: Added get and set illumination shape functions.
07 Jun 2013: Remodelled in parallel with SCASim illumination data model.
18 Jun 2013: data element renamed back to intensity, now jwst.datamodels can
             work with data structures without a data element.
             get_primary_array_name() method added.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
11 Apr 2014: Modified to define data units using the set_data_units
             method.
13 Aug 2015: Corrected MiriFilter import error. Corrected bug which caused
             first plane of intensity array to be incremented after each
             call to get_illumination. Added MiriIlluminationFringingModel
             class and miri_illumination_full.schema, to contain a fuller
             description of the illumination needed for fringing models.
08 Sep 2015: Added apply_scale function.
11 Sep 2015: Removed duplicated plot method.
18 Nov 2015: Added truncate() method.
30 Nov 2015: Removed the explicit data typing from the
             get_illumination_enlarged method. Corrected assert statement in
             truncate method.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
03 Jun 2016: Use _isvalid function rather than checking wavelength array
             against None. Check for scalar wavelength special case.
18 Jul 2016: Remove single-dimensional higher axes from the illumination
             data returned by get_illumination.
18 Aug 2016: Guard against testing an undefined data attribute.
24 Aug 2016: Don't try to adjust blank units.
06 Sep 2016: Check for negative or invalid wavelengths before interpolating.
             Replace NaNs within the intensity or wavelength arrays.
20 Jan 2017: Corrected a bug where the maximum wavelength calculation fails
             if the wavelength array is completely full of NaNs.
15 Jun 2017: TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
21 Sep 2017: Test World Coordinates.
21 Jun 2018: Define the FILETYPE keyword.
28 Jun 2018: Switch to using get_title_and_metadata() to display data model
             information.
13 Dec 2019: get_illumination_enlarged function documented (Bug 612).

@author: Steven Beard (UKATC)

"""

import sys
import numpy as np
#import numpy.ma as ma

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel


# List all classes and global functions here.
__all__ = ['MiriIlluminationModel']

# Define a limit within which two floating point wavelengths or directions are
# considered the same.
EPS = 10 * sys.float_info.epsilon


class MiriIlluminationModel(MiriDataModel):
    """
    
    A data model for MIRI illumination information, based on the STScI
    base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the intensity parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    intensity: numpy array (optional)
        An array describing the intensity seen by the MIRI instrument.
        This can be a 2-D array containing a single intensity layer, or
        it can be a 3-D array describing several intensity layers which
        are combined to give the total illumination. Each layer can be
        described by a different layer in the wavelength array.
        Some typical variations:
        
        * A 2-D array describing total intensity. No wavelength information.
        * A 2-D array describing intensity. The wavelength varies across
          the array, as described in a 2-D wavelength array (e.g. for
          spectroscopy).
        * A 3-D array describing intensity. A 1-D wavelength array lists
          the wavelength corresponding to each layer in the intensity
          array (e.g. for imaging through a filter).
        * A 3-D array describing intensity, together with a 3-D array
          describing how the wavelength varies for each layer in the
          intensity array (e.g. for slitless LRS spectroscopy with
          overlapping objects).
        
        If an intensity parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    wavelength: numpy array (optional)
        An array describing the wavelength(s) associated with the layers
        in the intensity array. Its shape must be compatible with the
        intensity array.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_illumination.schema"
    
    def __init__(self, init=None, intensity=None, wavelength=None, **kwargs):
        """
        
        Initialises the MiriIlluminationModel class.
        
        Parameters: See class doc string.

        """
        super(MiriIlluminationModel, self).__init__(init=init, **kwargs)

        # Data type is illumination map.
        self.meta.filetype = 'ILLUMINATION'
        # Initialise the model type
        self._init_data_type()  
            
        if intensity is not None:
            self.intensity = intensity      
        # Replace NaNs or invalid values in the intensity array with zeros.
        if self.intensity is not None:
            self.intensity[np.isnan(self.intensity)] = 0.0

        if wavelength is not None:
            self.wavelength = wavelength
        # Replace NaNs or invalid values in the wavelength array with the maximum.
        if self.wavelength is not None:
            notnanarray = self.wavelength[~np.isnan(self.wavelength)]
            if notnanarray.size > 0:
                wavmax = notnanarray.max()
                self.wavelength[np.isnan(self.wavelength)] = wavmax

        # Only test the wavelength array against the intensity array if the arrays
        # have either been provided in arrays or have been initialised from init.
        if (init is not None) or (intensity is not None and wavelength is not None):
            if not self._check_wavelength(self.intensity, self.wavelength):
                strg = "Shape of wavelength array (%s) " % \
                    str(self.wavelength.shape)
                strg += "is incompatible with intensity array (%s)." % \
                    str(self.intensity.shape)
                raise TypeError(strg)

        # Copy the units of the these arrays, if defined.
        self.set_data_units('intensity')
        self.set_data_units('wavelength')
        
        self._set_illumination_shape()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriIlluminationModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()
        
    def truncate(self, maxshape):
        """
        
        Truncate the illumination map to the given maximum size
        
        """
        assert(isinstance(maxshape, (list,tuple)))
        maxrows = maxshape[0]
        maxcolumns = maxshape[1]
        if len(self.intensity.shape) > 2:
            # 3-D data
            if self.intensity.shape[2] > maxcolumns:
                self.intensity = self.intensity[:,:,:maxcolumns]
            if self.intensity.shape[1] > maxrows:
                self.intensity = self.intensity[:,:maxrows,:]
        else:
            # 2-D data
            if self.intensity.shape[1] > maxcolumns:
                self.intensity = self.intensity[:,:maxcolumns]
            if self.intensity.shape[0] > maxrows:
                self.intensity = self.intensity[:maxrows,:]

        if len(self.wavelength.shape) > 2:
            # 3-D data
            if self.wavelength.shape[2] > maxcolumns:
                self.wavelength = self.wavelength[:,:,:maxcolumns]
            if self.wavelength.shape[1] > maxrows:
                self.wavelength = self.wavelength[:,:maxrows,:]
        else:
            # 2-D data
            if self.wavelength.shape[1] > maxcolumns:
                self.wavelength = self.wavelength[:,:maxcolumns]
            if self.wavelength.shape[0] > maxrows:
                self.wavelength = self.wavelength[:maxrows,:]
 
        if not self._check_wavelength(self.intensity, self.wavelength):
            strg = "Shape of wavelength array (%s) " % \
                str(self.wavelength.shape)
            strg += "is incompatible with intensity array (%s)." % \
                str(self.intensity.shape)
            raise TypeError(strg)
        
        self._set_illumination_shape()

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "intensity"
        and not "data".
        
        """
        return 'intensity'
    
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

    def _nonscalar(self, data):
        """
        
        Helper function to verify that a given array, tuple or list
        has valid content and is also non-scalar.
        
        """
        if not self._isvalid(data):
            return False
        if len(data) <= 1:
            return False
        else:
            return True

    def _set_illumination_shape(self):
        """
        
        Helper function to update the shape of the illumination data.
        
        """
        # The shape of the illumination map (rows, columns) corresponds
        # to the last 2 values in the overall shape of the intensity array,
        # whether it is 2 or 3 dimensional. The number of slices will be
        # non zero only if the intensity map is 3 dimensional.
        if hasattr(self, 'intensity') and self.intensity is not None:
            self._illumination_shape = self.intensity.shape[-2:]
            if len(self.intensity.shape) > 2:
                self.slices = self.intensity.shape[0]
            else:
                self.slices = 0
        else:
            self._illumination_shape = []
            self.slices = 0
            
    def _check_wavelength(self, intensity, wavelength=None):
        """
    
        Verifies that a wavelength array is compatible with the shape
        of the intensity array.
    
        :Parameters:
    
        intensity: array-like
            The intensity data array.
        wavelength: array-like, optional
            If present, the wavelength data array.
    
        :Returns:
    
        valid: bool
            True if the arrays are compatible and False if not compatible.
        
        """
        # The wavelength array can be null or empty - this is ok.
        if intensity is None or wavelength is None:
            return True
        if intensity.size == 0 or wavelength.size == 0:
            return True
    
        # The intensity and wavelength arrays are compatible if their
        # shapes match exactly.
        if wavelength.shape == intensity.shape:
            return True
        
        elif intensity.ndim == 3:
            # Non-matching 3 dimensional data:
            #
            # The wavelength array can be 1-dimensional if that dimension
            # matches the 3rd dimension of the intensity array.
            if wavelength.ndim == 1:
                if wavelength.shape[0] == intensity.shape[0]:
                    return True
            
            # If the wavelength array is 2-dimensional it must match the
            # first 2 dimensions of the 3-D intensity array.
            elif wavelength.ndim == 2:
                if wavelength.shape == intensity.shape[1:]:
                    return True

            # The wavelength array can also be 3-dimensional, match the
            # first 2 dimensions of the 3-D intensity array. and have a
            # 3rd dimension of 1.
            elif wavelength.ndim == 3:
                if wavelength.shape[1:] == intensity.shape[1:] and wavelength.shape[0] == 1:
                    return True
                
                # Also valid is a wavelength array whose first 2 dimensions
                # are 1 but whose 3rd dimension matches that of the intensity
                # array.
                if wavelength.shape[1] == 1 and wavelength.shape[2] == 1:
                    if wavelength.shape[0] == intensity.shape[0]:
                        return True

        # Otherwise the wavelength array is not compatible.
        return False

    def apply_scale(self, scale=1.0):
        """
        
        Apply a scaling factor to the intensity data.
        Wavelength data is unaffacted.
        
        :Parameters:
        
        scale: number
            A scaling factor. The intensity data will be multiplied
            by this number. The default is 1.0.
        
        """
        # Only apply a scaling factor that isn't 1.0
        if abs(scale-1.0) > EPS:
            # Adjust the intensity.
            self.intensity = self.intensity * scale
        
            # Adjust the units to show the multiplication.
            units = self.get_data_units('intensity')
            if units:
                units += " * %.4g" % scale
                self.set_data_units('intensity', units)

    def get_illumination_shape(self):
        """
        
        Get the shape of the integrated illumination map.
        Note that this is the shape of the integrated illumination
        map that would be returned by a call to get_illumination.
        It isn't the shape of the intensity data array.
        
        :Returns:
        
        shape: tuple
            The shape of the illumination map (columns, rows)
            
        """
        return self._illumination_shape

    def get_illumination(self, usefilter=None):
        """
        
        Get the integrated illumination across all wavelengths. If more
        than one plane of illumination data are provided they are
        integrated together, modified by the filter provided. The result
        is a 2-D illumination array.
        
        The intensity data is returned in the same units as the original
        data.
                
        :Parameters:
        
        usefilter: Filter, optional (default=None)
            A Filter object describing the transmission or QE as a function
            of wavelength to be applied to the illumination data. If None,
            no filter is applied.
            A filter is only valid when an illumination data product
            includes a wavelength array.
            NOTE: Filtering by wavelength is only possible when the
            miri.datamodels.filters module is available. See the miri.datamodels
            documentation for a description of the Filter object.
                
        :Returns:
        
        illumination: array-like
            A 2-D array containing the integrated illumination.
            
        """
        if self.intensity.ndim == 2:
            # The intensity array is a single image which may or may not
            # have an associated wavelength.
            if usefilter is None or not self._isvalid(self.wavelength) \
               or self.wavelength.any() <= 0.0:
                # No filter provided or invalid wavelength information.
                return self.intensity
            else:
                # Apply the filter.
                return usefilter.apply_filter(self.intensity, self.wavelength)
        else:
            # The intensity array consists of several planes
            # which need to be combined together into a single
            # 2-D illumination map.
            if usefilter is None or not self._isvalid(self.wavelength) \
               or self.wavelength.any() <= 0.0:
                # No filter provided or invalid wavelength information.
                # Just sum the planes together.
                illumination = self.intensity[0,:,:].copy()
                for planenum in range(1,self.intensity.shape[0]):
                    illumination += self.intensity[planenum,:,:].copy()
                # Remove single-dimensional axes higher than the 3rd dimension,
                # to ensure the resulting illumination data is always 2-D.
                while illumination.ndim > 2:
                    illumination = illumination.squeeze(axis=0)
                return illumination
            else:
                # Apply the filter to each plane (assuming the wavelength
                # array is broadcastable) and sum the results together.
                img = self.intensity[0,:,:]
                if len(self.wavelength) > 1:
                    wav = self.wavelength[0,:,:]
                else:
                    # Scalar special case
                    wav = self.wavelength
                illumination = usefilter.apply_filter(img, wav)
                for planenum in range(1,self.intensity.shape[0]):
                    img = self.intensity[planenum,:,:]
                    if len(self.wavelength) > 1:
                        wav = self.wavelength[planenum,:,:]
                    else:
                        # Scalar special case
                        wav = self.wavelength
                    illumination += usefilter.apply_filter(img, wav)
                # Remove single-dimensional axes higher than the 3rd dimension,
                # to ensure the resulting illumination data is always 2-D.
                while illumination.ndim > 2:
                    illumination = illumination.squeeze(axis=0)
                return illumination

    def get_illumination_enlarged(self, newshape, usefilter=None,
                                  location=(1,1), leftcrop=0, rightcrop=0 ):
        """
        
        Get the integrated illumination across all wavelengths, as in
        method get_illumination, but enlarge the resulting illumination
        data to the specified shape and size, placing the old data at
        the given location within the enlarged data and padding the outside
        areas with zeros. 
        
        :Parameters:
        
        newshape: tuple of 2 ints
            The new size and shape for the data (rows, columns)
        usefilter: Filter
            A Filter object describing the transmission or QE as a function
            of wavelength to be applied to the illumination data. If None,
            no filter is applied.
            See the miri.datamodels.Filter documentation for a description of
            the Filter object.
        location: tuple of 2 ints, optional, default=(1,1)
            The location at which to place the bottom, left corner of
            the old data (row, column). Rows and columns start at 1.
            If not specified, the old data is placed at the bottom
            left corner.
        leftcrop: int, optional, default=0
            The number of cropped columns at the left hand edge.
            The location coordinates have this number of columns
            subtracted, and the illumination array is cropped if
            it is placed within this cropped zone.
        rightcrop: int, optional, default=0
            The number of cropped columns at the right hand edge.
            The illumination array is cropped if it is placed within
            this cropped zone.
                
        :Returns:
        
        illumination_enlarged: array_like
            A 2-D array containing the integrated illumination enlarged
            to the new shape.
        
        """
        # The new size must be larger than the original illumination map.
        try:
            newrows = int(newshape[0])
            newcols = int(newshape[1])
            # If more than 2 ints are given, the rest are ignored.
        except ValueError or TypeError or IndexError:
            raise TypeError("newshape should be a tuple of 2 ints.")
        
        if (newrows < self._illumination_shape[0]) or \
           (newcols < self._illumination_shape[1]):
            strg = "Enlarged illumination shape (%d x %d) " % newshape
            strg += "is smaller than the current shape (%d x %d)." % \
                self._illumination_shape
            raise ValueError(strg)
        
        # First get the original-sized illumination map.
        illumination = self.get_illumination(usefilter=usefilter)
        
        # Create a zeroed 2-D array to hold the new size of illumination map.
        # Make sure the new array has the same data type as the original one.
        new_illumination = np.zeros([newrows, newcols], dtype=illumination.dtype)
                
        # Determine the slice which the original illumination data
        # will occupy within the new data. 1 is subtracted because
        # locations start at 1 but array indices start at 0.
        # The coordinates are shifted to the left by leftcrop columns.
        # TODO: Bug 612. Check this placement is done correctly.
        try:
            r1 = int(location[0]) - 1
            c1 = int(location[1]) - 1 - leftcrop
            # If more than 2 ints are given, the rest are ignored.
        except ValueError or TypeError or IndexError:
            raise TypeError("location should be a tuple of 2 ints.")

        r2 = r1 + self._illumination_shape[0]
        c2 = c1 + self._illumination_shape[1]
        
        # If c1 is placed beyond the left edge of the area, the
        # illumination data must be cropped.
        if c1 < 0:
            lcrop = -c1
            c1 = 0
        else:
            lcrop = 0
            
        # If c2 is placed beyond the right edge of the area the
        # illumination data must again be cropped.
        if c2 > newcols:
            rcrop = newcols - c2
            c2 = newcols
        else:
            rcrop = 0
        
        # The original illumination data is only sliced when leftcrop
        # and rightcrop are non zero.
        if lcrop == 0 and rcrop == 0:
            # No cropping needed. Place the illumination data directly
            # into the flux array.
            new_illumination[r1:r2, c1:c2] = illumination
        elif rcrop == 0:
            # Left cropping only.
            new_illumination[r1:r2, c1:c2] = illumination[:, lcrop:]
        elif lcrop == 0:
            # Right cropping only.
            new_illumination[r1:r2, c1:c2] = illumination[:, :rcrop]
        else:
            # Both left and right cropping.
            new_illumination[r1:r2, c1:c2] = illumination[:, lcrop:rcrop]
                            
        del illumination
        return new_illumination

    def __str__(self):
        """
        
        Return the contents of the illumination map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()
        
        # Add the intensity and wavelength information.
        strg += self.get_data_str('intensity', underline=True, underchar="-")
        strg += self.get_data_str('wavelength', underline=True, underchar="-")
        return strg


class MiriIlluminationFringingModel(MiriIlluminationModel):
    """
    
    A data model for MIRI illumination information which includes the
    direction information needed to describe the fringing model, based
    on the STScI base model, DataModel.
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the intensity parameter, the shape is derived from the
          array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    intensity: numpy array (optional)
        An array describing the intensity seen by the MIRI instrument.
        This can be a 2-D array containing a single intensity layer, or
        it can be a 3-D array describing several intensity layers which
        are combined to give the total illumination. Each layer can be
        described by a different layer in the wavelength array.
        Some typical variations:
        
        * A 2-D array describing total intensity. No wavelength information.
        * A 2-D array describing intensity. The wavelength varies across
          the array, as described in a 2-D wavelength array (e.g. for
          spectroscopy).
        * A 3-D array describing intensity. A 1-D wavelength array lists
          the wavelength corresponding to each layer in the intensity
          array (e.g. for imaging through a filter).
        * A 3-D array describing intensity, together with a 3-D array
          describing how the wavelength varies for each layer in the
          intensity array (e.g. for slitless LRS spectroscopy with
          overlapping objects).
        
        If an intensity parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    wavelength: numpy array (optional)
        An array describing the wavelength(s) associated with the layers
        in the intensity array. Its shape must be compatible with the
        intensity array.
    direction: numpy array (optional)
        An array describing the directions associated with the layers
        in the intensity array. Its shape must be compatible with the
        intensity array.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_illumination_full.schema"
    
    def __init__(self, init=None, intensity=None, wavelength=None,
                 direction=None, **kwargs):
        """
        
        Initialises the MiriIlluminationFringingModel class.
        
        Parameters: See class doc string.

        """
        super(MiriIlluminationFringingModel, self).__init__(init=init,
                                                            intensity=intensity,
                                                            wavelength=wavelength,
                                                            **kwargs)

        if direction is not None:
            self.direction = direction

        if not self._check_direction(self.intensity, self.direction):
            strg = "Shape of direction array (%s) " % \
                str(self.direction.shape)
            strg += "is incompatible with intensity array (%s)." % \
                str(self.intensity.shape)
            raise TypeError(strg)

        # Copy the units of the these arrays, if defined.
        self.set_data_units('direction')
 
    def _check_direction(self, intensity, direction=None):
        """
    
        Verifies that a direction array is compatible with the shape
        of the intensity array.
    
        :Parameters:
    
        intensity: array-like
            The intensity data array.
        direction: array-like, optional
            If present, the direction data array.
    
        :Returns:
    
        valid: bool
            True if the arrays are compatible and False if not compatible.
        
        """
        #print("check_direction:", intensity, direction)
        # The direction array can be null or empty - this is ok.
        if intensity is None or direction is None:
            return True
        if intensity.size == 0 or direction.size == 0:
            return True

        # The intensity and direction arrays are also compatible if the
        # direction array is defaulted.
        if direction.shape == intensity.shape:
            if np.mean(direction) < EPS:
                return True

        # The tests depend on whether the intensity array is structured
        # in layers or not.
        if intensity.ndim == 2:
            # No layers. The direction array must be 3-D and define 2 angles,
            # and the number of rows and columns must match.
            if direction.ndim == 3:
                if direction.shape[-1] == intensity.shape[-1] and \
                   direction.shape[-2] == intensity.shape[-2] and \
                   direction.shape[-3] == 2:
                    return True
        
        elif intensity.ndim == 3:
            # Layers. Make the same tests as the single layer case, except
            # the number of layers must also match, or be 1.
            if direction.ndim == 4:
                if direction.shape[-1] == intensity.shape[-1] and \
                   direction.shape[-2] == intensity.shape[-2] and \
                   direction.shape[-3] == 2:
                    if direction.shape[-4] == intensity.shape[-3] or \
                       direction.shape[-4] == 1:
                        return True
            elif direction.ndim == 3:
                if direction.shape[-1] == intensity.shape[-1] and \
                   direction.shape[-2] == intensity.shape[-2] and \
                   direction.shape[-3] == 2:
                    return True

        # Otherwise the direction array is not compatible.
        return False

    def __str__(self):
        """
        
        Return the contents of the illumination map with fringing object
        as a readable string.
        
        """
        # Start with the underlying illumination map.
        strg = super(MiriIlluminationFringingModel, self).__str__()

        # Add the direction information
        strg += self.get_data_str('direction', underline=True, underchar="-")
        return strg

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriIlluminationModel module.")
    try:
        from miri.datamodels import miri_filters
    except ImportError:
        print("MIRI filters module could not be imported - " + \
            "there will be no filtering of IlluminationDataProduct.")
        miri_filters = None
    
    PLOTTING = False
    SAVE_FILES = False
 
    ii1 = [[10,20,30,40], [50,60,70,80], [90,100,110,120]]
    ii3 = [ii1,ii1,ii1]
    ww = [[1,2,3,4], [5,6,7,8], [9,10,11,12]]
    ww3 = [ww,ww,ww]   
    dd1 = [[45,45,45,45], [45,45,45,45], [45,45,45,45]]
    dd2 = [[135,135,135,135], [135,135,135,135], [135,135,135,135]]
    dd3 = [dd1, dd2]
    dd4 = [dd3,dd3,dd3]

    print("Illumination data with intensity but no wavelength")
    with MiriIlluminationModel( intensity=ii3 ) as ill1:
        print(ill1)
        if PLOTTING:
            ill1.plot()
        if SAVE_FILES:
            ill1.save("test_illumination_model0.fits", overwrite=True)

        # Attempt to filter the illumination data with wavelength.
        # Nothing should happen. There should be no errors.
        if miri_filters is not None:
            wavelengths = np.linspace(1.0, 20.0, 180)
            transmissions = 0.5*np.ones((180,))
            filter_table = []
            for w,t in zip(wavelengths, transmissions):
                filter_table.append( (w,t) )
            flt = miri_filters.MiriFilter(filter_table=filter_table,
                                          filter_name='F560W',
                                          filter_type='BandPass')
        else:
            print("No filtering with wavelength")
            flt = None

        ill1shape = ill1.get_illumination_shape()
        light1 = ill1.get_illumination(usefilter=None)
        print("Illumination of shape %s (without filter) =\n%s" % \
              (str(ill1shape), str(light1)))
        del light1
        if flt is not None:
            light2 = ill1.get_illumination(usefilter=flt)
            filter_title = flt.get_title()
            print("Illumination of shape %s (with filter \'%s\') =\n%s" % \
                  (str(ill1shape),filter_title, str(light2)))
            del light2
        del ill1

    print("Illumination data with intensity and wavelength")
    with MiriIlluminationModel( intensity=ii3, wavelength=ww3 ) as ill2:
        ill2.set_wcs_metadata( wcsaxes=2 )
        print(ill2)
        if PLOTTING:
            ill2.plot()
        if SAVE_FILES:
            ill2.save("test_illumination_model1.fits", overwrite=True)

        # Filter the illumination data with wavelength.
        if miri_filters is not None:
            wavelengths = np.linspace(1.0, 20.0, 180)
            transmissions = 0.5*np.ones((180,))
            filter_table = []
            for w,t in zip(wavelengths, transmissions):
                filter_table.append( (w,t) )
            flt = miri_filters.MiriFilter(filter_table=filter_table,
                                          filter_name='F560W',
                                          filter_type='BandPass')
        else:
            print("No filtering with wavelength")
            flt = None
        ill2shape = ill2.get_illumination_shape()
        light1 = ill2.get_illumination(usefilter=None)
        print("Illumination of shape %s (without filter) =\n%s" % \
              (str(ill2shape), str(light1)))
        del light1
        if flt is not None:
            light2 = ill2.get_illumination(usefilter=flt)
            filter_title = flt.get_title()
            print("Illumination of shape %s (with filter \'%s\') =\n%s" % \
                  (str(ill2shape),filter_title, str(light2)))
            del light2

        lightbig = ill2.get_illumination_enlarged(newshape=(6,8),
                                                 usefilter=None, location=(2,2))
        print("Enlarged illumination = without filter \n" + str(lightbig))

        del ill2
        if flt is not None:
            del flt

    with MiriIlluminationFringingModel( intensity=ii3, wavelength=ww3,
                                        direction=dd4 ) as ill:
        # Try non-default units.
        ill.set_data_units('intensity', 'electrons/s')
        print(ill)
        if PLOTTING:
            ill.plot()
        if SAVE_FILES:
            ill.save("test_illumination_model2.fits", overwrite=True)

        if miri_filters is not None:
            wavelengths = np.linspace(1.0, 20.0, 180)
            transmissions = 0.5*np.ones((180,))
            filter_table = []
            for w,t in zip(wavelengths, transmissions):
                filter_table.append( (w,t) )
            flt = miri_filters.MiriFilter(filter_table=filter_table,
                                          filter_name='F560W',
                                          filter_type='BandPass')
        else:
            print("No filtering with wavelength")
            flt = None
        illshape = ill.get_illumination_shape()
        light1 = ill.get_illumination(usefilter=None)
        print("Illumination of shape %s (without filter) =\n%s" % \
              (str(illshape), str(light1)))
        del light1
        if flt is not None:
            light2 = ill.get_illumination(usefilter=flt)
            filter_title = flt.get_title()
            print("Illumination of shape %s (with filter \'%s\') =\n%s" % \
                  (str(illshape),filter_title, str(light2)))
            del light2

        lightbig = ill.get_illumination_enlarged(newshape=(6,8),
                                                 usefilter=None, location=(2,2))
        print("Enlarged illumination = without filter \n" + str(lightbig))

        del ill
        if flt is not None:
            del flt

    print("Test finished")
