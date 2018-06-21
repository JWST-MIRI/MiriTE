#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model for MIRI/JWST point spread 
function data. Essentially the same as the MIRI measured image data 
model, with additional metadata and an optional lookup table.

The spectroscopy PSF models contain a data cube describing the change in PSF
with wavelength and field position and use world coordinates metadata to map
each PSF plane to a wavelength and field position.

The imaging PSF models contain a data cube of PSFs stacked by filter
and use a lookup table to map each PSF plane to a filter and field position.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

10 Jan 2013: Created
21 Jan 2013: Included data type, wavelength and focal plane X,Y in metadata.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: 3-D data is allowed. Replaced 'to_fits' with more generic
             'save' method.
13 Feb 2013: Added psf_lut table.
26 Feb 2013: Only update the wavelength, xfield and yfield if explicitly
             provided.
17 May 2013: Create a default table object explicitly to stop
             the underlying data model filling the table with
             junk by default.
22 Aug 2013: columnnames renamed to fieldnames.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
10 Oct 2013: psf_lut table removed from MiriPointSpreadFunctionModel. It now
             only applies to MiriImagingPointSpreadFunctionModel.
             Tables containing a single row are possible. Removed the default
             psf_lut table.
30 Oct 2013: Renamed the schemas so they have adjacent listings. Separate
             PSF models for imaging, LRS and MRS.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
06 Mar 2014: Changed type back to 'PSF' for IM and MRS.
28 Mar 2014: Wavelength extension removed MiriLrsPointSpreadFunctionModel.
             Changed type back to 'PSF' for LRS.
09 Jul 2014: field_def changed to dq_def.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
17 Oct 2014: Updated the reference flags.
30 Oct 2014: WAVECENT keyword and FILTER column removed from schema.
09 Dec 2014: Removed old commented out warnings.
16 Jan 2015: Ensure the main data array is 2-D or 3-D when specified
             explicitly or when read from a file.
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
20 Nov 2015: Added PSF-OOF reference type option added for imager PSFs.
03 Dec 2015: Added 4 columns to the imager PSF_LUT table: COL_FIELD, ROW_FIELD,
             XAN_FIELD, YAN_FIELD.
08 Dec 2015: LRS PSF allowed to be 3-D.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
15 Jul 2016: Removed obsolete flags.
18 Aug 2016: Guard against testing an undefined data attribute.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""
# This module is now converted to Python 3.


#import warnings
import numpy as np
#import numpy.ma as ma

# Import the MIRI measured data model.
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['MiriPointSpreadFunctionModel', 'MiriImagingPointSpreadFunctionModel',
           'MiriLrsPointSpreadFunctionModel', 'MiriMrsPointSpreadFunctionModel']

psf_reference_setup = \
            [(0, 'DO_NOT_USE',         'Bad pixel. Do not use.'),
             (1, 'CDP_LOW_QUAL',       'Data of low quality'),
             (2, 'CDP_UNRELIABLE_ERROR', 'Data without reliable error estimate')]
psf_reference_flags = insert_value_column( psf_reference_setup )


class MiriPointSpreadFunctionModel(MiriMeasuredModel):
    """
    
    A data model for MIRI generic point spread function data.
        
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
        An array containing the PSF data. Either a single 2-D array
        containing one PSF at a specific (wavelength, xfield, yfield)
        or a 3-D array containing a stack of PSFs identified by a
        lookup table.
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
        example: [(0, 'good','Good data'), (1, 'dead', 'Dead Pixel'),
        (2, 'hot', 'Hot Pixel')];
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to the MIRI reserved flags.
    wavelen: float, optional
        The wavelength for which these PSF data are valid.
        Only sensible if the PSF data all apply to the same wavelength
        (otherwise use a lookup table).
    pixelsize: float, optional
        The pixel size for which these PSF data are valid.
    xfield: float, optional
        The X focal plane position for which these PSF data are valid.
        Only sensible if the PSF data all apply to the same xfield
        (otherwise use a lookup table).
    yfield: float, optional
        The X focal plane position for which these PSF data are valid.
        Only sensible if the PSF data all apply to the same yfield
        (otherwise use a lookup table).
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
    
    """
    schema_url = "miri_psf.schema.yaml"
    _default_dq_def = psf_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 wavelen=None, pixelsize=None, xfield=None, yfield=None,
                 **kwargs):
        """
        
        Initialises the MiriPointSpreadFunctionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPointSpreadFunctionModel, self).__init__(init=init, data=data,
                                                    dq=dq, err=err,
                                                    dq_def=dq_def, **kwargs)

        # Data type is PSF.
        if not self.meta.reftype:
            self.meta.model_type = 'PSF'
            self.meta.reftype = 'PSF'
        
        # The default pedigree is 'GROUND'
        if not self.meta.pedigree:
            self.meta.pedigree = 'GROUND'
            
        # A USEAFTER date must exist. If not relevant, set it to an
        # impossibly early date.
        if not self.meta.useafter:
            self.meta.useafter = '2000-01-01T00:00:00'

        # The main data array should be 2-D or 3-D.
        # TODO: Can this check be included in the schema?
        if data is not None:
            if not hasattr(data, 'ndim'):
                data = np.asarray(data)
            if data.ndim < 2 or data.ndim > 3:
                strg = "The main data array in a PSF object must be "
                strg += "2-D or 3-D. %d-D data provided" % data.ndim
                raise TypeError(strg)
        elif init is not None and self.data is not None and len(self.data) > 0 and \
             hasattr(self.data, 'ndim'):
            if self.data.ndim < 2 or self.data.ndim > 3:
                strg = "The main data array in a PSF object must be "
                strg += "2-D or 3-D. %d-D data provided" % self.data.ndim
                raise TypeError(strg)
        
        # Wavelength, pixel size and (xfield,yfield) are added to the metadata if present.
        if wavelen is not None:
            self.meta.wavelength_center = wavelen
        if pixelsize is not None:
            self.meta.pixel_size = pixelsize
        if xfield is not None and yfield is not None:
            self.meta.xfield = xfield
            self.meta.yfield = xfield

    def __str__(self):
        """
        
        Return the contents of the PSF object as a readable
        string.
        
        """
        # First obtain a string describing the underlying measured
        # model.
        strg = super(MiriPointSpreadFunctionModel, self).__str__()
        
        # Add the extras
        if self.meta.wavelength_center is not None:
            strg += "Data valid for wavelength=%.2f microns " % self.meta.wavelength_center
        else:
            strg += "Data valid for a variety of wavelengths "
        if self.meta.xfield is not None and \
           self.meta.yfield is not None:
            strg += "and focal plane position=(%f,%f).\n" % \
                (self.meta.xfield, self.meta.yfield)
        else:
            strg += "and a variety of focal plane positions.\n"
        return strg


class MiriImagingPointSpreadFunctionModel(MiriPointSpreadFunctionModel):
    """
    
    A data model for MIRI imaging PSF data, based on
    MiriPointSpreadFunctionModel, with the PSFs looked up
    by filter name.
    
    :Parameters:
    
    The same as MiriPointSpreadFunctionModel with the addition of
    
    psf_lut:  list of tuples or numpy record array (optional)
        Either: A list of tuples containing (xfield:number, yfield:number,
        stack:integer, col_field:number, row_field:number, xan_field:number,
        yan_field:number), giving the mapping between field position and location in the PSF stack.
        Or: A numpy record array containing the same information as above.
    psftype: str (optional)
        A string giving the specific kind of imager PSF contained
        in the data object: 'PSF' or 'PSF-OOF' (out of field PSF).
        If not given, the type defaults to the generic term 'PSF'.
    wavelen: float, optional
        The wavelength for which these PSF data are valid.
        Only sensible if the PSF data all apply to the same wavelength
        (otherwise use a lookup table).
    pixelsize: float, optional
        The pixel size for which these PSF data are valid.
    
    """
    schema_url = "miri_psf_imaging.schema.yaml"
    fieldnames = ('XFIELD', 'YFIELD', 'STACK', 'COL_FIELD', 'ROW_FIELD', \
                  'XAN_FIELD', 'YAN_FIELD')

    def __init__(self, init=None, psf_lut=None, psftype='', wavelen=None,
                 pixelsize=None, **kwargs):
        """
        
        Initialises the MiriImagingPointSpreadFunctionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriImagingPointSpreadFunctionModel, self).__init__(init=init,
                                        wavelen=wavelen, pixelsize=pixelsize,
                                        **kwargs)

        # Data type is imager PSF.
        if not psftype:
            if not self.meta.reftype:
                self.meta.model_type = 'PSF (Imaging)'
                self.meta.reftype = 'PSF'
        else:
            self.meta.model_type = psftype
            psupper = psftype.upper()
            if 'OOF' in psupper:
                self.meta.reftype = 'PSF-OOF'
            else:
                self.meta.reftype = 'PSF'
        
        if psf_lut is not None:
            try:
                self.psf_lut = psf_lut
            except (ValueError, TypeError) as e:
                strg = "psf_lut must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
# 
#         # Copy the table column units, if defined.
#         psf_units = self.set_table_units('psf_lut')

    def __str__(self):
        """
        
        Return the contents of the imaging PSF object as a readable
        string.
        
        """
        # First obtain a string describing the underlying PSF model.
        strg = super(MiriImagingPointSpreadFunctionModel, self).__str__()
        
        if self.psf_lut is not None:
            strg += self.get_data_str('psf_lut', underline=True, underchar="-")
        else:
            strg += "No psf_lut."
        return strg


class MiriLrsPointSpreadFunctionModel(MiriPointSpreadFunctionModel):
    """
    
    A data model for MIRI LRS PSF data, based on
    MiriPointSpreadFunctionModel. 
    
    :Parameters:
    
    The same as MiriPointSpreadFunctionModel with the addition of
    parameters CRVAL1, CRVAL2, CRVAL3 denoting the reference pixel(s)
    of dispersion zero point and
    
    psftype: str (optional)
        A string giving the specific kind of LRS PSF contained
        in the data object: 'PSF' or 'PSF-MONOCHROM' (monochomatic PSF).
        If not given, the type defaults to the generic term 'PSF'.
    
    """
    schema_url = "miri_psf_lrs.schema.yaml"

    def __init__(self, init=None, psftype=None, **kwargs):
        """
        
        Initialises the MiriLrsPointSpreadFunctionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsPointSpreadFunctionModel, self).__init__(init=init,
                                                **kwargs)

        # Data type is LRS PSF.
        if not psftype:
            if not self.meta.reftype:
                self.meta.model_type = 'PSF (LRS)'
                self.meta.reftype = 'PSF'
        else:
            self.meta.model_type = psftype
            psupper = psftype.upper()
            if 'MONOCHROM' in psupper:
                self.meta.reftype = 'PSF-MONOCHROM'
            else:
                self.meta.reftype = 'PSF'

    def __str__(self):
        """
        
        Return the contents of the LRS PSF object as a readable
        string.
        
        """
        # First obtain a string describing the underlying PSF model.
        strg = super(MiriLrsPointSpreadFunctionModel, self).__str__()
        
        return strg


class MiriMrsPointSpreadFunctionModel(MiriPointSpreadFunctionModel):
    """
    
    A data model for MIRI MRS PSF data, based on
    MiriPointSpreadFunctionModel, with the wavelength information for
    the PSFs provided in world coordinates keywords.
    
    :Parameters:
    
    The same as MiriPointSpreadFunctionModel plus
    
    World coordinates?
    
    """
    schema_url = "miri_psf_mrs.schema.yaml"

    def __init__(self, init=None, wavelength=None, **kwargs):
        """
        
        Initialises the MiriMrsPointSpreadFunctionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsPointSpreadFunctionModel, self).__init__(init=init,
                                                **kwargs)
        # Data type is MRS PSF.
        self.meta.model_type = 'PSF (MRS)'
        self.meta.reftype = 'PSF'

    def __str__(self):
        """
        
        Return the contents of the MRS PSF object as a readable
        string.
        
        """
        # First obtain a string describing the underlying PSF model.
        strg = super(MiriMrsPointSpreadFunctionModel, self).__str__()
        
        # World coordinates ?
        
        return strg


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriPointSpreadFunctionModel module.")

    PLOTTING = False
    SAVE_FILES = False

    data3x3 = np.array([[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]])
    err3x3 = np.array([[1.,1.,1.],[2.,2.,2.],[1.,1.,1.]])
    dq3x3 = np.array([[0,1,0],[1,0,1],[0,1,0]])
        
    data4x3x3 = [data3x3,data3x3,data3x3,data3x3]

#     lut_im = [('F560W',  1.0,  0.0, 0),
#               ('F560W',  1.1,  0.0, 1),
#               ('F1000W', 1.0,  0.0, 2),
#               ('F1000W', 1.1,  0.0, 3)
#               ]

    lut_im = [(1.0,  0.0, 0, 0.7, 0.7, 0.4, 0.4),
              (1.1,  0.0, 1, 0.5, 0.5, 0.2, 0.2),
              (1.0,  0.0, 2, 0.1, 0.1, 0.8, 0.8),
              (1.1,  0.0, 3, 0.6, 0.6, 0.1, 0.1)
              ]

    print("2-D PSF data with data + err + dq:")
    with MiriPointSpreadFunctionModel(data=data3x3, err=err3x3, dq=dq3x3,
                                      dq_def=psf_reference_flags,
                                      wavelen=12.0, xfield=4.0, yfield=3.0) \
            as testdata1:
        print(testdata1)
        if PLOTTING:
            testdata1.plot(description="PSF testdata1")
        if SAVE_FILES:
            testdata1.save("test_psf_model1.fits", overwrite=True)
        del testdata1

    print("3-D PSF data with data + err + dq:")
    with MiriPointSpreadFunctionModel(data=data4x3x3, err=err3x3, dq=dq3x3,
                                      dq_def=psf_reference_flags) \
            as testdata2:
        print(testdata2)
        if PLOTTING:
            testdata2.plot(description="PSF testdata2")
        if SAVE_FILES:
            testdata2.save("test_psf_model2.fits", overwrite=True)
        del testdata2
        
    # Read back the file using the imaging PSF class and save it
    # to a new file
    if SAVE_FILES:
        print("3-D PSF data read back as an imaging PSF:")
        with MiriImagingPointSpreadFunctionModel(init="test_psf_model2.fits") \
                as testdata2a:
            print(testdata2a)
            testdata2a.save("test_psf_model2_copy.fits", overwrite=True)
            del testdata2a
    
    print("3-D imaging PSF data with data + err + dq + psf_lut:")
    with MiriImagingPointSpreadFunctionModel(data=data4x3x3, err=err3x3, dq=dq3x3,
                                      dq_def=psf_reference_flags, psf_lut=lut_im) \
            as testdata3:
        print(testdata3)
        if PLOTTING:
            testdata3.plot(description="Imaging PSF testdata3")
        if SAVE_FILES:
            testdata3.save("test_psf_model3.fits", overwrite=True)
        del testdata3

    print("3-D LRS PSF data with data + err + dq + wavelength:")
    with MiriLrsPointSpreadFunctionModel(data=data4x3x3, err=err3x3, dq=dq3x3,
                                      dq_def=psf_reference_flags) \
            as testdata3:
        print(testdata3)
        if PLOTTING:
            testdata3.plot(description="LRS PSF testdata4")
        if SAVE_FILES:
            testdata3.save("test_psf_model4.fits", overwrite=True)
        del testdata3

    print("3-D MRS PSF data with data + err + dq:")
    with MiriMrsPointSpreadFunctionModel(data=data4x3x3, err=err3x3, dq=dq3x3,
                                      dq_def=psf_reference_flags) \
            as testdata3:
        print(testdata3)
        if PLOTTING:
            testdata3.plot(description="MRS PSF testdata4")
        if SAVE_FILES:
            testdata3.save("test_psf_model5.fits", overwrite=True)
        del testdata3

    print("Test finished.")    
