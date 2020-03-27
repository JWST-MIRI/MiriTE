#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

An extension to the standard STScI data model, which defines a means of 
describing MIRI distortion coefficients.

NOTE: The contents of this data model might change, depending on the 
STScI implementation of distortion models.

:Reference:

The STScI jwst.datamodels documentation. See
https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html

:History:

21 Jan 2013: Created
23 Jan 2013: Added plotting.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
21 Feb 2013: Changed default order of MiriDistortionModel from "3" to 
             "None" (Vincent Geers, DIAS)
25 Feb 2013: Corrected typo.
26 Feb 2013: Changed the default fit type from 'POLY2D' to None.
01 Jul 2013: get_primary_array_name() method added.
12 Sep 2013: Change the way the row and column matrices are checked
             so the .copy() method works.
13 Sep 2013: Changed CMATRIX and RMATRIX to BMATRIX and AMATRIX, added 
             new TMATRIX and MMATRIX.
13 Sep 2013: Corrected some typos.
16 Sep 2013: Removed the ORDER parameter, since it is derivable from
             the size of the matrices. BUNIT parameters added to schema
             metadata.
30 Oct 2013: All MIRI distortion models (imaging, LRS, MRS) combined
             into one module. New model for LRS distortion and wavelength
             calibration.
31 Oct 2013: BETA array removed from MRS D2C model.
27 Nov 2013: Modified to match update to distortion table definition in
             miri.distortion.lrs.schema
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
10 Apr 2014: Modified for jsonschema draft 4: Functions made more
             independent of schema structure. Modified to define data
             units using the set_data_units method.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: TYPE and REFTYPE are no longer identical.
07 Oct 2014: Added new inverse matrices BIMATRIX, AIMATRIX, TIMATRIX, MIMATRIX
             and the new BORESIGHT_OFFSETS table. Removed need for fitref, 
             changed fitmodel to now contain reference to documentation.
10 Oct 2014: Restored fitref for consistency with the linearity model.
16 Oct 2014: REFTYPE of WCS changed to DISTORTION.
02 Jul 2015: Major change to MRS distortion model. MiriMrsD2CModel (data array
             containing looking tables) replaced by MiriMrsDistortionModel
             (tables containing polynomial coefficients).
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
             Separate data models created for channel 12 and channel 34
             distortion data. (Merge these models after CDP-4 delivery.)
11 Sep 2015: Removed duplicated plot method.
17 Nov 2015: Changed column names and HDU names to eliminate fitsverify
             problems. NEW DATA MODELS TO BE INSTATED AFTER CDP-5 DELIVERY.
10 Dec 2015: Old and new data models merged into one module.
11 Dec 2015: v2v3 changed to XANYAN in new MRS distortion models.
16 Feb 2016: Imager distortion matrices changed to float64.
09 Jun 2016: Added set_exposure_type() call to MiriImagingDistortionModel,
             MiriLrsD2WModel, MiriMrsDistortionModel12, and 
             MiriMrsDistortionModel34 to set the EXP_TYPE keyword (now 
             required for DISTORTION files).
16 Jun 2016: Added new metadata keywords to MRS distortion schemas.
             Old format MRS data models removed (as MIRISim no longer
             uses them).
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
10 Aug 2018: Updated MRS distortion models to reflect CDP-7 format.
03 Sep 2018: Old CDP-6 variants of the distortion models included.
             Updated units for imager distortion model.
14 Nov 2018: Explicitly set table column units based on the tunit definitions
             in the schema. Removed redundant function.
30 Jan 2019: self.meta.model_type now set to the name of the STScI data
             model this model is designed to match (skipped if there isn't
             a corresponding model defined in ancestry.py).
11 Feb 2019: Added missing C, D, E, F matrices to imager distortion model.
22 Mar 2019: Changed the reference type for LRS distortion
             from 'DISTORTION' to 'SPECWCS'.
12 Sep 2019: Added CDP8 version of MRS distortion models
             while keeping CDP7 versions the default.
07 Oct 2019: Removed '.yaml' suffix from schema references.
26 Mar 2020: Ensure the model_type remains as originally defined when saving
             to a file.
             
@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""

# import warnings
import numpy as np
#import numpy.ma as ma

# Import the MIRI base data model and utilities.
from miri.datamodels.ancestry import get_my_model_type
from miri.datamodels.miri_model_base import MiriDataModel

# The distortion model might be represented by one of these STScI models,
# for example
#from jwst.datamodels.models import Poly2DModel, ICheb2DModel, ILegend2DModel, ...

# List all classes and global functions here.
__all__ = ['MiriImagingDistortionModel', 'MiriLrsD2WModel', \
           'MiriMrsDistortionModel12', 'MiriMrsDistortionModel34',
           'MiriMrsDistortionModel12_CDP6', 'MiriMrsDistortionModel34_CDP6',
           'MiriMrsDistortionModel12_CDP8', 'MiriMrsDistortionModel34_CDP8']


class MiriImagingDistortionModel(MiriDataModel):
    """
    
    A data model for MIRI distortion coefficients, based on the STScI
    base model, DataModel.
    
    After a data model has been created, data arrays and data tables
    are available as attributes with the same names as their input
    parameters, below. Metadata items are available within a ".meta"
    attribute tree.
    
    See https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html
    
    :Parameters:
    
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
        
        * None: A default data model with no shape. (If a data array is
          provided in the cmatrix parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
        
    bmatrix: numpy array (optional)
        An array containing the elements of the B matrix,
        describing distortion fit coefficients. Must be 2-D.
        A 3rd order polynomial fit will result in a 4x4 matrix.
        If a bmatrix parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    amatrix: numpy array (optional)
        An array containing the elements of the A matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 3rd order polynomial fit will result in a 4x4 matrix.
    tmatrix: numpy array (optional)
        An array containing the elements of the T matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    mmatrix: numpy array (optional)
        An array containing the elements of the M matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    bimatrix: numpy array (optional)
        An array containing the elements of the inverse B matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    aimatrix: numpy array (optional)
        An array containing the elements of the inverse A matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    timatrix: numpy array (optional)
        An array containing the elements of the inverse T matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    mimatrix: numpy array (optional)
        An array containing the elements of the inverse M matrix, 
        describing distortion fit coefficients. Must be 2-D.
        A 2nd order polynomial fit will result in a 3x3 matrix.
    boresight_offsets: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (parameter:object, filter:string,
        col_offset:number, row_offset:number).
        Or: A numpy record array containing the same information as above.
        If not specified, it will default to dummy values and no
        boresight_offset table will be assumed.        
    fitref: str (optional)
        A string containing a human-readable reference to a document
        describing the distortion model.
    fitmodel: str (optional)
        If a recognised JWST fitting model has been used (e.g. one of the
        models in the astropy.modeling package) a unique, machine-readable
        string defining the model used. If the model is not known or
        doesn't match one of the standard JWST models, leave this keyword
        blank and describe the model using the fitref parameter (above).
        Some example strings from astropy.modeling: Chebyshev1D', 'Chebyshev2D',
        'InverseSIP', 'Legendre1D','Legendre2D', 'Polynomial1D',
        'Polynomial2D', etc...
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
                
    """
    schema_url = "miri_distortion_imaging.schema"
    fieldnames = ('FILTER', 'COL_OFFSET', 'ROW_OFFSET')
    
    def __init__(self, init=None,
                 bmatrix=None, amatrix=None, tmatrix=None, mmatrix=None,
                 bimatrix=None, aimatrix=None, timatrix=None, mimatrix=None,
                 dmatrix=None, cmatrix=None, fmatrix=None, ematrix=None,
                 dimatrix=None, cimatrix=None, fimatrix=None, eimatrix=None,
                 fitref=None, fitmodel=None, boresight_offsets=None, **kwargs):
        """
        
        Initialises the MiriImagingDistortionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriImagingDistortionModel, self).__init__(init=init, **kwargs)

        # Data type is distortion map.
        self.meta.reftype = 'DISTORTION'
        # Initialise the model type
        self._init_data_type()
        # This is a reference data model.
        self._reference_model()

        # Verify the matrices have the correct shape. They are already
        # constrained to be 2-D in the schema.
        if bmatrix is not None:
            bmatrix = np.asarray(bmatrix)
            if bmatrix.ndim == 2:
                if bmatrix.shape[0] != bmatrix.shape[1]:
                    strg = "B Matrix should be square: "
                    strg += "%dx%d matrix provided instead." % bmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "B matrix should be 2-D. %d-D array provided." % \
                    bmatrix.ndim
                raise TypeError(strg)
            self.bmatrix = bmatrix

        if amatrix is not None:
            amatrix = np.asarray(amatrix)
            if amatrix.ndim == 2:
                if amatrix.shape[0] != amatrix.shape[1]:
                    strg = "A matrix should be square: "
                    strg += "%dx%d matrix provided instead." % amatrix.shape
                    raise TypeError(strg)
            else:
                strg = "A matrix should be 2-D. %d-D array provided." % \
                    self.amatrix.ndim
                raise TypeError(strg)
            self.amatrix = amatrix
            
        if tmatrix is not None:
            tmatrix = np.asarray(tmatrix)
            if tmatrix.ndim == 2:
                if tmatrix.shape[0] != tmatrix.shape[1]:
                    strg = "T matrix should be square: "
                    strg += "%dx%d matrix provided instead." % tmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "T matrix should be 2-D. %d-D array provided." % \
                    self.tmatrix.ndim
                raise TypeError(strg)
            self.tmatrix = tmatrix
            
        if mmatrix is not None:
            mmatrix = np.asarray(mmatrix)
            if mmatrix.ndim == 2:
                if mmatrix.shape[0] != mmatrix.shape[1]:
                    strg = "M matrix should be square: "
                    strg += "%dx%d matrix provided instead." % mmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "M matrix should be 2-D. %d-D array provided." % \
                    self.mmatrix.ndim
                raise TypeError(strg)
            self.mmatrix = mmatrix

        if dmatrix is not None:
            dmatrix = np.asarray(dmatrix)
            if dmatrix.ndim == 2:
                if dmatrix.shape[0] != dmatrix.shape[1]:
                    strg = "D Matrix should be square: "
                    strg += "%dx%d matrix provided instead." % dmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "D matrix should be 2-D. %d-D array provided." % \
                    dmatrix.ndim
                raise TypeError(strg)
            self.dmatrix = dmatrix

        if cmatrix is not None:
            cmatrix = np.asarray(cmatrix)
            if cmatrix.ndim == 2:
                if cmatrix.shape[0] != cmatrix.shape[1]:
                    strg = "C Matrix should be square: "
                    strg += "%dx%d matrix provided instead." % cmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "C matrix should be 2-D. %d-D array provided." % \
                    cmatrix.ndim
                raise TypeError(strg)
            self.cmatrix = cmatrix

        if fmatrix is not None:
            fmatrix = np.asarray(fmatrix)
            if fmatrix.ndim == 2:
                if fmatrix.shape[0] != fmatrix.shape[1]:
                    strg = "F Matrix should be square: "
                    strg += "%dx%d matrix provided instead." % fmatrix.shape
                    raise TypeError(strg)
            else:
                strg = "F matrix should be 2-D. %d-D array provided." % \
                    fmatrix.ndim
                raise TypeError(strg)
            self.fmatrix = fmatrix     

        if ematrix is not None:
            ematrix = np.asarray(ematrix)
            if ematrix.ndim == 2:
                if ematrix.shape[0] != ematrix.shape[1]:
                    strg = "E Matrix should be square: "
                    strg += "%dx%d matrix provided instead." % ematrix.shape
                    raise TypeError(strg)
            else:
                strg = "E matrix should be 2-D. %d-D array provided." % \
                    ematrix.ndim
                raise TypeError(strg)
            self.ematrix = ematrix   
       
        if bimatrix is not None:
            bimatrix = np.asarray(bimatrix)
            if bimatrix.ndim == 2:
                if bimatrix.shape[0] != bimatrix.shape[1]:
                    strg = "BI matrix should be square: "
                    strg += "%dx%d matrix provided instead." % bimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "BI matrix should be 2-D. %d-D array provided." % \
                    self.bimatrix.ndim
                raise TypeError(strg)
            self.bimatrix = bimatrix

        if aimatrix is not None:
            aimatrix = np.asarray(aimatrix)
            if aimatrix.ndim == 2:
                if aimatrix.shape[0] != aimatrix.shape[1]:
                    strg = "AI matrix should be square: "
                    strg += "%dx%d matrix provided instead." % aimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "AI matrix should be 2-D. %d-D array provided." % \
                    self.aimatrix.ndim
                raise TypeError(strg)
            self.aimatrix = aimatrix

        if timatrix is not None:
            timatrix = np.asarray(timatrix)
            if timatrix.ndim == 2:
                if timatrix.shape[0] != timatrix.shape[1]:
                    strg = "TI matrix should be square: "
                    strg += "%dx%d matrix provided instead." % timatrix.shape
                    raise TypeError(strg)
            else:
                strg = "TI matrix should be 2-D. %d-D array provided." % \
                    self.timatrix.ndim
                raise TypeError(strg)
            self.timatrix = timatrix

        if mimatrix is not None:
            mimatrix = np.asarray(mimatrix)
            if mimatrix.ndim == 2:
                if mimatrix.shape[0] != mimatrix.shape[1]:
                    strg = "MI matrix should be square: "
                    strg += "%dx%d matrix provided instead." % mimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "MI matrix should be 2-D. %d-D array provided." % \
                    self.mimatrix.ndim
                raise TypeError(strg)
            self.mimatrix = mimatrix

        if dimatrix is not None:
            dimatrix = np.asarray(dimatrix)
            if dimatrix.ndim == 2:
                if dimatrix.shape[0] != dimatrix.shape[1]:
                    strg = "DI Matrix should be square: "
                    strg += "%dx%d imatrix provided instead." % dimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "DI imatrix should be 2-D. %d-D array provided." % \
                    dimatrix.ndim
                raise TypeError(strg)
            self.dimatrix = dimatrix

        if cimatrix is not None:
            cimatrix = np.asarray(cimatrix)
            if cimatrix.ndim == 2:
                if cimatrix.shape[0] != cimatrix.shape[1]:
                    strg = "CI Matrix should be square: "
                    strg += "%dx%d imatrix provided instead." % cimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "CI imatrix should be 2-D. %d-D array provided." % \
                    cimatrix.ndim
                raise TypeError(strg)
            self.cimatrix = cimatrix

        if fimatrix is not None:
            fimatrix = np.asarray(fimatrix)
            if fimatrix.ndim == 2:
                if fimatrix.shape[0] != fimatrix.shape[1]:
                    strg = "FI Matrix should be square: "
                    strg += "%dx%d imatrix provided instead." % fimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "FI imatrix should be 2-D. %d-D array provided." % \
                    fimatrix.ndim
                raise TypeError(strg)
            self.fimatrix = fimatrix     

        if eimatrix is not None:
            eimatrix = np.asarray(eimatrix)
            if eimatrix.ndim == 2:
                if eimatrix.shape[0] != eimatrix.shape[1]:
                    strg = "EI Matrix should be square: "
                    strg += "%dx%d imatrix provided instead." % eimatrix.shape
                    raise TypeError(strg)
            else:
                strg = "EI imatrix should be 2-D. %d-D array provided." % \
                    eimatrix.ndim
                raise TypeError(strg)
            self.eimatrix = eimatrix   

        if boresight_offsets is not None:
            try:
                self.boresight_offsets = boresight_offsets
            except (ValueError, TypeError) as e:
                strg = "boresight_offsets must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Copy the units of the these arrays from the schema, if defined.
        aunits = self.set_data_units('amatrix')
        bunits = self.set_data_units('bmatrix')
        tunits = self.set_data_units('tmatrix')
        munits = self.set_data_units('mmatrix')
        dunits = self.set_data_units('dmatrix')
        cunits = self.set_data_units('cmatrix')
        funits = self.set_data_units('fmatrix')
        eunits = self.set_data_units('ematrix')
        biunits = self.set_data_units('bimatrix')
        aiunits = self.set_data_units('aimatrix')
        tiunits = self.set_data_units('timatrix')
        miunits = self.set_data_units('mimatrix')
        
        # Copy the table column units from the schema, if defined.
        boresight_units = self.set_table_units('boresight_offsets')
           
        if fitref is not None:
            self.meta.fit.reference = fitref
        if fitmodel is not None:
            self.meta.fit.model = fitmodel

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriImagingDistortionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()
           
    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "bmatrix"
        and not "data".
        
        """
        return 'bmatrix'

    def __str__(self):
        """
        
        Return the contents of the distortion map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        if self.meta.fit.model is not None:
            strg += "Fit model is \'%s\'\n" % str(self.meta.fit.model)
        if self.meta.fit.reference is not None:
            strg += "See \'%s\' for a description of the fit.\n" % \
                str(self.meta.fit.reference)
        strg += self.get_history_str()
            
        strg += self.get_data_str('bmatrix', underline=True, underchar="-")
        strg += self.get_data_str('amatrix', underline=True, underchar="-")
        strg += self.get_data_str('tmatrix', underline=True, underchar="-")
        strg += self.get_data_str('mmatrix', underline=True, underchar="-")
        strg += self.get_data_str('dmatrix', underline=True, underchar="-")
        strg += self.get_data_str('cmatrix', underline=True, underchar="-")
        strg += self.get_data_str('fmatrix', underline=True, underchar="-")
        strg += self.get_data_str('ematrix', underline=True, underchar="-")
        strg += self.get_data_str('bimatrix', underline=True, underchar="-")
        strg += self.get_data_str('aimatrix', underline=True, underchar="-")
        strg += self.get_data_str('timatrix', underline=True, underchar="-")
        strg += self.get_data_str('mimatrix', underline=True, underchar="-")
        strg += self.get_data_str('dimatrix', underline=True, underchar="-")
        strg += self.get_data_str('cimatrix', underline=True, underchar="-")
        strg += self.get_data_str('fimatrix', underline=True, underchar="-")
        strg += self.get_data_str('eimatrix', underline=True, underchar="-")        

        if self.boresight_offsets is not None:
            strg += self.get_data_str('boresight_offsets', underline=True, underchar="-")
        else:
            strg += "No boresight_offsets."
        
        return strg


class MiriLrsD2WModel(MiriDataModel):
    """
    
    A generic data model for a MIRI LRS distortion and wavelength
    calibration table.

    See MIRI-TR-10020-MPI for a detailed description of the content
    of the data model.

    After a data model has been created, the wavelength table is available
    within the attribute .wavelength_table. Metadata items are available
    within a ".meta" attribute tree.
    
    See https://jwst-pipeline.readthedocs.io/en/latest/jwst/datamodels/index.html
   
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
        
    wavelength_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (parameter:object, factor:number,
        uncertainty:number), giving the wavelength calibration factors valid for
        different generic parameters.
        Or: A numpy record array containing the same information as above.
        A wavelength table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_distortion_lrs.schema"
    fieldnames = ('X_CENTER', 'Y_CENTER', 'WAVELENGTH', 'X0', 'Y0', 'X1', 'Y1', \
                  'X2', 'Y2', 'X3', 'Y3')
    
    def __init__(self, init=None, wavelength_table=None, **kwargs):
        """
        
        Initialises the MiriLrsD2WModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsD2WModel, self).__init__(init=init, **kwargs)

        # Data type is wavelength calibration world coordinates.
        #self.meta.reftype = 'DISTORTION'
        self.meta.reftype = 'SPECWCS'
        # Initialise the model type
        self._init_data_type()       
        # This is a reference data model.
        self._reference_model()

        if wavelength_table is not None:
            try:
                self.wavelength_table = wavelength_table
            except (ValueError, TypeError) as e:
                strg = "wavelength_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # Copy the table column units from the schema, if defined.
        wavelength_units = self.set_table_units('wavelength_table')
        
        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriImagingDistortionModel, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

# TODO: Over-complicated data structure needs to be simplified.
class MiriMrsDistortionModel12(MiriDataModel):
    """
    
    A data model for a MIRI MRS distortion model - CHANNEL 34 VARIANT,
    based on the STScI base model, DataModel. Old CDP-7 version.
     
    :Parameters:
     
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
         
        * None: A default data model with no shape. (If a data array is
          provided in the lambda parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
         
    slicenumber: numpy array (optional)
        An array containing the elements of the slice array, which
        describes the mapping of pixel corners to slice number.
        Must be 2-D.
    fov_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    fov_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    albe_xanyan: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    xanyan_albe: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.
    bzero1: float (optional)
        Beta coordinate of the centre of slice 1 of channel 1
    bdel1: float (optional)
        Slice width (delta beta) for channel 1
    bzero2: float (optional)
        Beta coordinate of the centre of slice 1 of channel 2
    bdel2: float (optional)
        Slice width (delta beta) for channel 2
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
            
    """
    schema_url = "miri_distortion_mrs12.schema"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))
    
    def __init__(self, init=None, slicenumber=None, fov_ch1=None, fov_ch2=None,
                 alpha_ch1=None, lambda_ch1=None, alpha_ch2=None, lambda_ch2=None,
                 x_ch1=None, y_ch1=None, x_ch2=None, y_ch2=None,
                 albe_xanyan=None, xanyan_albe=None, bzero1=None, bdel1=None,
                 bzero2=None, bdel2=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel12 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel12, self).__init__(init=init, **kwargs)

        # Data type is MRS DISTORTION.
        self.meta.reftype = 'DISTORTION'
        # Initialise the model type
        self._init_data_type()      
        # This is a reference data model.
        self._reference_model()

        if slicenumber is not None:
            self.slicenumber = slicenumber
 
        # Define the beta coordinates and slice widths, if given
        if bzero1 is not None:
            self.meta.instrument.bzero1 = bzero1
        if bdel1 is not None:
            self.meta.instrument.bdel1 = bdel1
        if bzero2 is not None:
            self.meta.instrument.bzero2 = bzero2
        if bdel2 is not None:
            self.meta.instrument.bdel2 = bdel2
 
        if fov_ch1 is not None:
            try:
                self.fov_ch1 = fov_ch1
            except (ValueError, TypeError) as e:
                strg = "fov_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if fov_ch2 is not None:
            try:
                self.fov_ch2 = fov_ch2
            except (ValueError, TypeError) as e:
                strg = "fov_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        if alpha_ch1 is not None:
            try:
                self.alpha_ch1 = alpha_ch1
            except (ValueError, TypeError) as e:
                strg = "alpha_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch1 is not None:
            try:
                self.lambda_ch1 = lambda_ch1
            except (ValueError, TypeError) as e:
                strg = "lambda_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if alpha_ch2 is not None:
            try:
                self.alpha_ch2 = alpha_ch2
            except (ValueError, TypeError) as e:
                strg = "alpha_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch2 is not None:
            try:
                self.lambda_ch2 = lambda_ch2
            except (ValueError, TypeError) as e:
                strg = "lambda_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
          
        if x_ch1 is not None:
            try:
                self.x_ch1 = x_ch1
            except (ValueError, TypeError) as e:
                strg = "x_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch1 is not None:
            try:
                self.y_ch1 = y_ch1
            except (ValueError, TypeError) as e:
                strg = "y_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if x_ch2 is not None:
            try:
                self.x_ch2 = x_ch2
            except (ValueError, TypeError) as e:
                strg = "x_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch2 is not None:
            try:
                self.y_ch2 = y_ch2
            except (ValueError, TypeError) as e:
                strg = "y_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if albe_xanyan is not None:
            try:
                self.albe_to_xanyan = albe_xanyan
            except (ValueError, TypeError) as e:
                strg = "albe_xanyan must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if xanyan_albe is not None:
            try:
                self.xanyan_to_albe = xanyan_albe
            except (ValueError, TypeError) as e:
                strg = "xanyan_albe must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Copy the table column units from the schema, if defined.
        fov_ch1_units = self.set_table_units('fov_ch1')
        fov_ch2_units = self.set_table_units('fov_ch2')
        alpha_ch1_units = self.set_table_units('alpha_ch1')
        alpha_ch2_units = self.set_table_units('alpha_ch2')
        lambda_ch1_units = self.set_table_units('lambda_ch1')
        lambda_ch2_units = self.set_table_units('lambda_ch2')
        x_ch1_units = self.set_table_units('x_ch1')
        x_ch2_units = self.set_table_units('x_ch2')
        y_ch1_units = self.set_table_units('y_ch1')
        y_ch2_units = self.set_table_units('y_ch2')
        albe_to_xanyan_units = self.set_table_units('albe_to_xanyan')
        xanyan_to_albe_units = self.set_table_units('xanyan_to_albe')

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "slicenumber"
        and not "data".
        
        """
        return 'slicenumber'

    def __str__(self):
        """
        
        Return the contents of the D2C map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        strg += self.get_history_str()
            
        strg += self.get_data_str('slicenumber', underline=True, underchar="-")

        strg += self.get_data_str('fov_ch1', underline=True, underchar="-")
        strg += self.get_data_str('fov_ch2', underline=True, underchar="-")
 
        strg += self.get_data_str('alpha_ch1', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch1', underline=True, underchar="-")
        strg += self.get_data_str('alpha_ch2', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch2', underline=True, underchar="-")
  
        strg += self.get_data_str('x_ch1', underline=True, underchar="-")
        strg += self.get_data_str('y_ch1', underline=True, underchar="-")
        strg += self.get_data_str('x_ch2', underline=True, underchar="-")
        strg += self.get_data_str('y_ch2', underline=True, underchar="-")
 
        strg += self.get_data_str('albe_to_xanyan', underline=True, underchar="-")
        strg += self.get_data_str('xanyan_to_albe', underline=True, underchar="-")
        return strg

class MiriMrsDistortionModel12_CDP8(MiriDataModel):
    """
    
    A data model for a MIRI MRS distortion model - CHANNEL 12 VARIANT,
    based on the STScI base model, DataModel.
     
    :Parameters:
     
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
         
        * None: A default data model with no shape. (If a data array is
          provided in the lambda parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
         
    slicenumber: numpy array (optional)
        An array containing the elements of the slice array, which
        describes the mapping of pixel corners to slice number.
        Must be 2-D.
    fov_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    fov_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch1: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch2: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    albe_v2v3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    v2v3_albe: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.
    bzero1: float (optional)
        Beta coordinate of the centre of slice 1 of channel 1
    bdel1: float (optional)
        Slice width (delta beta) for channel 1
    bzero2: float (optional)
        Beta coordinate of the centre of slice 1 of channel 2
    bdel2: float (optional)
        Slice width (delta beta) for channel 2
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
            
    """
    schema_url = "miri_distortion_mrs12_CDP8.schema"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))
    
    def __init__(self, init=None, slicenumber=None, fov_ch1=None, fov_ch2=None,
                 alpha_ch1=None, lambda_ch1=None, alpha_ch2=None, lambda_ch2=None,
                 x_ch1=None, y_ch1=None, x_ch2=None, y_ch2=None,
                 albe_v2v3=None, v2v3_albe=None, bzero1=None, bdel1=None,
                 bzero2=None, bdel2=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel12 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel12_CDP8, self).__init__(init=init, **kwargs)

        # Data type is MRS DISTORTION.
        self.meta.reftype = 'DISTORTION'
        # Initialise the model type
        self._init_data_type()     
        # This is a reference data model.
        self._reference_model()

        if slicenumber is not None:
            self.slicenumber = slicenumber
 
        # Define the beta coordinates and slice widths, if given
        if bzero1 is not None:
            self.meta.instrument.bzero1 = bzero1
        if bdel1 is not None:
            self.meta.instrument.bdel1 = bdel1
        if bzero2 is not None:
            self.meta.instrument.bzero2 = bzero2
        if bdel2 is not None:
            self.meta.instrument.bdel2 = bdel2
 
        if fov_ch1 is not None:
            try:
                self.fov_ch1 = fov_ch1
            except (ValueError, TypeError) as e:
                strg = "fov_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if fov_ch2 is not None:
            try:
                self.fov_ch2 = fov_ch2
            except (ValueError, TypeError) as e:
                strg = "fov_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        if alpha_ch1 is not None:
            try:
                self.alpha_ch1 = alpha_ch1
            except (ValueError, TypeError) as e:
                strg = "alpha_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch1 is not None:
            try:
                self.lambda_ch1 = lambda_ch1
            except (ValueError, TypeError) as e:
                strg = "lambda_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if alpha_ch2 is not None:
            try:
                self.alpha_ch2 = alpha_ch2
            except (ValueError, TypeError) as e:
                strg = "alpha_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch2 is not None:
            try:
                self.lambda_ch2 = lambda_ch2
            except (ValueError, TypeError) as e:
                strg = "lambda_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
          
        if x_ch1 is not None:
            try:
                self.x_ch1 = x_ch1
            except (ValueError, TypeError) as e:
                strg = "x_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch1 is not None:
            try:
                self.y_ch1 = y_ch1
            except (ValueError, TypeError) as e:
                strg = "y_ch1 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if x_ch2 is not None:
            try:
                self.x_ch2 = x_ch2
            except (ValueError, TypeError) as e:
                strg = "x_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch2 is not None:
            try:
                self.y_ch2 = y_ch2
            except (ValueError, TypeError) as e:
                strg = "y_ch2 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if albe_v2v3 is not None:
            try:
                self.albe_to_v2v3 = albe_v2v3
            except (ValueError, TypeError) as e:
                strg = "albe_v2v3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if v2v3_albe is not None:
            try:
                self.v2v3_to_albe = v2v3_albe
            except (ValueError, TypeError) as e:
                strg = "v2v3_albe must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        
        # Copy the table column units from the schema, if defined.
        fov_ch1_units = self.set_table_units('fov_ch1')
        fov_ch2_units = self.set_table_units('fov_ch2')
        alpha_ch1_units = self.set_table_units('alpha_ch1')
        alpha_ch2_units = self.set_table_units('alpha_ch2')
        lambda_ch1_units = self.set_table_units('lambda_ch1')
        lambda_ch2_units = self.set_table_units('lambda_ch2')
        x_ch1_units = self.set_table_units('x_ch1')
        x_ch2_units = self.set_table_units('x_ch2')
        y_ch1_units = self.set_table_units('y_ch1')
        y_ch2_units = self.set_table_units('y_ch2')
        albe_to_v2v3_units = self.set_table_units('albe_to_v2v3')
        v2v3_to_albe_units = self.set_table_units('v2v3_to_albe')

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriMrsDistortionModel12_CDP8, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "slicenumber"
        and not "data".
        
        """
        return 'slicenumber'

    def __str__(self):
        """
        
        Return the contents of the D2C map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title(underline=True, underchar="=") + "\n"
        strg += self.get_meta_str(underline=True, underchar='-')
        strg += self.get_history_str()
            
        strg += self.get_data_str('slicenumber', underline=True, underchar="-")

        strg += self.get_data_str('fov_ch1', underline=True, underchar="-")
        strg += self.get_data_str('fov_ch2', underline=True, underchar="-")
 
        strg += self.get_data_str('alpha_ch1', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch1', underline=True, underchar="-")
        strg += self.get_data_str('alpha_ch2', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch2', underline=True, underchar="-")
  
        strg += self.get_data_str('x_ch1', underline=True, underchar="-")
        strg += self.get_data_str('y_ch1', underline=True, underchar="-")
        strg += self.get_data_str('x_ch2', underline=True, underchar="-")
        strg += self.get_data_str('y_ch2', underline=True, underchar="-")
 
        strg += self.get_data_str('albe_to_xanyan', underline=True, underchar="-")
        strg += self.get_data_str('xanyan_to_albe', underline=True, underchar="-")
        return strg

class MiriMrsDistortionModel12_CDP6(MiriMrsDistortionModel12):
    """
    
    This class can be used to access the old CDP-6 version of the
    MiriMrsDistortionModel12 data model.
    
    See the MiriMrsDistortionModel12 class for full documentation.
    
    """
    schema_url = "miri_distortion_mrs12_CDP6.schema"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))

    def __init__(self, init=None, slicenumber=None, fov_ch1=None, fov_ch2=None,
                 alpha_ch1=None, lambda_ch1=None, alpha_ch2=None, lambda_ch2=None,
                 x_ch1=None, y_ch1=None, x_ch2=None, y_ch2=None,
                 albe_xanyan=None, xanyan_albe=None, bzero1=None, bdel1=None,
                 bzero2=None, bdel2=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel12_CDP6 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel12_CDP6, self).__init__(init=init,
                slicenumber=None, fov_ch1=fov_ch1, fov_ch2=fov_ch2,
                 alpha_ch1=alpha_ch1, lambda_ch1=lambda_ch1,
                 alpha_ch2=alpha_ch2, lambda_ch2=lambda_ch2,
                 x_ch1=x_ch1, y_ch1=y_ch1, x_ch2=x_ch2, y_ch2=y_ch2,
                 albe_xanyan=albe_xanyan, xanyan_albe=xanyan_albe,
                 bzero1=bzero1, bdel1=bdel1,
                 bzero2=bzero2, bdel2=bdel2, **kwargs)


# TODO: Over-complicated data structure needs to be simplified.
class MiriMrsDistortionModel34(MiriDataModel):
    """
    
    A data model for a MIRI MRS distortion model - CHANNEL 34 VARIANT,
    based on the STScI base model, DataModel. Old CDP-7 version.
     
    :Parameters:
     
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
         
        * None: A default data model with no shape. (If a data array is
          provided in the lambda parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
         
    slicenumber: numpy array (optional)
        An array containing the elements of the slice array, which
        describes the mapping of pixel corners to slice number.
        Must be 2-D.
    fov_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    fov_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    albe_xanyan: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    xanyan_albe: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.
    bzero3: float (optional)
        Beta coordinate of the centre of slice 1 of channel 3
    bdel3: float (optional)
        Slice width (delta beta) for channel 3
    bzero4: float (optional)
        Beta coordinate of the centre of slice 1 of channel 4
    bdel4: float (optional)
        Slice width (delta beta) for channel 4
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
            
    """
    schema_url = "miri_distortion_mrs34.schema"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))
    
    def __init__(self, init=None, slicenumber=None, fov_ch3=None, fov_ch4=None,
                 alpha_ch3=None, lambda_ch3=None, alpha_ch4=None, lambda_ch4=None,
                 x_ch3=None, y_ch3=None, x_ch4=None, y_ch4=None,
                 albe_xanyan=None, xanyan_albe=None, bzero3=None, bdel3=None,
                 bzero4=None, bdel4=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel34 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel34, self).__init__(init=init, **kwargs)

        # Data type is MRS DISTORTION.
        self.meta.reftype = 'DISTORTION'
        # Initialise the model type
        self._init_data_type()     
        # This is a reference data model.
        self._reference_model()

        if slicenumber is not None:
            self.slicenumber = slicenumber
 
        # Define the beta coordinates and slice widths, if given
        if bzero3 is not None:
            self.meta.instrument.bzero3 = bzero3
        if bdel3 is not None:
            self.meta.instrument.bdel3 = bdel3
        if bzero4 is not None:
            self.meta.instrument.bzero4 = bzero4
        if bdel4 is not None:
            self.meta.instrument.bdel4 = bdel4
 
        if fov_ch3 is not None:
            try:
                self.fov_ch3 = fov_ch3
            except (ValueError, TypeError) as e:
                strg = "fov_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if fov_ch4 is not None:
            try:
                self.fov_ch4 = fov_ch4
            except (ValueError, TypeError) as e:
                strg = "fov_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        if alpha_ch3 is not None:
            try:
                self.alpha_ch3 = alpha_ch3
            except (ValueError, TypeError) as e:
                strg = "alpha_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch3 is not None:
            try:
                self.lambda_ch3 = lambda_ch3
            except (ValueError, TypeError) as e:
                strg = "lambda_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if alpha_ch4 is not None:
            try:
                self.alpha_ch4 = alpha_ch4
            except (ValueError, TypeError) as e:
                strg = "alpha_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch4 is not None:
            try:
                self.lambda_ch4 = lambda_ch4
            except (ValueError, TypeError) as e:
                strg = "lambda_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
          
        if x_ch3 is not None:
            try:
                self.x_ch3 = x_ch3
            except (ValueError, TypeError) as e:
                strg = "x_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch3 is not None:
            try:
                self.y_ch3 = y_ch3
            except (ValueError, TypeError) as e:
                strg = "y_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if x_ch4 is not None:
            try:
                self.x_ch4 = x_ch4
            except (ValueError, TypeError) as e:
                strg = "x_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch4 is not None:
            try:
                self.y_ch4 = y_ch4
            except (ValueError, TypeError) as e:
                strg = "y_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if albe_xanyan is not None:
            try:
                self.albe_to_xanyan = albe_xanyan
            except (ValueError, TypeError) as e:
                strg = "albe_xanyan must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if xanyan_albe is not None:
            try:
                self.xanyan_to_albe = xanyan_albe
            except (ValueError, TypeError) as e:
                strg = "xanyan_albe must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # Copy the table column units from the schema, if defined.
        fov_ch3_units = self.set_table_units('fov_ch3')
        fov_ch4_units = self.set_table_units('fov_ch4')
        alpha_ch3_units = self.set_table_units('alpha_ch3')
        alpha_ch4_units = self.set_table_units('alpha_ch4')
        lambda_ch3_units = self.set_table_units('lambda_ch3')
        lambda_ch4_units = self.set_table_units('lambda_ch4')
        x_ch3_units = self.set_table_units('x_ch3')
        x_ch4_units = self.set_table_units('x_ch4')
        y_ch3_units = self.set_table_units('y_ch3')
        y_ch4_units = self.set_table_units('y_ch4')
        albe_to_xanyan_units = self.set_table_units('albe_to_xanyan')
        xanyan_to_albe_units = self.set_table_units('xanyan_to_albe')

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "slicenumber"
        and not "data".
        
        """
        return 'slicenumber'

    def __str__(self):
        """
        
        Return the contents of the D2C map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()
            
        strg += self.get_data_str('slicenumber', underline=True, underchar="-")

        strg += self.get_data_str('fov_ch3', underline=True, underchar="-")
        strg += self.get_data_str('fov_ch4', underline=True, underchar="-")
 
        strg += self.get_data_str('alpha_ch3', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch3', underline=True, underchar="-")
        strg += self.get_data_str('alpha_ch4', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch4', underline=True, underchar="-")
  
        strg += self.get_data_str('x_ch3', underline=True, underchar="-")
        strg += self.get_data_str('y_ch3', underline=True, underchar="-")
        strg += self.get_data_str('x_ch4', underline=True, underchar="-")
        strg += self.get_data_str('y_ch4', underline=True, underchar="-")
 
        strg += self.get_data_str('albe_to_xanyan', underline=True, underchar="-")
        strg += self.get_data_str('xanyan_to_albe', underline=True, underchar="-")
        return strg

class MiriMrsDistortionModel34_CDP8(MiriDataModel):
    """
    
    A data model for a MIRI MRS distortion model - CHANNEL 34 VARIANT,
    based on the STScI base model, DataModel.
     
    :Parameters:
     
    init: shape tuple, file path, file object, pyfits.HDUList, numpy array
        An optional initializer for the data model, which can have one
        of the following forms:
         
        * None: A default data model with no shape. (If a data array is
          provided in the lambda parameter, the shape is derived from
          the array.)
        * Shape tuple: Initialize with empty data of the given shape.
        * File path: Initialize from the given file.
        * Readable file object: Initialize from the given file object.
        * pyfits.HDUList: Initialize from the given pyfits.HDUList.
         
    slicenumber: numpy array (optional)
        An array containing the elements of the slice array, which
        describes the mapping of pixel corners to slice number.
        Must be 2-D.
    fov_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    fov_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (alpha_min:value, beta_min:value)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    alpha_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    lambda_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    x_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    y_ch4: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    albe_v2v3: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.        
    v2v3_albe: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (...)
        Or: A numpy record array containing the same information as above.
        If not specified, no table will be defined.
    bzero3: float (optional)
        Beta coordinate of the centre of slice 1 of channel 3
    bdel3: float (optional)
        Slice width (delta beta) for channel 3
    bzero4: float (optional)
        Beta coordinate of the centre of slice 1 of channel 4
    bdel4: float (optional)
        Slice width (delta beta) for channel 4
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
            
    """
    schema_url = "miri_distortion_mrs34_CDP8.schema"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))
    
    def __init__(self, init=None, slicenumber=None, fov_ch3=None, fov_ch4=None,
                 alpha_ch3=None, lambda_ch3=None, alpha_ch4=None, lambda_ch4=None,
                 x_ch3=None, y_ch3=None, x_ch4=None, y_ch4=None,
                 albe_v2v3=None, v2v3_albe=None, bzero3=None, bdel3=None,
                 bzero4=None, bdel4=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel34 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel34_CDP8, self).__init__(init=init, **kwargs)

        # Data type is MRS DISTORTION.
        self.meta.reftype = 'DISTORTION'
        # Initialise the model type
        self._init_data_type()      
        # This is a reference data model.
        self._reference_model()

        if slicenumber is not None:
            self.slicenumber = slicenumber
 
        # Define the beta coordinates and slice widths, if given
        if bzero3 is not None:
            self.meta.instrument.bzero3 = bzero3
        if bdel3 is not None:
            self.meta.instrument.bdel3 = bdel3
        if bzero4 is not None:
            self.meta.instrument.bzero4 = bzero4
        if bdel4 is not None:
            self.meta.instrument.bdel4 = bdel4
 
        if fov_ch3 is not None:
            try:
                self.fov_ch3 = fov_ch3
            except (ValueError, TypeError) as e:
                strg = "fov_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if fov_ch4 is not None:
            try:
                self.fov_ch4 = fov_ch4
            except (ValueError, TypeError) as e:
                strg = "fov_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
         
        if alpha_ch3 is not None:
            try:
                self.alpha_ch3 = alpha_ch3
            except (ValueError, TypeError) as e:
                strg = "alpha_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch3 is not None:
            try:
                self.lambda_ch3 = lambda_ch3
            except (ValueError, TypeError) as e:
                strg = "lambda_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if alpha_ch4 is not None:
            try:
                self.alpha_ch4 = alpha_ch4
            except (ValueError, TypeError) as e:
                strg = "alpha_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if lambda_ch4 is not None:
            try:
                self.lambda_ch4 = lambda_ch4
            except (ValueError, TypeError) as e:
                strg = "lambda_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
          
        if x_ch3 is not None:
            try:
                self.x_ch3 = x_ch3
            except (ValueError, TypeError) as e:
                strg = "x_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch3 is not None:
            try:
                self.y_ch3 = y_ch3
            except (ValueError, TypeError) as e:
                strg = "y_ch3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if x_ch4 is not None:
            try:
                self.x_ch4 = x_ch4
            except (ValueError, TypeError) as e:
                strg = "x_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if y_ch4 is not None:
            try:
                self.y_ch4 = y_ch4
            except (ValueError, TypeError) as e:
                strg = "y_ch4 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if albe_v2v3 is not None:
            try:
                self.albe_to_v2v3 = albe_v2v3
            except (ValueError, TypeError) as e:
                strg = "albe_v2v3 must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
        if v2v3_albe is not None:
            try:
                self.v2v3_to_albe = v2v3_albe
            except (ValueError, TypeError) as e:
                strg = "v2v3_albe must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)

        # Copy the table column units from the schema, if defined.
        fov_ch3_units = self.set_table_units('fov_ch3')
        fov_ch4_units = self.set_table_units('fov_ch4')
        alpha_ch3_units = self.set_table_units('alpha_ch3')
        alpha_ch4_units = self.set_table_units('alpha_ch4')
        lambda_ch3_units = self.set_table_units('lambda_ch3')
        lambda_ch4_units = self.set_table_units('lambda_ch4')
        x_ch3_units = self.set_table_units('x_ch3')
        x_ch4_units = self.set_table_units('x_ch4')
        y_ch3_units = self.set_table_units('y_ch3')
        y_ch4_units = self.set_table_units('y_ch4')
        albe_to_v2v3_units = self.set_table_units('albe_to_v2v3')
        v2v3_to_albe_units = self.set_table_units('v2v3_to_albe')

        # Define the exposure type (if not already contained in the data model)
        # NOTE: This will only define an exposure type when a valid detector
        # is defined in the metadata.
        if not self.meta.exposure.type:
            self.set_exposure_type()

    def _init_data_type(self):
        # Initialise the data model type
        model_type = get_my_model_type( self.__class__.__name__ )
        self.meta.model_type = model_type        

    def on_save(self, path):
       super(MiriMrsDistortionModel34_CDP8, self).on_save(path)
        # Re-initialise data type on save
       self._init_data_type()

    def get_primary_array_name(self):
        """
        
        Returns the name "primary" array for this model, which controls
        the size of other arrays that are implicitly created.
        For this data structure, the primary array's name is "slicenumber"
        and not "data".
        
        """
        return 'slicenumber'

    def __str__(self):
        """
        
        Return the contents of the D2C map object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()
            
        strg += self.get_data_str('slicenumber', underline=True, underchar="-")

        strg += self.get_data_str('fov_ch3', underline=True, underchar="-")
        strg += self.get_data_str('fov_ch4', underline=True, underchar="-")
 
        strg += self.get_data_str('alpha_ch3', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch3', underline=True, underchar="-")
        strg += self.get_data_str('alpha_ch4', underline=True, underchar="-")
        strg += self.get_data_str('lambda_ch4', underline=True, underchar="-")
  
        strg += self.get_data_str('x_ch3', underline=True, underchar="-")
        strg += self.get_data_str('y_ch3', underline=True, underchar="-")
        strg += self.get_data_str('x_ch4', underline=True, underchar="-")
        strg += self.get_data_str('y_ch4', underline=True, underchar="-")
 
        strg += self.get_data_str('albe_to_v2v3', underline=True, underchar="-")
        strg += self.get_data_str('v2v3_to_albe', underline=True, underchar="-")
        return strg

class MiriMrsDistortionModel34_CDP6(MiriMrsDistortionModel34):
    """
    
    This class can be used to access the old CDP-6 version of the
    MiriMrsDistortionModel34 data model.
    
    See the MiriMrsDistortionModel34 class for full documentation.
    
    """
    schema_url = "miri_distortion_mrs34_CDP6.schema.yaml"
    fieldnames_fov = ('alpha_min', 'alpha_max')
    fieldnames_d2c = ['VAR1']
    for i in (0,1,2,3,4):
        for j in (0,1,2,3,4):
            fieldnames_d2c.append('VAR2_%d_%d' % (i,j))
    fieldnames_trans = ['Label']
    for i in (0,1):
        for j in (0,1):
            fieldnames_trans.append('COEFF_%d_%d' % (i,j))
            
    def __init__(self, init=None, slicenumber=None, fov_ch3=None, fov_ch4=None,
                 alpha_ch3=None, lambda_ch3=None, alpha_ch4=None, lambda_ch4=None,
                 x_ch3=None, y_ch3=None, x_ch4=None, y_ch4=None,
                 albe_xanyan=None, xanyan_albe=None, bzero3=None, bdel3=None,
                 bzero4=None, bdel4=None, **kwargs):
        """
        
        Initialises the MiriMrsDistortionModel34_CDP6 class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsDistortionModel34_CDP6, self).__init__(init=init,
                 slicenumber=slicenumber, fov_ch3=fov_ch3, fov_ch4=fov_ch4,
                 alpha_ch3=alpha_ch3, lambda_ch3=lambda_ch3,
                 alpha_ch4=alpha_ch4, lambda_ch4=lambda_ch4,
                 x_ch3=x_ch3, y_ch3=y_ch3, x_ch4=x_ch4, y_ch4=y_ch4,
                 albe_xanyan=albe_xanyan, xanyan_albe=xanyan_albe,
                 bzero3=bzero3, bdel3=bdel3,
                 bzero4=bzero4, bdel4=bdel4, **kwargs)

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MIRI distortion models module.")

    PLOTTING = False
    SAVE_FILES = False
    
    print("Testing the MiriImagingDistortionModel class.")
    bmatrix = [[0.1,0.2,0.3,0.4],
               [0.5,0.6,0.7,0.8],
               [0.8,0.7,0.6,0.5],
               [0.4,0.3,0.2,0.1]
               ]
    amatrix = [[0.1,0.0,0.0,0.0],
               [0.0,0.1,0.0,0.0],
               [0.0,0.0,0.1,0.0],
               [0.0,0.0,0.0,0.1]
               ]
    tmatrix = [[0.1,0.0,0.0],
               [0.0,0.1,0.0],
               [0.0,0.0,0.1],
               ]
    mmatrix = [[0.1,0.0,0.0],
               [0.0,0.1,0.0],
               [0.0,0.0,0.1],
               ]
    boffsets = [('One', 0.1, 0.2),
                ('Two', 0.2, 0.3)]
    with MiriImagingDistortionModel( bmatrix=bmatrix, amatrix=amatrix,
                              tmatrix=tmatrix, mmatrix=mmatrix,
                              dmatrix=amatrix, cmatrix=bmatrix,
                              fmatrix=amatrix, ematrix=bmatrix,
                              boresight_offsets=boffsets,
                              fitref='MIRI-TN-00070-ATC version 3',
                              fitmodel='Polynomial2D' ) as testdata1:
        # This is how to set the matrix units (if not obtained from a file).
        testdata1.meta.bmatrix.units = 'mm ** (1-ij)'
        testdata1.meta.amatrix.units = 'mm ** (1-ij)'
        print(testdata1)
        if PLOTTING:
            testdata1.plot(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_imaging_distortion_model1.fits", overwrite=True)
        del testdata1
 
    print("Testing the MiriLrsD2WModel class.")
    wavedata = [
      (61.13976, 80.25328, 5.12652, 78.78318, 80.00265, 43.53634, 81.47993, 43.49634, 80.50392, 78.74317, 79.02664),
      (61.09973, 79.27727, 5.15676, 78.74317, 79.02664, 43.49634, 80.50392, 43.45628, 79.52791, 78.70311, 78.05063),
      (61.05965, 78.30126, 5.18700, 78.70311, 78.05063, 43.45628, 79.52791, 43.41617, 78.55190, 78.66300, 77.07462),
      (61.01951, 77.32526, 5.21725, 78.66300, 77.07462, 43.41617, 78.55190, 43.37601, 77.57589, 78.62285, 76.09861),
      (60.97933, 76.34925, 5.24749, 78.62285, 76.09861, 43.37601, 77.57589, 43.33580, 76.59988, 78.58264, 75.12260),
      (60.93910, 75.37324, 5.27773, 78.58264, 75.12260, 43.33580, 76.59988, 43.29554, 75.62387, 78.54238, 74.14659),
      (60.89881, 74.39723, 5.30797, 78.54238, 74.14659, 43.29554, 75.62387, 43.25523, 74.64786, 78.50207, 73.17058),
      (60.85848, 73.42122, 5.33821, 78.50207, 73.17058, 43.25523, 74.64786, 43.21487, 73.67186, 78.46171, 72.19458)
      ]
  
    print("\nWavelength calibration with table derived from list of tuples:")
    with MiriLrsD2WModel( wavelength_table=wavedata ) as testwave1:
        print(testwave1)
        if PLOTTING:
            testwave1.plot(description="testwave1")
        if SAVE_FILES:
            testwave1.save("test_lrs_d2w_model1.fits", overwrite=True)
        del testwave1
  
    print("Testing the MiriMrsDistortionModel classes.")    
    slicenumber = [[1,2,3,4],
                   [1,2,3,4],
                   [1,2,3,4],
                   [1,2,3,4]
                  ]
    slicenumber3 = [slicenumber, slicenumber]
    fovdata = [(-2.95, 3.09),
               (-2.96, 3.00)]
    d2cdata = [(100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                21.0, 22.0, 323.0, 24.0, 25.0), 
               (101.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                21.0, 22.0, 323.0, 24.0, 25.0)]
    c2ddata = [(99.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                21.0, 22.0, 323.0, 24.0, 25.0), 
               (98.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                21.0, 22.0, 323.0, 24.0, 25.0)]
    transform = [('T_CH3C,V2', 0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.71, 0.81, 0.91),
                 ('T_CH3C,V3', 0.12, 0.22, 0.32, 0.42, 0.52, 0.62, 0.72, 0.82, 0.92)]
     
    with MiriMrsDistortionModel12_CDP8( slicenumber=slicenumber3,
                                 fov_ch1=fovdata, fov_ch2=fovdata,
                                 alpha_ch1=d2cdata, lambda_ch1=d2cdata,
                                 alpha_ch2=d2cdata, lambda_ch2=d2cdata,
                                 x_ch1=c2ddata, y_ch1=c2ddata,
                                 x_ch2=c2ddata, y_ch2=c2ddata,
                                 albe_v2v3=transform, v2v3_albe=transform,
                                 bzero1=-1.772, bdel1=0.177,
                                 bzero2=-2.238, bdel2=0.280
                                 ) as testdata1:
        
        print(testdata1)
        print("Data arrays=", testdata1.list_data_arrays())
        print("Data tables=", testdata1.list_data_tables())
        if PLOTTING:
            testdata1.plot(description="testdata1")
        if SAVE_FILES:
            testdata1.save("test_mrs_distortion_model1.fits", overwrite=True)    
#             newmodel = MiriMrsDistortionModel12("test_mrs_distortion_model1.fits")
#             print(newmodel)
        del testdata1

    with MiriMrsDistortionModel34_CDP8( slicenumber=slicenumber3,
                                 fov_ch3=fovdata, fov_ch4=fovdata,
                                 alpha_ch3=d2cdata, lambda_ch3=d2cdata,
                                 alpha_ch4=d2cdata, lambda_ch4=d2cdata,
                                 x_ch3=c2ddata, y_ch3=c2ddata,
                                 x_ch4=c2ddata, y_ch4=c2ddata,
                                 albe_v2v3=transform, v2v3_albe=transform,
                                 bzero3=-1.772, bdel3=0.177,
                                 bzero4=-2.238, bdel4=0.280
                                 ) as testdata2:
        
        print(testdata2)
        print("Data arrays=", testdata2.list_data_arrays())
        print("Data tables=", testdata2.list_data_tables())
        if PLOTTING:
            testdata2.plot(description="testdata2")
        if SAVE_FILES:
            testdata2.save("test_mrs_distortion_model2.fits", overwrite=True)    
#             newmodel = MiriMrsDistortionModel34("test_mrs_distortion_model2.fits")
#             print(newmodel)
        del testdata2

    print("Test finished.")
