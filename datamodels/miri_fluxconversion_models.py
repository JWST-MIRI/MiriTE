#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Some extensions to the standard STScI data model which 
define MIRI flux conversion models.

NOTE: After CDP-4, the MIRI imaging mode data models have been
replaced by new photometric models contained in miri_photometric_models.py.
The old models are defined here for backwards compatibility with
CDPs released before CDP-4.

The MRS flux conversion model, MiriMrsFluxconversionModel, is still in use.

:Reference:

The STScI jwst.datamodels documentation. See
http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

:History:

10 Jan 2013: Created
21 Jan 2013: Included data type in metadata. Added plotting.
05 Feb 2013: Reformatted test code using "with" context manager.
             Modified to use functions from MiriDataModel.
08 Feb 2013: Replaced 'to_fits' with more generic 'save' method.
28 Feb 2013: Corrected mismatch in colour correction column names.
04 Mar 2013: Absolute spectral response keywords, REFSRF, REFERR and
             REFWAVEL, added to metadata.
17 May 2013: Do not allow a blank table to be created.
13 Aug 2013: Distinguish imaging and spectroscopy flux conversion models.
20 Aug 2013: Added FluxModel
22 Aug 2013: columnnames renamed to fieldnames.
02 Sep 2013: Pass the responsibility for creating record arrays to jwst.datamodels
             - a solution to the "Types in column 0 do not match" problem
             suggested by Michael Droettboom at STScI.
12 Sep 2013: Make the creation of an empty table a warning rather than an
             error.
30 Oct 2013: Renamed the schemas so they have adjacent listings. LRS and MRS
             now use different flux conversion model classes.
10 Dec 2013: Delimiter in MIRI schema names changed from "." to "_".
04 Mar 2014: Updated MiriLrsFluxconversionModel: changed from RSRF to SRF;
             removed keywords srf_at_reference, srf_error, srf_wavelength.
09 Jul 2014: field_def changed to dq_def.
29 Aug 2014: Included new reference file keywords (REFTYPE, AUTHOR, PEDIGREE)
25 Sep 2014: Updated the reference flags. insert_value_column function
             used to convert between 3 column and 4 column flag tables.
             TYPE and REFTYPE are no longer identical.
30 Sep 2014: Superflous flags commented out.
09 Dec 2014: Removed old commented out warnings.
30 Jun 2015: Refactored the JSON schemas so they are no longer hierarchical,
             which was causing the table column units to be lost, and
             removed FluxModel.
09 Jul 2015: Removed duplication of table units between schema and metadata.
             Units are now only defined in the metadata.
             Use of the fieldnames class variable removed from the code and
             deprecated. It is now used only by a few conversion scripts.
10 Sep 2015: Added a logging facility and a note that the imaging flux
             conversion model is deprecated.
11 Sep 2015: Removed duplicated plot method. Added plot_srf method.
03 Dec 2015: Added MiriPowerlawColourCorrectionModel.
10 Dec 2015: TYPE and REFTYPE strings rationalised.
31 Aug 2015: Python logger object removed from class.
31 Oct 2016: New data format for MiriMrsFluxConversionModel. PIXSIZ array added.
             The data arrays can be 2-D or 3-D (for backwards compatibility).
             Default flux units changed from Jy to mJy.
06 Dec 2016: Worked around data model problem with the pixsiz attribute.
15 Jun 2017: meta.reffile schema level removed to match changes in the
             JWST build 7.1 data models release. meta.reffile.type also
             changed to meta.reftype. TYPE keyword replaced by DATAMODL.
12 Jul 2017: Replaced "clobber" parameter with "overwrite".
10 Oct 2018: MiriLrsFluxconversionModel marked as deprecated. The JWST pipeline
             will use a combined MiriPhotometricModel after CDP-7, but
             MiriLrsFluxconversionModel will remain to support the existing
             LRS-only files. Corrected __all__.
14 Nov 2018: Added PHOTMJSR and PHOTUJA2 keywords to the schema.

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""
# This module is now converted to Python 3.


# Python logging facility.
import logging
logging.basicConfig(level=logging.INFO)  # Choose ERROR, WARN, INFO or DEBUG 
LOGGER = logging.getLogger("miri.datamodels") # Get a logger

import numpy as np

# Import the MIRI base data model and utilities.
from jwst.datamodels.model_base import DataModel
from miri.datamodels.dqflags import insert_value_column
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.miri_measured_model import MiriMeasuredModel

# List all classes and global functions here.
__all__ = ['MiriFluxconversionModel', \
           'MiriImagingFluxconversionModel', \
           'MiriImagingColourCorrectionModel', \
           'MiriPowerlawColourCorrectionModel', \
           'MiriLrsFluxconversionModel', \
           'MiriMrsFluxconversionModel']

# class FluxModel(DataModel):
#     """
#     
#     A data model for a trivial flux table. Used for testing.
#     
#     """
#     schema_url = "flux.schema.yaml"
#     fieldnames = ('parameter', 'factor', 'uncertainty')
# 
#     def __init__(self, init=None, flux_table=None, **kwargs):
#         super(FluxModel, self).__init__(init=init, **kwargs)
# 
#         if flux_table is not None:
#             try:
#                 self.flux_table = flux_table
#             except (ValueError, TypeError) as e:
#                 strg = "flux_table must be a numpy record array or list of records"
#                 strg += "\n   %s" % str(e)
#                 raise TypeError(strg)


class MiriFluxconversionModel(MiriDataModel):
    """
    
    A generic data model for a MIRI flux conversion table, based on the STScI
    base model, DataModel.
    
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
        
    flux_table: list of tuples or numpy record array (optional)
        Either: A list of tuples containing (parameter:object, factor:number,
        uncertainty:number), giving the flux conversion factors valid for
        different generic parameters.
        Or: A numpy record array containing the same information as above.
        A flux table must either be defined in the initializer or in
        this parameter. A blank table is not allowed.
    \*\*kwargs:
        All other keyword arguments are passed to the DataModel initialiser.
        See the jwst.datamodels documentation for the meaning of these keywords.
        
    """
    schema_url = "miri_fluxconversion.schema.yaml"
    fieldnames = ('PARAMETER', 'FACTOR', 'UNCERTAINTY')
    
    def __init__(self, init=None, flux_table=None, **kwargs):
        """
        
        Initialises the MiriFluxconversionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriFluxconversionModel, self).__init__(init=init, **kwargs)

        # Data type is flux conversion.
        self.meta.model_type = 'PHOTOM'
        self.meta.reftype = 'PHOTOM'
        
        # This is a reference data model.
        self._reference_model()

        if flux_table is not None:
            try:
                self.flux_table = flux_table
            except (ValueError, TypeError) as e:
                strg = "flux_table must be a numpy record array or list of records."
                strg += "\n   %s" % str(e)
                raise TypeError(strg)
# 
#         # Copy the table column units from the schema, if defined.
#         tableunits = self.set_table_units('flux_table')
        
    # TODO: Is this function needed?
    def __str__(self):
        """
        
        Return the contents of the flux conversion object as a readable
        string.
        
        """
        # Start with the data object title, metadata and history
        strg = self.get_title_and_metadata()

        # Describe the flux conversion table
        if self.flux_table is not None:
            strg += self.get_data_str('flux_table', underline=True, underchar="-")
        return strg


class MiriImagingFluxconversionModel(MiriFluxconversionModel):
    """
    
    A data model for MIRI imaging mode flux conversion data, based on
    MiriFluxconversionModel, but with the flux factors looked up
    by filter name.
    
    NOTE: OBSOLETE DATA MODEL. VALID FOR PRE-CDP-4 DATA ONLY.
    
    :Parameters:
    
    Exactly the same as MiriFluxconversionModel.
    
    """
    schema_url = "miri_fluxconversion_imaging.schema.yaml"
    fieldnames = ('FILTER', 'FACTOR', 'UNCERTAINTY')
    info_logged = False

    def __init__(self, init=None, flux_table=None, **kwargs):
        """
        
        Initialises the MiriImagingFluxconversionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriImagingFluxconversionModel, self).__init__(init=init,
                                                flux_table=flux_table,
                                                **kwargs)
        # Log an information string just once.
        if not MiriImagingFluxconversionModel.info_logged:
            LOGGER.info(self.__class__.__name__ + ": This class is valid for pre-CDP-4 data only")
            MiriImagingFluxconversionModel.info_logged = True
        # Data type is imaging flux conversion.
        self.meta.model_type = 'PHOTOM (OLD Imaging)'
        self.meta.reftype = 'PHOTOM'


class MiriImagingColourCorrectionModel(MiriFluxconversionModel):
    """
    
    A data model for MIRI imaging mode colour correction data, based on
    MiriFluxconversionModel, but with the flux factors looked up
    by black body temperature and filter name.
    
    :Parameters:
    
    Exactly the same as MiriFluxconversionModel.
    
    """
    schema_url = "miri_colourcorrection_imaging.schema.yaml"
    fieldnames = ('BBTEMP', 'FILTER', 'FACTOR', 'UNCERTAINTY')

    def __init__(self, init=None, flux_table=None, **kwargs):
        """
        
        Initialises the MiriImagingColourCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriImagingColourCorrectionModel, self).__init__(init=init,
                                                flux_table=flux_table,
                                                **kwargs)
        # Data type is colour correction.
        self.meta.model_type = 'COLCORR (Black Body)'
        self.meta.reftype = 'COLCORR'


class MiriPowerlawColourCorrectionModel(MiriFluxconversionModel):
    """
    
    A data model for MIRI imaging mode colour correction data, based on
    MiriFluxconversionModel, but with the flux factors looked up
    by black body temperature and filter name.
    
    :Parameters:
    
    Exactly the same as MiriFluxconversionModel.
    
    """
    schema_url = "miri_colourcorrection_powerlaw.schema.yaml"
    fieldnames = ('EXPONENT', 'FILTER', 'FACTOR', 'UNCERTAINTY')

    def __init__(self, init=None, flux_table=None, **kwargs):
        """
        
        Initialises the MiriImagingColourCorrectionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriPowerlawColourCorrectionModel, self).__init__(init=init,
                                                flux_table=flux_table,
                                                **kwargs)
        # Data type is colour correction.
        self.meta.model_type = 'COLCORR (Power Law)'
        self.meta.reftype = 'COLCORRPL'


class MiriLrsFluxconversionModel(MiriFluxconversionModel):
    """
    
    A data model for MIRI LRS mode flux conversion data, based on
    MiriFluxconversionModel, but with the flux factors looked up
    by wavelength.
    
    NOTE: DEPRECATED DATA MODEL. VALID FOR INDIVIDUAL PRE-CDP-7
    LRS DATA ONLY. After CDP-7, the imager and LRS data are merged
    into a single MiriPhotometricModel.
    
    :Parameters:
    
    Exactly the same as MiriFluxconversionModel.
    
    """
    schema_url = "miri_fluxconversion_lrs.schema.yaml"
    fieldnames = ('WAVELENGTH', 'SRF', 'SRF_ERROR')

    def __init__(self, init=None, flux_table=None, **kwargs):
        """
        
        Initialises the MiriLrsFluxconversionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriLrsFluxconversionModel, self).__init__(init=init,
                                                    flux_table=flux_table,
                                                    **kwargs)

        # Data type is LRS flux conversion.
        self.meta.model_type = 'PHOTOM (LRS)'
        self.meta.reftype = 'PHOTOM'

    def plot_srf(self, description=''):
        """
        
        Plot a graph showing the SRF as a function of wavelength.
        
        :Parameters:
        
        description: string, optional
            Additional description to be shown on the plot, if required. 
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # TODO: TO BE COMPLETED PROPERLY USING VISITOR CLASS
        
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt

        # Extract wavelength and SRF arrays from the flux_table.
        wavlist = []
        srflist = []
        errlist = []
        for (wavelength, srf, srf_error) in self.flux_table:
            wavlist.append(wavelength)
            srflist.append(srf)
            errlist.append(srf_error)

        # Create a new figure
        tstrg = self.__class__.__name__
        if description:
            tstrg += " - " + description
        fig = mplt.new_figure(1, stitle=tstrg)
        ax = mplt.add_subplot(fig, 1, 1, 1)

        # Plot SRF against wavelength
        xlabel = "Wavelength"
        ylabel = 'SRF'
        title = self.get_data_title('flux_table')
        # Plot each ramp as an XY plot overlaid on the same axis.
        mplt.plot_xy(wavlist, srflist, xlabel=xlabel, ylabel=ylabel,
                     title=title)
        mplt.show_plot()
        mplt.close()


# The new JWST SRF reference flags
mrssrf_reference_setup = \
            [(0, 'DO_NOT_USE',         'Bad pixel. Do not use.'),
             (1, 'NON_SCIENCE',        'Pixel not on science portion of detector'),
             (2, 'CDP_LOW_QUAL',       'Data of low quality')]
#              (3, 'CDP_DO_NOT_USE_ERR', 'Data without reliable error estimate')]
mrssrf_reference_flags = insert_value_column( mrssrf_reference_setup )

class MiriMrsFluxconversionModel(MiriMeasuredModel):
    """
    
    A data model for MIRI LRS mode flux conversion data, based on
    MiriFluxconversionModel, but with the flux factors looked up
    by wavelength.
    
    :Parameters:
    
    The same as MiriMeasuredModel plus
    
    pixsiz: numpy array (optional)
        An array containing the pixel size data.
        If a data parameter is provided, its contents overwrite the
        data initialized by the init parameter.
    srf_at_reference: number, optional
        The absolute spectral response at the reference wavelength
        [in (DN/s)/Jy].
    srf_error: number, optional
        The uncertainty in the absolute spectral response at the
        reference wavelength [in (DN/s)/Jy].
    srf_wavelength: number, optional
        The reference wavelength [in microns].
    
    """
    schema_url = "miri_fluxconversion_mrs.schema.yaml"
    _default_dq_def = mrssrf_reference_flags

    def __init__(self, init=None, data=None, dq=None, err=None, dq_def=None,
                 pixsiz=None,
                 srf_at_reference=None, srf_error=None, srf_wavelength=None,
                 **kwargs):
        """
        
        Initialises the MiriMrsFluxconversionModel class.
        
        Parameters: See class doc string.

        """
        super(MiriMrsFluxconversionModel, self).__init__(init=init, data=data,
                                                         dq=dq, err=err,
                                                         dq_def=dq_def,
                                                         **kwargs)

        # Data type is MRS flux conversion.
        self.meta.model_type = 'PHOTOM (MRS)'
        self.meta.reftype = 'PHOTOM'
        
        # This is a reference data model.
        self._reference_model()
        
        # Set the metadata, if provided.
        if srf_at_reference is not None:
            self.meta.srf_at_reference = srf_at_reference
        if srf_error is not None:
            self.meta.srf_error = srf_error
        if srf_wavelength is not None:
            self.meta.srf_wavelength = srf_wavelength
            
        if pixsiz is not None:
            # Work around the "list object does not have shape attribute"
            # data model problem. Why does the pixsiz attribute need to
            # be converted with np.asarray but the data and err attributes
            # do not?
            self.pixsiz = np.asarray(pixsiz)

#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miri/datamodels/tests.
#
if __name__ == '__main__':
    print("Testing the MiriFluxconversionModel module.")
    
    PLOTTING = False
    SAVE_FILES = False
    
    flux_im = [('F560W',  1.0,  0.0),
               ('F770W',  1.1,  0.0),
               ('F1000W', 1.2,  0.01),
               ('F1130W', 1.3,  0.0),
               ('F1280W', 1.4,  0.0),
               ('F1500W', 1.5,  0.02),
               ('F1800W', 1.6,  0.0),
               ('F2100W', 1.7,  0.03),
               ('F2550W', 1.8,  0.0),
               ]
    cc_im =   [(1.0, 'F560W',  1.0,  0.0),
               (1.0, 'F770W',  1.1,  0.0),
               (1.0, 'F1000W', 1.2,  0.01),
               (1.0, 'F1130W', 1.3,  0.0),
               (1.0, 'F1280W', 1.4,  0.0),
               (1.0, 'F1500W', 1.5,  0.02),
               (1.0, 'F1800W', 1.6,  0.0),
               (1.0, 'F2100W', 1.7,  0.03),
               (1.0, 'F2550W', 1.8,  0.0),
               (5.0, 'F560W',  1.0,  0.0),
               (5.0, 'F770W',  1.1,  0.0),
               (5.0, 'F1000W', 1.2,  0.01),
               (5.0, 'F1130W', 1.3,  0.0),
               (5.0, 'F1280W', 1.4,  0.0),
               (5.0, 'F1500W', 1.5,  0.02),
               (5.0, 'F1800W', 1.6,  0.0),
               (5.0, 'F2100W', 1.7,  0.03),
               (5.0, 'F2550W', 1.8,  0.0),
               ]
    cc_pl =   [(1.0, 'F560W',  1.0,  0.0),
               (1.1, 'F770W',  1.1,  0.0),
               (1.2, 'F1000W', 1.2,  0.01),
               (1.3, 'F1130W', 1.3,  0.0),
               (1.4, 'F1280W', 1.4,  0.0),
               (1.5, 'F1500W', 1.5,  0.02),
               (1.6, 'F1800W', 1.6,  0.0),
               (1.7, 'F2100W', 1.7,  0.03),
               (1.8, 'F2550W', 1.8,  0.0)
               ]
    flux_lrs = [( 2.0, 1.0,  0.0),
               ( 4.0, 1.1,  0.0),
               ( 6.0, 1.2,  0.01),
               ( 8.0, 1.3,  0.0),
               (10.0, 1.4,  0.0),
               (12.0, 1.5,  0.02),
               (14.0, 1.6,  0.0),
               (16.0, 1.7,  0.03),
               (18.0, 1.8,  0.0),
               ]
    flux_plane = [(1.0, 1.1, 1.2),
                  (1.3, 1.4, 1.5),
                  (1.6, 1.7, 1.8)]
    err_plane  = [(0.0, 0.01, 0.0),
                  (0.02, 0.0, 0.03),
                  (0.01, 0.04, 0.0)]
    pixsiz     = [(0.09, 0.11, 0.07),
                  (0.08, 0.10, 0.08),
                  (0.09, 0.07, 0.12)]
    flux_mrs = [flux_plane, flux_plane, flux_plane]
    err_mrs = [err_plane, err_plane, err_plane]

    print("\nImaging flux with factors derived from list of tuples:")
    with MiriImagingFluxconversionModel( flux_table=flux_im ) as testflux1:
        print(testflux1)
        if PLOTTING:
            testflux1.plot(description="testflux1")
        if SAVE_FILES:
            testflux1.save("test_im_flux_model1.fits", overwrite=True)
        del testflux1

    print("\nLRS flux with factors derived from list of tuples:")
    with MiriLrsFluxconversionModel( flux_table=flux_lrs ) as testflux3:
        print(testflux3)
        if PLOTTING:
            testflux3.plot(description="testflux3")
            testflux3.plot_srf(description="testflux3")
        if SAVE_FILES:
            testflux3.save("test_lrs_flux_model1.fits", overwrite=True)
        del testflux3

    print("\nMRS flux with RSRF defined in a data cube (old format):")
    with MiriMrsFluxconversionModel( data=flux_mrs, err=err_mrs ) as testflux3:
        print(testflux3)
        if PLOTTING:
            testflux3.plot(description="testflux3")
        if SAVE_FILES:
            testflux3.save("test_mrs_flux_model1.fits", overwrite=True)
        del testflux3

    print("\nMRS flux with RSRF defined in a the detector plane (new format):")
    with MiriMrsFluxconversionModel( data=flux_plane, err=err_plane,
                                     pixsiz=pixsiz ) as testflux4:
        print(testflux4)
        if PLOTTING:
            testflux4.plot(description="testflux4")
        if SAVE_FILES:
            testflux4.save("test_mrs_flux_model2.fits", overwrite=True)
        del testflux4

    print("\nImaging black body colour correction with factors derived from list of tuples:")
    with MiriImagingColourCorrectionModel( flux_table=cc_im ) as testflux5:
        print(testflux5)
        if PLOTTING:
            testflux5.plot(description="testflux1")
        if SAVE_FILES:
            testflux5.save("test_im_ccbb_model1.fits", overwrite=True)
        del testflux5

    print("\nImaging power law colour correction with factors derived from list of tuples:")
    with MiriPowerlawColourCorrectionModel( flux_table=cc_pl ) as testflux6:
        print(testflux6)
        if PLOTTING:
            testflux6.plot(description="testflux1")
        if SAVE_FILES:
            testflux6.save("test_im_ccpl_model1.fits", overwrite=True)
        del testflux6
        
    print("Test finished.")
