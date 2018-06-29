#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Import the CDP models. Adding this module to datamodels allows
external software to import all the MIRI CDP models with one import line,
for example

   from datamodels.cdp import MiriBadPixelMaskModel, MiriFlatfieldModel
   
instead of needing to know the names of all the individual source files.

:Reference:

MIRI Calibration Data Products Documentation:


JWST Calibration Data Products Documentation:
http://ssb.stsci.edu/doc/jwst/jwst/introduction.html#crds-reference-files

:History:

23 May 2013: Created
29 Jul 2013: CDP_DICT added.
13 Aug 2013: Distinguish imaging and spectroscopy flux conversion models.
02 Oct 2013: Added MiriJumpModel.
16 Oct 2013: Read a PSF as a MiriImagingPointSpreadFunctionModel.
30 Oct 2013: Repackaged the distortion classes into one module. RSRF and
             PSF models split into separate classes for IM/LRS/MRS.
11 Dec 2013: Modified the dictionary to distinguish CDPs by detector and
             filter as well as by data type.
11 Feb 2014: Modified to keep up with change made by schreibe on 10 Feb:
             LRSRSRF --> LRSSRF.
31 Mar 2014: Added MiriMrsStraylightModel and MiriLastFrameModel to CDP_DICT.
04 Jul 2014: Modified to assign LRS absolute flux calibration CDP to the 
             appropriate datamodel, based on filter.
16 Jul 2014: Added MiriGainModel and MiriReadnoiseModel to CDP_DICT.
21 Jul 2014: Detector names changed to MIRIMAGE, MIRIFUSHORT and MIRIFULONG.
             Both sets of names are included in CDP_DICT for backwards
             compatibility.
25 Sep 2014  Added MiriResetModel. Added new data types to CDP_DICT.
             NOTE: CDP_DICT is now getting very messy.
16 Oct 2014: REFTYPE of WCS changed to DISTORTION.
30 Oct 2014: Renamed REFTYPE = BAD to REFTYPE = MASK, for bad pixel mask. Added
             SKYFLAT to data types and associated it with MiriFlatfieldModel.
31 Oct 2014: Added MiriIPCModel. Added corresponding data type IPC to CDP_DICT.
25 Jun 2015: Updated the expected set of REFTYPEs and data products for
             CDP-4. TCORR and WSHIFT still need to be added.
26 Jun 2015: Added MiriMrsTransmissionCorrectionModel
30 Jun 2015: Added MiriWavelengthCorrectionModel
02 Jul 2015: Removed MiriMrsD2CModel and added MiriMrsDistortionModel.
             Duplicate dictionary entries removed by converting all data
             type keywords to uppercase.
09 Jul 2015: Restored legacy data types which are still being used.
20 Aug 2015: MiriNonlinearityModel changed to MiriLinearityModel.
06 Nov 2015: Added MiriImagingPhotometricModel and MiriPixelAreaModel.
             Changed CDP_DICT so that 'PHOTOM' with 'IM' or 'MIRIMAGE' refers
             to MiriImagingPhotometricModel.
20 Nov 2015: Added new data models and reftypes for CDP-5 delivery.
03 Dec 2015: Added MiriPowerlawColourCorrectionModel.
10 Dec 2015: Included old and new MRS distortion models.
16 Jun 2016: Old format MRS data models removed (as MIRISim no longer
             uses them).
23 Jun 2016: Added MiriPceModel.
01 Apr 2018: Separated the legacy models from the supported ones.

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""
# This module is now converted to Python 3.


# Import CDP data products
from miri.datamodels.miri_badpixel_model import \
    MiriBadPixelMaskModel
from miri.datamodels.miri_dark_reference_model import \
    MiriDarkReferenceModel
from miri.datamodels.miri_distortion_models import \
    MiriImagingDistortionModel, MiriLrsD2WModel, MiriMrsDistortionModel12, \
    MiriMrsDistortionModel34
from miri.datamodels.miri_droop_model import MiriDroopModel
from miri.datamodels.miri_flatfield_model import \
    MiriFlatfieldModel
from miri.datamodels.miri_fringe_frequencies_model import \
    MiriMrsFringeFrequenciesModel
from miri.datamodels.miri_pce_model import MiriPceModel
from miri.datamodels.miri_fluxconversion_models import \
    MiriFluxconversionModel, MiriImagingFluxconversionModel, \
    MiriImagingColourCorrectionModel, MiriPowerlawColourCorrectionModel, \
    MiriLrsFluxconversionModel, MiriMrsFluxconversionModel
from miri.datamodels.miri_photometric_models import \
    MiriPhotometricModel, MiriImagingPhotometricModel, MiriPixelAreaModel
from miri.datamodels.miri_transmission_correction_model import \
    MiriMrsTransmissionCorrectionModel
from miri.datamodels.miri_wavelength_correction_model import \
    MiriMrsWavelengthCorrectionModel
from miri.datamodels.miri_spectral_spatial_resolution_model import \
    MiriMrsResolutionModel
from miri.datamodels.miri_aperture_correction_model import \
    MiriMrsApertureCorrectionModel
from miri.datamodels.miri_reset_switch_charge_decay_model import \
    MiriResetSwitchChargeDecayModel
from miri.datamodels.miri_gain_model import MiriGainModel
from miri.datamodels.miri_ipc_model import MiriIPCModel
from miri.datamodels.miri_jump_model import MiriJumpModel
from miri.datamodels.miri_lastframe_model import MiriLastFrameModel
from miri.datamodels.miri_latent_model import \
    MiriLatentDecayModel
from miri.datamodels.miri_linearity_model import \
    MiriLinearityModel
from miri.datamodels.miri_pixel_saturation_model import \
    MiriPixelSaturationModel
from miri.datamodels.miri_psf_models import \
    MiriPointSpreadFunctionModel, MiriImagingPointSpreadFunctionModel, \
    MiriLrsPointSpreadFunctionModel, MiriMrsPointSpreadFunctionModel
from miri.datamodels.miri_readnoise_model import MiriReadnoiseModel
from miri.datamodels.miri_reset_model import MiriResetModel
from miri.datamodels.miri_straylight_model import \
    MiriMrsStraylightModel
from miri.datamodels.miri_telescope_emission_model import \
    MiriTelescopeEmissionModel

# Define a dictionary which distinguishes the different kinds of
# models, which are distinguished by detector name and filter.

# Define a dictionary giving the data type code for each of the above
# CDP data models. Some of the models need to be distinguished by
# detector and filter in addition to the type code.
CDP_DICT = { \
            # -----------------------------------------------------------
            # Currently Supported Calibration Data Products
            # -----------------------------------------------------------
            'MASK'    : MiriBadPixelMaskModel, \
            # BAD is an alias of MASK
            'BAD'     : MiriBadPixelMaskModel, \
            'DARK'    : MiriDarkReferenceModel, \
            'DISTORTION' : {'IM'  : {'P750L' : MiriLrsD2WModel, \
                                     'ANY'   : MiriImagingDistortionModel}, \
                         'SW'  : MiriMrsDistortionModel12, \
                         'LW'  : MiriMrsDistortionModel34, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsD2WModel, \
                                        'ANY'   : MiriImagingDistortionModel}, \
                         'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                         'MIRIFULONG'  : MiriMrsDistortionModel34,
                         'ANY' : MiriImagingDistortionModel }, \
            'DROOP'   : MiriDroopModel, \
            'FLAT'    : MiriFlatfieldModel, \
            # FRINGE, PIXELFLAT and SKYFLAT are all kinds of FLAT
            'FRINGE' : MiriFlatfieldModel,  \
            'PIXELFLAT' : MiriFlatfieldModel,  \
            'SKYFLAT' : MiriFlatfieldModel,  \
            'FRINGEFREQ' : MiriMrsFringeFrequenciesModel, \
            'RESET'   : MiriResetModel, \
            'RSCD'    : MiriResetSwitchChargeDecayModel, \
            'GAIN'    : MiriGainModel, \
            'IPC'     : MiriIPCModel, \
            'READNOISE' : MiriReadnoiseModel, \
            'TRACORR' : MiriMrsTransmissionCorrectionModel, \
            'WAVCORR' : MiriMrsWavelengthCorrectionModel, \
            'RESOL' : MiriMrsResolutionModel, \
            'APERCORR' : MiriMrsApertureCorrectionModel, \
            'PCE'     : MiriPceModel, \
            'PHOTOM'  : {'IM'  : {'P750L' : MiriLrsFluxconversionModel, \
                                  'ANY'   : MiriImagingPhotometricModel},
                         'SW'  : MiriMrsFluxconversionModel, \
                         'LW'  : MiriMrsFluxconversionModel, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsFluxconversionModel, \
                                        'ANY'   : MiriImagingPhotometricModel},
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel}, \
            'AREA' : MiriPixelAreaModel, \
            'COLCORR' : MiriImagingColourCorrectionModel, \
            'COLCORRPL' : MiriPowerlawColourCorrectionModel, \
            'JUMP'    : MiriJumpModel, \
            'LASTFRAME' : MiriLastFrameModel, \
            'LATENT'  : MiriLatentDecayModel, \
            'LINEARITY' : MiriLinearityModel, \
            # LIN is an alias for LINEARITY
            'LIN'     : MiriLinearityModel,  \
            'SATURATION' : MiriPixelSaturationModel, \
            # SAT is an alias for SATURATION
            'SAT'     : MiriPixelSaturationModel, \
            'STRAYMASK' : {'SW'  : MiriMrsStraylightModel, \
                         'LW'  : MiriMrsStraylightModel,
                         'MIRIFUSHORT' : MiriMrsStraylightModel, \
                         'MIRIFULONG'  : MiriMrsStraylightModel }, \
            # STRAY is an alias for STRAYMASK
            'STRAY'   : {'SW'  : MiriMrsStraylightModel, \
                         'LW'  : MiriMrsStraylightModel,
                         'MIRIFUSHORT' : MiriMrsStraylightModel, \
                         'MIRIFULONG'  : MiriMrsStraylightModel }, \
            'PSF'     : {'IM'  : {'P750L' : MiriLrsPointSpreadFunctionModel, \
                                  'ANY'   : MiriImagingPointSpreadFunctionModel}, \
                         'SW'  : MiriMrsPointSpreadFunctionModel, \
                         'LW'  : MiriMrsPointSpreadFunctionModel, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsPointSpreadFunctionModel, \
                                        'ANY'   : MiriImagingPointSpreadFunctionModel}, \
                         'MIRIFUSHORT' : MiriMrsPointSpreadFunctionModel, \
                         'MIRIFULONG'  : MiriMrsPointSpreadFunctionModel,
                         'ANY' : MiriPointSpreadFunctionModel }, \
            'PSF-OOF' : MiriImagingPointSpreadFunctionModel, \
            'PSF-MONOCHROM' : MiriLrsPointSpreadFunctionModel, \
            # -----------------------------------------------------------
            # Legacy Keywords and Data Products - Backwards Compatibility Only
            # The following section could be removed without affecting the
            # JWST pipeline.
            # -----------------------------------------------------------
            'DISTORT' : MiriImagingDistortionModel,  \
            'D2W'     : MiriLrsD2WModel,  \
            'D2C'     : {'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                         'MIRIFULONG'  : MiriMrsDistortionModel34},  \
            'WCS'     : {'IM'  : {'P750L' : MiriLrsD2WModel, \
                                  'ANY'   : MiriImagingDistortionModel}, \
                         'SW'  : MiriMrsDistortionModel12, \
                         'LW'  : MiriMrsDistortionModel12, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsD2WModel, \
                                        'ANY'   : MiriImagingDistortionModel}, \
                         'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                         'MIRIFULONG'  : MiriMrsDistortionModel34,
                         'ANY' : MiriImagingDistortionModel }, \
            'PIXFLAT' : MiriFlatfieldModel,  \
            'FRINGEFLAT' : MiriFlatfieldModel,  \
            'FLUX'    : {'IM'  : {'P750L' : MiriLrsFluxconversionModel, \
                                  'ANY'   : MiriImagingFluxconversionModel}, \
                         'SW'  : MiriMrsFluxconversionModel, \
                         'LW'  : MiriMrsFluxconversionModel, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsFluxconversionModel, \
                                        'ANY'   : MiriImagingFluxconversionModel}, \
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel,
                         'ANY' : MiriFluxconversionModel }, \
            'ABSFLUX' : {'IM'  : {'P750L' : MiriLrsFluxconversionModel, \
                                  'ANY'   : MiriImagingFluxconversionModel},
                         'SW' : MiriMrsFluxconversionModel, \
                         'LW' : MiriMrsFluxconversionModel, \
                         'MIRIMAGE'  : {'P750L' : MiriLrsFluxconversionModel, \
                                        'ANY'   : MiriImagingFluxconversionModel},
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel}, \
            'SRF'     : {'IM' : MiriLrsFluxconversionModel, \
                         'SW' : MiriMrsFluxconversionModel, \
                         'LW' : MiriMrsFluxconversionModel, \
                         'MIRIMAGE' : MiriLrsFluxconversionModel, \
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel}, \
            'IMPSF'  : MiriImagingPointSpreadFunctionModel, \
            'LRSPSF' : MiriLrsPointSpreadFunctionModel, \
            'MRSPSF' : MiriMrsPointSpreadFunctionModel, \
            'TelEm'  : MiriTelescopeEmissionModel, \
            'TEL_EMISSION'  : MiriTelescopeEmissionModel, \
            # -----------------------------------------------------------
            }
