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
03 Sep 2018: Obsolete detector names removed. There should no longer be any
             need for backwards compatibility to these old names.
             Dictionary extended to include CDP-6 variant of the MRS
             distortion models.
26 Sep 2018: Added a REFTYPE of 'FLAT-TA' for an on-board TA flat field.
18 Oct 2018: Dictionary extended to include CDP-5 variant of the MIRI
             photometric models.
30 Oct 2018: There are now different classes for different types of flat-field.
09 Nov 2018: Added 'FRINGE_CALSOURCE' to the dictionary. Moved 'PSF-MONOCHROM',
             'RESET' and 'IPC' into the OBSOLETE section. Reordered the CDPs
             in the dictionary to reflect the order on the wiki.
15 Nov 2018: Dictionary modified to include CDP-3 variant of the MRS
             straylight model.
12 Mar 2019: Added 'SPECWCS' as an alias for 'DISTORTION'.

@author: Steven Beard (UKATC), Vincent Geers (DIAS)

"""

# Import CDP data products
from miri.datamodels.miri_badpixel_model import \
    MiriBadPixelMaskModel
from miri.datamodels.miri_dark_reference_model import \
    MiriDarkReferenceModel
from miri.datamodels.miri_distortion_models import \
    MiriImagingDistortionModel, MiriLrsD2WModel, MiriMrsDistortionModel12, \
    MiriMrsDistortionModel34, MiriMrsDistortionModel12_CDP6, \
    MiriMrsDistortionModel34_CDP6
from miri.datamodels.miri_droop_model import MiriDroopModel
from miri.datamodels.miri_flatfield_model import \
    MiriFlatfieldModel, MiriSkyFlatfieldModel, MiriFringeFlatfieldModel, \
    MiriTargetFlatfieldModel
from miri.datamodels.miri_fringe_frequencies_model import \
    MiriMrsFringeFrequenciesModel
from miri.datamodels.miri_pce_model import MiriPceModel
from miri.datamodels.miri_fluxconversion_models import \
    MiriFluxconversionModel, MiriImagingFluxconversionModel, \
    MiriImagingColourCorrectionModel, MiriPowerlawColourCorrectionModel, \
    MiriLrsFluxconversionModel, MiriMrsFluxconversionModel
from miri.datamodels.miri_photometric_models import \
    MiriPhotometricModel, MiriPhotometricModel_CDP5, MiriImagingPhotometricModel, \
    MiriLrsPhotometricModel, MiriPixelAreaModel
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
    MiriMrsStraylightModel, MiriMrsStraylightModel_CDP3
from miri.datamodels.miri_telescope_emission_model import \
    MiriTelescopeEmissionModel

# Define a dictionary giving the data type code for each of the above
# CDP data models. Some of the models need to be distinguished by
# detector and filter in addition to the data type.
CDP_DICT = { \
            # -----------------------------------------------------------
            # Currently Supported Calibration Data Products
            # -----------------------------------------------------------
            'MASK'    : MiriBadPixelMaskModel, \
            # BAD is an alias of MASK
            'BAD'     : MiriBadPixelMaskModel, \
            'DARK'    : MiriDarkReferenceModel, \
            'LINEARITY' : MiriLinearityModel, \
            # LIN is an alias for LINEARITY
            'LIN'     : MiriLinearityModel,  \
            'SATURATION' : MiriPixelSaturationModel, \
            # SAT is an alias for SATURATION
            'SAT'     : MiriPixelSaturationModel, \
            'DROOP'   : MiriDroopModel, \
            'LATENT'  : MiriLatentDecayModel, \
            'JUMP'    : MiriJumpModel, \
            'GAIN'    : MiriGainModel, \
            'READNOISE' : MiriReadnoiseModel, \
            'LASTFRAME' : MiriLastFrameModel, \
            # TODO: Remove this cdprelease complexity after CDP-7 release
#             'RSCD'    : {'MIRIMAGE'    : {'ANY' : {'6'   : MiriResetSwitchChargeDecayModel_CDP6, \
#                                                    '7'   : MiriResetSwitchChargeDecayModel, \
#                                                    'ANY' : MiriResetSwitchChargeDecayModel}}, \
#                          'MIRIFULONG'  : {'ANY' : {'6'   : MiriResetSwitchChargeDecayModel_CDP6, \
#                          'MIRIFUSHORT' : {'ANY' : {'6'   : MiriResetSwitchChargeDecayModel_CDP6, \
#                                                    '7'   : MiriResetSwitchChargeDecayModel, \
#                                                    'ANY' : MiriResetSwitchChargeDecayModel}}, \
#                          'MIRIFULONG'  : {'ANY' : {'6'   : MiriResetSwitchChargeDecayModel_CDP6, \
#                                                    '7'   : MiriResetSwitchChargeDecayModel, \
#                                                    'ANY' : MiriResetSwitchChargeDecayModel}}, \
#                          'ANY'        : MiriResetSwitchChargeDecayModel }, \
            'RSCD'    : MiriResetSwitchChargeDecayModel, \
            'PCE'     : MiriPceModel, \
            'PSF'     : {'MIRIMAGE'    : {'P750L' : MiriLrsPointSpreadFunctionModel, \
                                          'ANY'   : MiriImagingPointSpreadFunctionModel}, \
                         'MIRIFUSHORT' : MiriMrsPointSpreadFunctionModel, \
                         'MIRIFULONG'  : MiriMrsPointSpreadFunctionModel,
                         'ANY' : MiriPointSpreadFunctionModel }, \
            'PSF-OOF' : MiriImagingPointSpreadFunctionModel, \
            # TODO: Remove this cdprelease complexity after CDP-7 release
            'DISTORTION' : {'MIRIMAGE'    : {'P750L' : MiriLrsD2WModel, \
                                             'ANY'   : MiriImagingDistortionModel}, \
                            'MIRIFUSHORT' : {'ANY'   : {'6'   : MiriMrsDistortionModel12_CDP6, \
                                                        '7'   : MiriMrsDistortionModel12, \
                                                        '8'   : MiriMrsDistortionModel12_CDP8, \
                                                        '8B'  : MiriMrsDistortionModel12_CDP8, \
                                                        'ANY' : MiriMrsDistortionModel12}}, \
                            'MIRIFULONG'  : {'ANY'   : {'6'   : MiriMrsDistortionModel34_CDP6, \
                                                        '7'   : MiriMrsDistortionModel34, \
                                                        '8'   : MiriMrsDistortionModel34_CDP8, \
                                                        '8B'  : MiriMrsDistortionModel34_CDP8, \
                                                        'ANY' : MiriMrsDistortionModel34}}, \
                            'ANY'         : MiriImagingDistortionModel }, \
# Previous code without the extra level for CDP version
#             'DISTORTION' : {'MIRIMAGE'  : {'P750L' : MiriLrsD2WModel, \
#                                            'ANY'   : MiriImagingDistortionModel}, \
#                             'MIRIFUSHORT' : MiriMrsDistortionModel12, \
#                             'MIRIFULONG'  : MiriMrsDistortionModel34, \
#                             'ANY'         : MiriImagingDistortionModel }, \
            # The REFTYPE for spectroscopy distortion has changed to SPECWCS
            # so make this an alias for DISTORTION.
            'SPECWCS' :    {'MIRIMAGE'  : {'P750L' : MiriLrsD2WModel, \
                                           'ANY'   : MiriImagingDistortionModel}, \
                            'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                            'MIRIFULONG'  : MiriMrsDistortionModel34, \
                            'ANY'         : MiriImagingDistortionModel }, \
            'AREA' : MiriPixelAreaModel, \
            # TODO: Remove this cdprelease complexity after CDP-7 release
            'PHOTOM'  : {'MIRIMAGE'    : {'P750L' : {'5'   : MiriLrsFluxconversionModel, \
                                                     '6'   : MiriLrsFluxconversionModel, \
                                                     '7'   : MiriLrsPhotometricModel, \
                                                     'ANY' : MiriLrsPhotometricModel}, \
                                          'ANY'   : {'5'   : MiriPhotometricModel_CDP5, \
                                                     '6'   : MiriPhotometricModel_CDP5, \
                                                     '7'   : MiriImagingPhotometricModel, \
                                                     'ANY' : MiriImagingPhotometricModel}}, \
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel }, \
# Previous code without the extra level for CDP version
#             'PHOTOM'  : {'MIRIMAGE'    : {'P750L' : MiriLrsFluxconversionModel, \
#                                           'ANY'   : MiriImagingPhotometricModel},
#                          'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
#                          'MIRIFULONG'  : MiriMrsFluxconversionModel }, \
            'COLCORR' : MiriImagingColourCorrectionModel, \
            'COLCORRPL' : MiriPowerlawColourCorrectionModel, \
            'FLAT'    : MiriFlatfieldModel, \
            # FRINGE, PIXELFLAT, SKYFLAT and FLAT-TA are all kinds of FLAT
            'FRINGE' : MiriFringeFlatfieldModel,  \
            'FRINGE_CALSOURCE' : MiriFringeFlatfieldModel,  \
            'PIXELFLAT' : MiriFlatfieldModel,  \
            'SKYFLAT' : MiriSkyFlatfieldModel,  \
            'FLAT-TA' : MiriTargetFlatfieldModel,  \
            'FRINGEFREQ' : MiriMrsFringeFrequenciesModel, \
            # STRAYMASK corresponds to the new format CDP and STRAY to the old one.
            'STRAYMASK' : {'MIRIFUSHORT' : MiriMrsStraylightModel, \
                           'MIRIFULONG'  : MiriMrsStraylightModel }, \
            'STRAY'   : {'MIRIFUSHORT' : MiriMrsStraylightModel_CDP3, \
                         'MIRIFULONG'  : MiriMrsStraylightModel_CDP3 }, \
            'TRACORR' : MiriMrsTransmissionCorrectionModel, \
            'WAVCORR' : MiriMrsWavelengthCorrectionModel, \
            'RESOL'   : MiriMrsResolutionModel, \
            'APERCORR' : MiriMrsApertureCorrectionModel, \
            # -----------------------------------------------------------
            # Legacy Keywords and Data Products - Backwards Compatibility Only
            # The following section could be removed without affecting the
            # JWST pipeline.
            # -----------------------------------------------------------
            'RESET'   : MiriResetModel, \
            'IPC'     : MiriIPCModel, \
            'DISTORT' : MiriImagingDistortionModel,  \
            'D2W'     : MiriLrsD2WModel,  \
            'D2C'     : {'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                         'MIRIFULONG'  : MiriMrsDistortionModel34},  \
            'WCS'     : {'MIRIMAGE'    : {'P750L' : MiriLrsD2WModel, \
                                          'ANY'   : MiriImagingDistortionModel}, \
                         'MIRIFUSHORT' : MiriMrsDistortionModel12, \
                         'MIRIFULONG'  : MiriMrsDistortionModel34,
                         'ANY' : MiriImagingDistortionModel }, \
            'PIXFLAT' : MiriFlatfieldModel,  \
            'FRINGEFLAT' : MiriFringeFlatfieldModel,  \
            'FLUX'    : {'MIRIMAGE'    : {'P750L' : MiriLrsFluxconversionModel, \
                                          'ANY'   : MiriImagingFluxconversionModel}, \
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel,
                         'ANY' : MiriFluxconversionModel }, \
            'ABSFLUX' : {'MIRIMAGE'    : {'P750L' : MiriLrsFluxconversionModel, \
                                          'ANY'   : MiriImagingFluxconversionModel},
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel }, \
            'SRF'     : {'MIRIMAGE'    : MiriLrsFluxconversionModel, \
                         'MIRIFUSHORT' : MiriMrsFluxconversionModel, \
                         'MIRIFULONG'  : MiriMrsFluxconversionModel }, \
            'IMPSF'   : MiriImagingPointSpreadFunctionModel, \
            'LRSPSF'  : MiriLrsPointSpreadFunctionModel, \
            'PSF-MONOCHROM' : MiriLrsPointSpreadFunctionModel, \
            'MRSPSF'  : MiriMrsPointSpreadFunctionModel, \
            'TelEm'   : MiriTelescopeEmissionModel, \
            'TEL_EMISSION' : MiriTelescopeEmissionModel, \
            # -----------------------------------------------------------
            }
