#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Define the relationship between the MIRI data models and the JWST
data models defined by STScI. This relationship is useful for populating
the DATAMODL metadata.

:History:

30 Jan 2019: Created
22 Mar 2019: Return a blank string when a data model does not have an
             STScI equivalent.
13 Dec 2019: MIRIRampModel replaced with RampModel.
05 Aug 2021: JWST PsfMaskModel uncommented as experiment, then recommented.
02 Sep 2021: Corrections from ticket jwst #34 (from Misty Cracraft):
             The MiriFringeFlatfieldModel should map toFringeModel,
             not FlatModel.
             MiriMrsApertureCorrectionModel should map to MirMrsApcorrModel.
             MiriImagingPhotometricModel should be MirImgPhotomModel,
             not MiriImgPhotomModel (extra i).
             MiriMrsFluxconversionModel should be MirMrsPhotomModel,
             not MiriMrsPhotomModel (extra i).
28 Sep 2021: Added missing MiriLrsPhotometricModel/MirLrsPhotomModel and
             MiriLrsD2WModel/SpecwcsModel
22 Nov 2021: MiriLrsPathlossCorrectionModel added.

@author: Steven Beard (UKATC)

"""

# Define a dictionary translating each MIRI data model name into its
# equivalent JWST data model. Data models without a JWST equivalent are
# not listed here.
#
# TODO: The gaps in this mapping need to be filled in.
#       Uncertain mappings are commented out.
#       MIRI-only data models are deliberately matched to an empty string.
#
PARENT_DICT = {'MiriDataModel' : 'DataModel', \
               'MiriRampModel' : 'RampModel', \
               'MiriExposureModel' : 'RampModel', \
               'MiriSlopeModel' : 'ImageModel', \
               'MiriImagingApertureCorrectionModel' : 'MirImgApcorrModel', \
               'MiriLrsApertureCorrectionModel': 'MirLrsApcorrModel',
               'MiriMrsApertureCorrectionModel': 'MirMrsApcorrModel',
               'MiriLrsPathlossCorrectionModel': 'MirLrsPathlossModel',
               'MiriBadPixelMaskModel' :    'MaskModel', \
               'MiriDarkReferenceModel' : 'DarkMIRIModel', \
               'MiriImagingDistortionModel' : 'DistortionModel', \
               'MiriLrsD2WModel' : 'SpecwcsModel', \
               'MiriMrsDistortionModel12' : 'DistortionMRSModel', \
               'MiriMrsDistortionModel34' : 'DistortionMRSModel', \
#                'MiriDroopModel' : '', \
               'MiriFlatfieldModel' : 'FlatModel', \
               'MiriFringeFlatfieldModel' : 'FringeModel', \
               'MiriSkyFlatfieldModel' : 'FlatModel', \
               'MiriTargetFlatfieldModel' : 'FlatModel', \
               'MiriMrsFluxconversionModel' : 'MirMrsPhotomModel', \
#                'MiriMrsFringeFrequenciesModel' : 'FringeModel', \
               'MiriGainModel' : 'GainModel', \
               'MiriIPCModel' : 'IPCModel', \
#                'MiriJumpModel' : '', \
               'MiriLastFrameModel' : 'LastFrameModel', \
#                'MiriLatentDecayModel' : '', \
               'MiriLinearityModel' : 'LinearityModel', \
#                'MiriPceModel' : 'PhotomModel', \
               'MiriImagingPhotometricModel' : 'MirImgPhotomModel', \
               'MiriLrsPhotometricModel' : 'MirLrsPhotomModel', \
               'MiriPixelAreaModel' : 'PixelAreaModel',
               'MiriPixelSaturationModel' : 'SaturationModel', \
               'MiriReadnoiseModel' : 'ReadnoiseModel', \
               'MiriPointSpreadFunctionModel' : 'PsfMaskModel', \
#                'MiriImagingPointSpreadFunctionModel' : 'PsfMaskModel', \
#                'MiriLrsPointSpreadFunctionModel' : 'PsfMaskModel', \
#                'MiriMrsPointSpreadFunctionModel' : 'PsfMaskModel', \
               'MiriResetModel' : 'ResetModel', \
               'MiriResetSwitchChargeDecayModel' : 'RSCDModel', \
               'MiriMrsResolutionModel' : 'MiriResolutionModel', \
               'MiriMrsStraylightModel' : 'StrayLightModel', \
#                'MiriTelescopeEmissionModel' : '', \
#                'MiriMrsTransmissionCorrectionModel' : '', \
               'MiriMrsWavelengthCorrectionModel' : 'WaveCorrModel', \
               }

def get_my_model_type( input_model_name ):
    """
    
    Global function which compares a data model class with the above
    dictionary and returns the equivalent STScI model name
    
    If the model name doesn't match, a blank string is returned
    
    """
    output_model_name = ' '
    if input_model_name:
        if input_model_name in list(PARENT_DICT.keys()):
            output_model_name = PARENT_DICT[input_model_name]
    #print("get_my_model_type: given \'%s\', returned \'%s\'" % (input_model_name, output_model_name ))
    return output_model_name
