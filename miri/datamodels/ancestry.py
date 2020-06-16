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


@author: Steven Beard (UKATC)

"""

# Define a dictionary translating each MIRI data model name into its
# equivalent JWST data model. Data models without a JWST equivalent are
# not listed here.
#
# TODO: The gaps in this mapping need to be filled in.
#       Uncertain mappings are commented out.
#
PARENT_DICT = {'MiriDataModel' : 'DataModel', \
               'MiriRampModel' : 'RampModel', \
               'MiriExposureModel' : 'RampModel', \
               'MiriSlopeModel' : 'ImageModel', \
               'MiriImagingApertureCorrectionModel' : 'MirImgApcorrModel', \
               'MiriMrsApertureCorrectionModel' : '', \
               'MiriBadPixelMaskModel' :    'MaskModel', \
               'MiriDarkReferenceModel' : 'DarkMIRIModel', \
               'MiriImagingDistortionModel' : 'DistortionModel', \
               'MiriMrsDistortionModel12' : 'DistortionMRSModel', \
               'MiriMrsDistortionModel34' : 'DistortionMRSModel', \
#                'MiriDroopModel' : '', \
               'MiriFlatfieldModel' : 'FlatModel', \
               'MiriFringeFlatfieldModel' : 'FlatModel', \
               'MiriSkyFlatfieldModel' : 'FlatModel', \
               'MiriTargetFlatfieldModel' : 'FlatModel', \
               'MiriMrsFluxconversionModel' : 'MiriMrsPhotomModel', \
#                'MiriMrsFringeFrequenciesModel' : 'FringeModel', \
               'MiriGainModel' : 'GainModel', \
               'MiriIPCModel' : 'IPCModel', \
#                'MiriJumpModel' : '', \
               'MiriLastFrameModel' : 'LastFrameModel', \
#                'MiriLatentDecayModel' : '', \
               'MiriLinearityModel' : 'LinearityModel', \
#                'MiriPceModel' : 'PhotomModel', \
               'MiriImagingPhotometricModel' : 'MiriImgPhotomModel', \
               'MiriPixelAreaModel' : 'PixelAreaModel',
               'MiriPixelSaturationModel' : 'SaturationModel', \
               'MiriReadnoiseModel' : 'ReadnoiseModel', \
#                'MiriPointSpreadFunctionModel' : 'PsfMaskModel', \
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
