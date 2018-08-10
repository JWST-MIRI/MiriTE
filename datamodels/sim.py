#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Import the simulator support models. Adding this module to datamodels
allows external software to import all the MIRI simulator models with one
import line, for example

   from datamodels.sim import MiriIlluminationModel, MiriExposureModel
   
instead of needing to know the names of all the individual source files.

:History:

23 May 2013: Created
29 Jul 2013: SIM_DICT added.
06 Jul 2015: Duplicate dictionary entries removed by converting all data
             type keywords to uppercase.
20 May 2016: Added ramp and slope data to SIM_DICT.

@author: Steven Beard (UKATC)

"""

# Import simulator support data products
from miri.datamodels.miri_illumination_model import MiriIlluminationModel
from miri.datamodels.miri_exposure_model import MiriExposureModel
from miri.datamodels.miri_measured_model import MiriRampModel, MiriSlopeModel

# Define a dictionary giving the data type code for each of the above
# simulator data models.
SIM_DICT = {'ILLUMINATION' : MiriIlluminationModel, \
            'RAW' : MiriExposureModel, \
            'RAMP' : MiriRampModel, \
            'RAMP (LEVEL 1A)' : MiriRampModel, \
            'RAMP (LEVEL 1B)' : MiriRampModel, \
            'SLOPE' : MiriSlopeModel, \
            'SLOPE (LEVEL 2))' : MiriSlopeModel \
            }
