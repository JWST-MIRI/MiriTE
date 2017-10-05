MIRI Simulation Exposure Data Module (:mod:`datamodels.miri_exposure_model`)
============================================================================

.. module:: miri.datamodels.miri_exposure_model

Description
~~~~~~~~~~~
This module contains the MiriExposureModel class, which is designed to be
used by simulators to store simulated JWST level 1 exposure data.
It is based on the MiriRampDataModel.

Objects
~~~~~~~
.. autoclass:: MiriExposureModel
   :members:

Functions
~~~~~~~~~
TBD.

Data formats
~~~~~~~~~~~~
MIRI exposure data has the same format as JWST level 1 exposure data.
A FITS file of exposure data contains the following HDUs

* Primary - metadata
* SCI - ramp data (integrations x groups x rows x columns)
* ERR - ramp data uncertainty
* PIXEL_DQ - ramp data pixel quality (rows x columns)
* PIXELDQ_DEF - description of the contents of the PIXEL_DQ array
* GROUP_DQ - ramp data group quality (integrations x groups x rows x columns)
* GROUPDQ_DEF - description of the contents of the GROUP_DQ array
* REFOUT - reference output data

PIXEL_DQ and GROUP_DQ are optional (and are normally missing from
raw JWST ramp data).
