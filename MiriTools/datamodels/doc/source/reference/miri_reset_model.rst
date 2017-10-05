Detector Reset Model (:mod:`datamodels.miri_reset_model`)
=========================================================

.. module:: miri.datamodels.miri_reset_model

Description
~~~~~~~~~~~
This module contains the MiriResetModel class, which describes
MIRI detector reset level.

Compare with the MIRI Reset Switch Charge Decay (RSCD) model.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriResetModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
reset_reference_flags - Reference flags for the reset data model.

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a RESET CDP contains the following HDUs

* Primary - metadata
* SCI - reset data (rows x columns)
* ERR - reset data uncertainty
* DQ - reset data quality
* DQ_DEF - description of the contents of the DQ array
