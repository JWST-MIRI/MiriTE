Transmission Correction (:mod:`datamodels.miri_transmission_correction__model`)
================================================================================

.. module:: miri.datamodels.miri_transmission_correction_model

Description
~~~~~~~~~~~
This module contains the MiriMrsTransmissionCorrectionModel class, which describes
transmission corrections to be applied to MIRI MRS data.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsTransmissionCorrectionModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a transmission correction CDP contains the
following HDUs

* Primary - metadata
* TRA_CORR - transmission correction table
