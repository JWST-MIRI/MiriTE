Detector Read noise Data (:mod:`datamodels.miri_readnoise_model`)
=================================================================

.. module:: miri.datamodels.miri_readnoise_model

Description
~~~~~~~~~~~
This module contains the MiriReadnoiseModel class, which describes
MIRI detector amplifier read noise.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriReadnoiseModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a READNOISE CDP contains the following HDUs

* Primary - metadata
* SCI - read noise (rows x columns)
* ERR - read noise uncertainty
* DQ - read noise quality
* DQ_DEF - description of the contents of the DQ array
