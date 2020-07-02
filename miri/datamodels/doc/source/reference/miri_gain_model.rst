Detector Gain Data (:mod:`datamodels.miri_gain_model`)
======================================================

.. module:: miri.datamodels.miri_gain_model

Description
~~~~~~~~~~~
This module contains the MiriGainModel class, which describes
MIRI detector amplifier gains.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriGainModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a GAIN CDP contains the following HDUs

* Primary - metadata
* SCI - gain (rows x columns)
* ERR - gain uncertainty
* DQ - gain quality
* DQ_DEF - description of the contents of the DQ array
