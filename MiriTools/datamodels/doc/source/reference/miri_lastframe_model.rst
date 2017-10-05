Last Frame Correction Model (:mod:`datamodels.miri_lastframe_model`)
====================================================================

.. module:: miri.datamodels.miri_lastframe_model

Description
~~~~~~~~~~~
This module contains the MiriLastFrameModel class, which describes
MIRI last frame correction coefficients.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriLastFrameModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
lastframe_reference_flags - Reference flags used by the MiriLastFrameModel.

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a last frame CDP contains the following HDUs

* Primary - metadata
* SCI - last frame correction data (rows x columns)
* ERR - last frame uncertainty
* DQ - last frame quality
* DQ_DEF - description of the contents of the DQ array
