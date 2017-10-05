Pixel Saturation Model (:mod:`datamodels.miri_pixel_saturation_model`)
======================================================================

.. module:: miri.datamodels.miri_pixel_saturation_model

Description
~~~~~~~~~~~
This module contains the MiriPixelSaturationModel class, which contains
pixel saturation data.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriPixelSaturationModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
saturation_reference_flags - Defines the MIRI/JWST saturation reference flags

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a pixel saturation CDP contains the following HDUs

* Primary - metadata
* ERR - pixel saturation error (rows x columns)
