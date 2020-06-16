Bad Pixel Data (:mod:`datamodels.miri_badpixel_model`)
======================================================

.. module:: miri.datamodels.miri_badpixel_model

Description
~~~~~~~~~~~
This module contains the MiriBadPixelMaskModel class, which describes a 
MIRI bad pixel mask.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriBadPixelMaskModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
mask_reference_flags - Defines the MIRI/JWST pixel mask reference flags

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a bad pixel mask contains the following HDUs

* Primary - metadata
* DQ - bad pixel mask (rows x columns)
* DQ_DEF - description of the contents of the bad pixel mask

