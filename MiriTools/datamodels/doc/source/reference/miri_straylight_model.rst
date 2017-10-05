Stray Light Model (:mod:`datamodels.miri_straylight_model`)
===========================================================

.. module:: miri.datamodels.miri_straylight_model

Description
~~~~~~~~~~~
This module contains the MiriMrsStraylightModel class, which describes
MIRI MRS stray light properties.

Further stray light models may be added here.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsStraylightModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
straylight_reference_flags - Reference flags for the straylight model.

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a stray light CDP contains the following HDUs

* Primary - metadata
* DQ - stray light mask (rows x columns)
* DQ_DEF - description of the contents of the stray light mask
