Flat-Field Data (:mod:`datamodels.miri_flatfield_model`)
========================================================

.. module:: miri.datamodels.miri_flatfield_model

Description
~~~~~~~~~~~
This module contains the MiriFlatfieldModel class, which describes
MIRI flat-field data. The same model is used for pixel flats, fringe
flats and sky flats (differing only in the data type stored in the
metadata).

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriFlatfieldModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
flat_reference_flags - Defines the MIRI/JWST flat-field reference flags

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a flat-field CDP contains the following HDUs

* Primary - metadata
* SCI - flat-field (rows x columns)
* ERR - flat-field uncertainty
* DQ - flat-field quality
* DQ_DEF - description of the contents of the DQ array
