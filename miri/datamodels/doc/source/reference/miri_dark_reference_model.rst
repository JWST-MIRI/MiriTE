Dark Reference Data (:mod:`datamodels.miri_dark_reference_model`)
=================================================================

.. module:: miri.datamodels.miri_dark_reference_model

Description
~~~~~~~~~~~
This module contains the MiriDarkReferenceModel class, which describes
MIRI dark reference data.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

WARNING: File containing dark reference data tend to be very large.
Use slice operations to read only the relevant parts of the data
structure. See the memory management part of the
jwst.datamodels documentation for details.

Objects
~~~~~~~
.. autoclass:: MiriDarkReferenceModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
dark_reference_flags - Defines the MIRI/JWST dark reference flags

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a DARK CDP contains the following HDUs

* Primary - metadata
* SCI - dark current (integrations x groups x rows x columns)
* ERR - dark current uncertainty
* DQ - dark current quality
* DQ_DEF - description of the contents of the DQ array
