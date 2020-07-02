Linearity Data (:mod:`datamodels.miri_linearity_model`)
=======================================================

.. module:: miri.datamodels.miri_linearity_model

Description
~~~~~~~~~~~
This module contains the MiriLinearityModel class, which describes
MIRI nonlinearity coefficents.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriLinearityModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
linearity_reference_flags - Defines the MIRI/JWST linearity reference flags

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a non-linearity CDP contains the following HDUs

* Primary - metadata
* SCI - nonlinearity coefficients array (order x rows x columns)
* ERR - nonlinearity coefficients uncertainty
* DQ - nonlinearity coefficients quality
* DQ_DEF - description of the contents of the DQ array
