Jump Detection Coefficients (:mod:`datamodels.miri_jump_model`)
===============================================================

.. module:: miri.datamodels.miri_jump_model

Description
~~~~~~~~~~~
This module contains the MiriJumpModel class, which describes
MIRI jump detection coefficients.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriJumpModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a jump detection CDP contains the
following HDUs

* Primary - metadata
* FINE - fine jump thresholds table
