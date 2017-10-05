Droop Correction Coefficients (:mod:`datamodels.miri_droop_model`)
==================================================================

.. module:: miri.datamodels.miri_droop_model

Description
~~~~~~~~~~~
This module contains the MiriDroopModel class, which describes
MIRI droop correction coefficients.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriDroopModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a DROOP CDP contains the following HDUs

* Primary - metadata
* DROOP - droop correction table
