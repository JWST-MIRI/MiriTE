Telescope Emission Data (:mod:`datamodels.miri_telescope_emission_model`)
=========================================================================

.. module:: miri.datamodels.miri_telescope_emission_model

Description
~~~~~~~~~~~
This module contains the MiriTelescopeEmissionModel class, which describes
MIRI telescope emission data.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriTelescopeEmissionModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a telescope emission CDP contains the following HDUs

* Primary - metadata
* SCI - telescope emission data (rows x columns)
* ERR - telescope emission uncertainty
* DQ - telescope emission quality
* DQ_DEF - description of the contents of the DQ array
