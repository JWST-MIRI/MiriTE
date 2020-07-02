Latent Decay Model (:mod:`datamodels.miri_latent_model`)
========================================================

.. module:: miri.datamodels.miri_latent_model

Description
~~~~~~~~~~~
This module contains the MiriLatentDecay and MiriLatentTable classes,
which describe the MIRI detector latent decay model.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriLatentTable
   :members:

.. autoclass:: MiriLatentDecayModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a LATENT CDP contains the following HDUs

* Primary - metadata
* LATENT1 - latent correction table 1
* LATENT2 - latent correction table 2
* LATENT3 - latent correction table 3
* LATENT4 - latent correction table 4
