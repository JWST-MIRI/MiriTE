Aperture Correction Data (:mod:`datamodels.miri_aperture_correction_model`)
===========================================================================

.. module:: miri.datamodels.miri_aperture_correction_model

Description
~~~~~~~~~~~
This module contains the MiriMrsApertureCorrectionModel class, which 
describes a MIRI aperture correction model.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsApertureCorrectionModel
   :members:

Functions
~~~~~~~~~
When stored as a FITS file, an aperture correction CDP contains the following HDUs

* Primary - metadata
* APER_CORR - aperture correction table
