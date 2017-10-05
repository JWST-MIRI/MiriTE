Spectral and Spatial Resolution Data (:mod:`datamodels.miri_spectral_spatial_resolution_model`)
===============================================================================================

.. module:: miri.datamodels.miri_spectral_spatial_resolution_model

Description
~~~~~~~~~~~
This module contains the MiriMrsResolutionModel class, which describes a 
MIRI MRS spectral and spatial resolution model.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsResolutionModel
   :members:

Functions
~~~~~~~~~
When stored as a FITS file, a spectral spatial resolution CDP contains the
following HDUs

* Primary - metadata
* RESOLVING_POWER - resolving power table
* PSF_FWHM_ALPHA - alpha spectral resolution table
* PSF_FWHM_BETA - beta spectral resolution table
