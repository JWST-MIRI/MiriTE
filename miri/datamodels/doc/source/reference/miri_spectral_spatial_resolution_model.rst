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
* PSF_FWHM_ALPHA - Alpha spectral resolution table
* PSF_FWHM_BETA - Beta spectral resolution table
* RESOL_DATA - Resolving power table
* MSLF_DATA - Parameters describing the MLSF profile
* PHASE1_DATA - Samples for constructing the Phase1 spline
* PHASE2_DATA - Coefficients for constructing the Phase2 polynomials
* PHASE3_DATA - Coefficients for constructing the Phase3 polynomials
* ETALON_DATA - Description of the etalon line fit
