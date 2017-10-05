Wavelength Correction (:mod:`datamodels.miri_wavelength_correction__model`)
===========================================================================

.. module:: miri.datamodels.miri_wavelength_correction_model

Description
~~~~~~~~~~~
This module contains the MiriMrsWavelengthCorrectionModel class, which describes
across slice wavelength corrections to be applied to MIRI MRS data.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsWavelengthCorrectionModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a wavelength correction CDP contains the
following HDUs

* Primary - metadata
* WAVCORR_OPTICAL - Wavelength correction optical parameters table
* WAVCORR_XSLICE - Wavelength correction cross slice table
* WAVCORR_SHIFT - Wavelength correction spectral shift table
