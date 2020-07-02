Point Spread Function Data (:mod:`datamodels.miri_psf_models`)
==============================================================

.. module:: miri.datamodels.miri_psf_models

Description
~~~~~~~~~~~
This module contains the MiriPointSpreadFunctionModel and related classes,
which describe the JWST/MIRI point spread function measurements.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriPointSpreadFunctionModel
   :members:

.. autoclass:: MiriImagingPointSpreadFunctionModel
   :members:

.. autoclass:: MiriLrsPointSpreadFunctionModel
   :members:

.. autoclass:: MiriMrsPointSpreadFunctionModel
   :members:

Functions
~~~~~~~~~
None

Global Data
~~~~~~~~~~~
psf_reference_flags - Reference flags for the PSF data models.

Data formats
~~~~~~~~~~~~
Each of the PSF CDPs has a different file structure (see the YAML
schema files for details).

The imager PSF model contains a primary HDU followed by a SCI
array containing PSF data (stacks x rows x columns) and a PSF_LUT
lookup table.

The LRS PSF model contains primary HDU followed by a SCI array
containing PSF data (planes x rows x columns).

The MRS PSF mode contains a primary HDU followed by a SCI array
containing PSF data (wavelength x rows x columns).
