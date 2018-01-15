Flux Conversion Data (:mod:`datamodels.miri_fluxconversion_models`)
===================================================================

.. module:: miri.datamodels.miri_fluxconversion_models

Description
~~~~~~~~~~~

This module contains the MiriFluxconversionModel, 
MiriImagingFluxconversionModel, MiriImagingColourCorrectionModel, 
MiriPowerlawColourCorrectionModel, MiriLrsFluxconversionModel and 
MiriLrsFluxconversionModel classes, which describe MIRI flux conversion 
factors. The imaging and spectroscopy models differ in the way the 
wavelength dependency is parameterized.

NOTE: See also the MIRI PCE and photometric data models.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriFluxconversionModel
   :members:

.. autoclass:: MiriImagingFluxconversionModel
   :members:

.. autoclass:: MiriImagingColourCorrectionModel
   :members:

.. autoclass:: MiriPowerlawColourCorrectionModel
   :members:

.. autoclass:: MiriLrsFluxconversionModel
   :members:

.. autoclass:: MiriMrsFluxconversionModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
Each of the flux conversion CDPs has a different file structure (see the
YAML schema files for details).

The imager flux conversion model contains a primary HDU followed by a
flux table (in (Jy/arcsec2)/(DN/s/pixel) units).

The LRS flux conversion model contains primary HDU followed by a
flux table (in DN/s/Jy/spaxel units).

The MRS flux conversion mode contains a primary HDU followed by a 
flux conversion data array, flux conversion uncertainty array and
pixel size array.
