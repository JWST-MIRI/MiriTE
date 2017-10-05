Flux Conversion Data (:mod:`datamodels.miri_photmetric_models`)
===============================================================

.. module:: miri.datamodels.miri_photometric_models

Description
~~~~~~~~~~~

This module contains the MiriPhotometricModel classes, which describe MIRI
photmetric conversion factors. It also contains the MiriPixelAreaModel
classes.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriPhotometricModel
   :members:

.. autoclass:: MiriImagingPhotometricModel
   :members:

.. autoclass:: MiriPixelAreaModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a MIRI photometric model contains the following HDUs

* Primary - metadata
* PHOTOM - Photometric flux conversion table
