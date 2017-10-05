Photometric Conversion Data (:mod:`datamodels.miri_pce_model`)
==============================================================

.. module:: miri.datamodels.miri_pce_model

Description
~~~~~~~~~~~

This module contains the MiriPceModel class, which describes the 
reflectance of the optics, the transmission of the spectral filters and 
the QE of the detectors as a function of wavelength.

NOTE: The PCE data model is expected to be used by the JWST exposure time
calculator (ETC). Compare with the MIRI photometric models data models.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriPceModel
   :members:

Functions
~~~~~~~~~
.. autofunction:: ascii_to_pce

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a MIRI PCE model contains the following HDUs

* Primary - metadata (including ETCNAME=ETC file name,
  COMPNAME=ETC component name, LITREF=ETC literary reference and
  SYSTEM=Is data intended for use in ETC?)
* DATA_TABLE - Photon conversion efficiency table
