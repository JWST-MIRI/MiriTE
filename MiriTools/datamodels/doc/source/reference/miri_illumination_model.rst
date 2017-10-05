MIRI Simulation Illumination Data Module (:mod:`datamodels.miri_illumination_model`)
====================================================================================

.. module:: miri.datamodels.miri_illumination_model

Description
~~~~~~~~~~~
This module contains the MiriIlluminationModel and 
MiriIlluminationFringingModel classes, which are designed to be used by 
simulators to describe a detector illumination pattern.

MiriIlluminationModel is the basic data model, which describes the intensity
and wavelength of the illumination. MiriIlluminationFringingModel adds
direction information, which can be used to simulate detector fringing
effects.

Objects
~~~~~~~
.. autoclass:: MiriIlluminationModel
   :members:

.. autoclass:: MiriIlluminationFringingModel
   :members:

Functions
~~~~~~~~~
None.

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a detector illumination model contains the
following HDUs

* Primary - metadata
* INTENSITY - flux intensity (layers x rows x columns)
* WAVELENGTH - flux wavelength (layers x rows x columns)

The INTENSITY and WAVELENGTH arrays must broadcast onto each other.
