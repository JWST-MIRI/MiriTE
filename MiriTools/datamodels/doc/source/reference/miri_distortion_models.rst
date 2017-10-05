Distortion Coefficients (:mod:`datamodels.miri_distortion_models`)
===================================================================

.. module:: miri.datamodels.miri_distortion_models

Description
~~~~~~~~~~~
This module contains the MiriImagingDistortionModel class and related classes,
which describe MIRI distortion coefficients.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriImagingDistortionModel
   :members:

.. autoclass:: MiriLrsD2WModel
   :members:

.. autoclass:: MiriMrsDistortionModel12
   :members:

.. autoclass:: MiriMrsDistortionModel34
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
Each of the distortion CDPs has a different file structure (see the YAML
schema files for details).

The imager distortion model contains a primary HDU followed by B, A,
T, M , BI, AI, TI matrix arrays, followed by a table of boresight offsets.

The LRS distortion model contains primary HDU followed by a wavelength table.

The MRS distortion mode contains a primary HDU followed by a slicenumber
array, followed by a field of view array for each channel, followed by
alpha and lambda coefficient tables for each channel, followed by x and y
tables for each channel, followed by albe_to_xanyan and xanyan_to_albe
conversion tables.

