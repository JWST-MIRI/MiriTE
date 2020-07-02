Detector Interpixel Capacitance (:mod:`datamodels.miri_ipc_model`)
==================================================================

.. module:: miri.datamodels.miri_ipc_model

Description
~~~~~~~~~~~
This module contains the MiriIPCModel class, which describes
MIRI detector Inter-Pixel Capacitance (IPC) kernel.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriIPCModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, an IPC CDP contains the following HDUs

* Primary - metadata
* SCI - IPC deconvolution kernel (rows x columns)
