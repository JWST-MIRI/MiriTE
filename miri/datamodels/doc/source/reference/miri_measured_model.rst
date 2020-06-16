MIRI Measured Data Module (:mod:`datamodels.miri_measured_model`)
=================================================================

.. module:: miri.datamodels.miri_measured_model

Description
~~~~~~~~~~~

This module contains the MiriSimpleModel, MiriMeasuredModel, 
MiriRampModel and MiriSlopeModel, classes, which describe generic MIRI 
data models containing a data array, error array and data quality array.
(The ramp model separates the pixel and group data quality into two arrays.)
The classes also include functions to support simple arithmetic operations,
mainly for interactive use.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriSimpleModel
   :members:

.. autoclass:: MiriMeasuredModel
   :members:

.. autoclass:: MiriRampModel
   :members:

.. autoclass:: MiriSlopeModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
All data models derived from MiriMeasuredModel contain the following
data arrays and tables:

   * data: The main data array, which is saved to a FITS HDU named 'SCI'.

   * err: An array containing the errors associated with the data array,
     which is saved to a FITS HDU named 'ERR'.

   * dq: A data quality array associated with the data array, which is
     saved to a FITS HDU named 'DQ'.

   * dq_def: A table describing the meaning of the flags contained in
     the dq array, which is saved to a FITS binary table HDU named
     DQ_DEF.

The err and dq arrays must be broadcastable onto the data array.
