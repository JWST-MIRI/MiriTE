MIRI Data Model Operations Module (:mod:`datamodels.operations`)
================================================================

.. module:: miri.datamodels.operations

Description
~~~~~~~~~~~
This module contains the HasMask, HasData and HasDataErrAndDq classes,
which provide functions for performing arithmetical operations on
data structures containing onr or more data arrays or bitwise operations
on data structures containing a mask array.

These classes are not meant to be used on their own. They are designed
to supplement the functionality of other data model classes.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

See also the description of the dqflags.py module.

Objects
~~~~~~~
.. autoclass:: HasMask
   :members:

.. autoclass:: HasData
   :members:

.. autoclass:: HasDataErrAndDq
   :members:

Functions
~~~~~~~~~
None
