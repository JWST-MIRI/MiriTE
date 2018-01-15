MIRI Data Model Compatibility Module (:mod:`datamodels.compatibility`)
======================================================================

.. module:: miri.datamodels.compatibility

Description
~~~~~~~~~~~
This module contains the WithCompatibility class, which can be used to
help convert existing software using the old MIRI data products to use
the new data products. The class provides some bolt-on attributes and
functions which mimic some features of the old model which are not
included in the MiriDataModel class.

NOTE: There are not many functions in this module. Long term use of
this module is not recommended. It should be used temporarily
while converting existing software. The source code shows how the
old functionality can be implemented in the new data model.

Objects
~~~~~~~
.. autoclass:: WithCompatibility
   :members:

Functions
~~~~~~~~~
None
