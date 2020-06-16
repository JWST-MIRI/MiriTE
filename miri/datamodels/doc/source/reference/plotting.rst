MIRI Data Model Plotting Module (:mod:`datamodels.plotting`)
============================================================

.. module:: miri.datamodels.plotting

Description
~~~~~~~~~~~
This module contains the DataModelVisitor and DataModelPlotVisitor classes,
which provide model-independent plotting functions. These classes implement
the visitor design pattern, and are invoked by calling their visit() method
with the data object to be plotted passed as a parameter.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: DataModelVisitor
   :members:

.. autoclass:: DataModelPlotVisitor
   :members:

Functions
~~~~~~~~~
None
