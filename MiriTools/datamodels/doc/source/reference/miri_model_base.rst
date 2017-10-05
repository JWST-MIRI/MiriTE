MIRI Base Data Module (:mod:`datamodels.miri_model_base`)
=========================================================

.. module:: miri.datamodels.miri_model_base

Description
~~~~~~~~~~~
This module contains the MiriDataModel class, which is the base class
for all the MIRI data models. The class builds on the functions
provided by the STScI data model class DataModel and adds extra functions
specific to MIRI. Convenience functions are also included which implement
some of the features of the old MIRI data products.

All MIRI data models depend on the STScI data model, found in the
jwst.datamodels package. See

https://aeon.stsci.edu/ssb/svn/jwst/trunk/jwst_lib/models/doc/source/devel/models.rst

Objects
~~~~~~~
.. autoclass:: MIRIExtension
   :members:

.. autoclass:: MiriDataModel
   :members:

Functions
~~~~~~~~~
.. autofunction:: get_exp_type

Data formats
~~~~~~~~~~~~
The data formats for the JWST and MIRI data models are described in YAML
schemas contained in the "schemas" subdirectory. The following schema
files are particularly vital, since they are common to all the data models:

   * core.schema.yaml - Contains metadata common to all JWST data models.

   * referencefile.schema.yaml - Contains additional metadata used for
     calibration reference files.

   * miri_metadata.schema.yaml - Contains MIRI-specific metadata definitions.

Metadata may be found within the .meta attribute of a data model.

