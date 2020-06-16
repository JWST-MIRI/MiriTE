MIRI Data Quality Flags Module (:mod:`datamodels.dqflags`)
===========================================================

.. module:: miri.datamodels.dqflags

Description
~~~~~~~~~~~
This module contains the FlagsTable class, together with the functions
for managing and combining data quality flags.

Objects
~~~~~~~
.. autoclass:: FlagsTable
   :members:

Functions
~~~~~~~~~
.. autofunction:: convert_dq
.. autofunction:: flags_table_to_metadata
.. autofunction:: make_mask
.. autofunction:: format_mask
.. autofunction:: raise_mask
.. autofunction:: lower_mask
.. autofunction:: test_mask_all
.. autofunction:: test_mask_any
.. autofunction:: test_mask_none
.. autofunction:: mask_array_all
.. autofunction:: mask_array_any
.. autofunction:: mask_array_none
.. autofunction:: combine_quality

Global Data
~~~~~~~~~~~
master_flags - Defines the JWST pipeline master flags table
groupdq_flags - Defines the flags used to describe the GROUPDQ table
pixeldq_flags - Defines the flags used to desribe the PIXELDQ table

