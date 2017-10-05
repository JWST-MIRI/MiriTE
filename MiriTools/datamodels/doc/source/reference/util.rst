General Utility Functions (:mod:`datamodels.util`)
==================================================

.. module:: miri.datamodels.util

Description
~~~~~~~~~~~
This module contains a collection of general purpose utility functions
used to access, query, convert and test the MIRI data models.

Objects
~~~~~~~
None.

Functions
~~~~~~~~~
.. autofunction:: open
.. autofunction:: read_fits_header
.. autofunction:: get_data_class
.. autofunction:: assert_recarray_equal
.. autofunction:: assert_products_equal
.. autofunction:: add_subarray_metadata
.. autofunction:: verify_metadata
.. autofunction:: verify_cdp_file
.. autofunction:: verify_fits_file
.. autofunction:: convert_detector
.. autofunction:: convert_band
.. autofunction:: convert_cdp_2to3

Global Data
~~~~~~~~~~~
MIRI_MODELS - List of recognised MIRI models
MIRI_DETECTORS - List of recognised MIRI detectors
MIRI_SETTINGS - List of recognised MIRI detector settings
MIRI_READPATTS - List of recognised MIRI readout patterns
MIRI_SUBARRAYS - List of recognised MIRI subarrays
MIRI_CHANNELS - List of recognised MIRI channels
MIRI_BANDS - List of recognised MIRI bands
MIRI_FILTERS - List of recognised MIRI filters
CDP_METADATA - Rules for verifying compulsory CDP metadata
CDP_SUBARRAY - Rules for verifying CDP subarray metadata
CDP_HISTORY - Rules for verifying CDP HISTORY metadata
