Calibration Data Product Access (:mod:`datamodels.cdplib`)
==========================================================

.. module:: miri.datamodels.cdplib

Description
~~~~~~~~~~~
This module contains the MiriCDPInterface class and the global function 
used to access the MIRI Calibration Data Products (CDPs). CDP files are 
copied from the MIRI sftp repository to a local cache.

Objects
~~~~~~~
.. autoclass:: MiriCDPInterface
   :members:

Functions
~~~~~~~~~
.. autofunction:: get_cdp
.. autofunction:: cdp_version_decode

Environment Variables
~~~~~~~~~~~~~~~~~~~~~
CDP_DIR - Expected to contain the location of the local cache directory.

Global Data
~~~~~~~~~~~
CDP_DICT - A dictionary translating CDP data types into classes.
