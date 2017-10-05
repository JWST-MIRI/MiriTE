Installation utilities (:mod:`miri.setup_utils`)
================================================

.. module:: miri.setup_utils

Description
~~~~~~~~~~~
This module contains general purpose functions for installing
and checking the MIRI software.


Functions
~~~~~~~~~
.. autofunction:: check_pythonpath
.. autofunction:: check_deps
.. autofunction:: clean_install
.. autofunction:: rmtree
.. autofunction:: parse_defsetup
.. autofunction:: prefix_dir_list
.. autofunction:: file_copy
.. autofunction:: file_append
.. autofunction:: file_add_package_imports
.. autofunction:: file_script_list
.. autofunction:: file_collection

Scripts
~~~~~~~
The script miri_installation_check.py can be run to check
the installation of the MIRI software and report the
version numbers of the Python libraries installed.
