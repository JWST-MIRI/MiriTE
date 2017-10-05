MIRI Software |release| Installation Instructions
=================================================

:Release: |release|
:Date: |today|
:Author: Steven Beard <steven.beard@stfc.ac.uk>

Notes
~~~~~

* The source file for this page is available at: \
  https:https://aeon.stsci.edu/ssb/svn/jwst/trunk/teams/miri/doc/source/install.rst

* The recommended way of installing this software is using the MIRICLE
  script described at: \
  http://miri.ster.kuleuven.be/bin/view/Internal/MiricleSoftwareSystem
  
  The MIRICLE script can be downloaded from \
  http://miri.ster.kuleuven.be/pub/Internal/MiricleSoftwareSystem/MIRICLE_install.bash

* MIRICLE is designed for operating systems which understand the bash script.
  For installation on Windows workstations, or anything else without bash,
  the following MIRI project wiki pages provide alternative instructions:

  - Installation instructions for Python and stsci_python may be found at: \
    http://miri.ster.kuleuven.be/bin/view/Internal/Software/PythonInstall
  - Information about the configuration of the MIRI software may be found at: \
    http://miri.ster.kuleuven.be/bin/view/Internal/Software/PythonCommonCodeRepository
  - Information on the installation and use of SCASim may be found at: \
    http://miri.ster.kuleuven.be/bin/view/Internal/Software/ScaSim

Building the MIRI Software
--------------------------
The software is built by executing the command

python setup.py install

from the top level teams/miri directory.

Test the MIRI software installation
-----------------------------------
The "teams/miri/miri_installation_check.py" script can be run to check that
all the necessary modules have been installed.

The "teams/miri/miri_run_tests.py" script can be run to execute all the MIRI
module unit tests.
