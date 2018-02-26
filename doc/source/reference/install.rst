MIRI TE Software |release| Installation Instructions
====================================================

:Release: |release|
:Date: |today|
:Author: MIRI Software Team <mirisim@roe.ac.uk>

Notes
~~~~~

* The source file for this page is available at: \
  https://github.com/JWST-MIRI/MiriTE/tree/master/doc/source/install.rst

* For more information on installing this software, see: \
  https://github.com/JWST-MIRI/MiriTE
  
* The installation script is designed for operating systems which understand
  the bash script.

Building the MIRI Software
--------------------------
The software is built by executing the command

python setup.py install

from the top level directory.

Test the MIRI software installation
-----------------------------------
The "miri_installation_check.py" script can be run to check that
all the necessary modules have been installed.

The "run_nosetests.sh" or "run_pytests.sh" scripts can be run to execute
all the MIRI module unit tests. "python setup.py nosetests" should also
work.
