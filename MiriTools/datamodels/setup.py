#!/usr/bin/env python
#
# Build the package using the STScI distutils extensions.
#
# The default setup parameters are contained in the script defsetup.py.
#
from __future__ import division, print_function

import stsci.tools.stsci_distutils_hack as stsci_distutils_hack
stsci_distutils_hack.run()
# This is an alternative build if you want to check that a particular 
# version of pytools is being used. 
#stsci_distutils_hack.run(pytools_version = "3.1")
print("setup.py completed.")
 
