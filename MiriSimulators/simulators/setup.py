#!/usr/bin/env python
#
# Build the package using the STScI distutils extensions.
#
# The default setup parameters are contained in the script defsetup.py.
#
from __future__ import division, print_function

import stsci.tools.stsci_distutils_hack as stsci_distutils_hack
stsci_distutils_hack.run()
print("setup.py completed.")
