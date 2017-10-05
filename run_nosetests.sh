#!/bin/bash
#
# If this script complains about the --with-xcoverage parameter
# not being recognised, try the command
#
#   conda install nosexcover
#
# Remove __init__.py and .pyc so the top level directory is no longer a package.
# Necessary to work around nosetests import problem.
\rm -f __init*
# Run nosetests with the default configuration parameters contained in setup.cfg
nosetests
#
# Reverse the removal of the __init__.py.
svn revert __init__.py
