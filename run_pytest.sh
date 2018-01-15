#!/bin/bash
#
# Remove __init__.py and .pyc so the top level directory is no longer a package.
# Necessary to work around nosetests import problem.
cp __init__.py __init__.py_BACKUP
\rm -f __init*
# Run pytest with the default configuration parameters contained in setup.cfg
pytest
#
# Reverse the removal of the __init__.py.
git checkout -- __init__.py
