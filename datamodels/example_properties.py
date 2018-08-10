#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module example_properties - an example MIRI parameters file.

:History:

25 Oct 2011: Created
28 Mar 2012: Recursion comment added.
08 Sep 2015: Made compatible with Python 3.
01 Jul 2018: Final conversion to Python 3.

@author: Steven Beard (UKATC)

"""

# The file searching function, find_file, can be used to find a parameter
# file with the specified name somewhere within the tree of directories
# within the current working directory, in the module data directories or
# in the PYTHONPATH.
# BEWARE! This module is included as an example by the __main__ of
# filesearching.py. The __main__ test within filesearching.py prevents an
# infinite recursion.
from miri.tools.filesearching import find_file

# Define a few simple scalar parameters
INTEGER_PARAMETER = 42
FLOAT_PARAMETER = 3.14159265
STRING_PARAMETER = 'This is a string'

# Define a tuple parameter
TUPLE_PARAMETER = (1,2,3,4,5,6,7,8,9,10)

# Define a dictionary parameter
DICT_PARAM = {'one': 1,
              'two': 2,
              'three': 3}

#Define a calibration file parameter, located using the find_file function.
FILE_PARAM = find_file('example_qe.fits')

# Now define set of parameters contained in a nested dictionary
_detector1 = {'NAME': 'Holmes', 'ROWS': 1024, 'COLUMNS': 1024}
_detector2 = {'NAME': 'Marple', 'ROWS': 2048, 'COLUMNS': 2048}
_detector3 = {'NAME': 'Maigret', 'ROWS': 512, 'COLUMNS': 512}
_detector4 = {'NAME': 'Frost', 'ROWS': 1600, 'COLUMNS': 1200}

DETECTORS_DICT = {'Holmes':  _detector1,
                  'Marple':  _detector2,
                  'Maigret': _detector3,
                  'Frost':   _detector4}

