#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module find_simulator_file - Contains function to search for simulator files.


:History:

24 Feb 2014: Created
08 Sep 2015: Made compatible with Python 3.
23 Mar 2016: Added a Python logger so a message can be displayed if a
             simulator file is not found.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

import os

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

from miri.tools.filesearching import find_file_in_path

def find_simulator_file( filename, logger=LOGGER ):
    """
    
    Find a file matching the given name first in a given search path of
    directories and (failing that) within the PYTHONPATH.
    By default, this function will search for a file within the current
    directory the MIRI module data directories and then try PYTHONPATH.
    
    :Parameters:
    
    filename: string
        The name of the file to be located.
    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
        
    :Returns:
    
    filepath: string
        The full path and name of the matching file.
        If no file is found an empty string is returned.    
    
    """
    import miri.simulators.data
    import miri.simulators.scasim.data
    import miri.simulators.mirimsim.data
    
    search_path = os.pathsep + './data'
    search_path += os.pathsep + miri.simulators.data.__path__[0]
    search_path += os.pathsep + miri.simulators.scasim.data.__path__[0]
    search_path += os.pathsep + miri.simulators.mirimsim.data.__path__[0]
    
    firsttry = find_file_in_path( filename, search_path=search_path,
                                  pathsep=os.pathsep, walkdir=True)
    if firsttry:
        logger.debug("Found simulator file: %s\n\tat %s" % (filename, firsttry))
        return firsttry
    else:
        logger.error("***Could not find simulator file: %s" % filename)
        return ''
        
