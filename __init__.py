#!/usr/bin/env/python

"""
MIRI Tools and Environment Software
===================================

Package MiriTE contains the tools, utilities, data models and common analysis
functions developed by the MIRI Software Team

The entire package can be built by executing the setup.py
script in the top level directory with the "install" argument.

Available packages
------------------
miri.tools
    Contained in MiriTE/tools
    Common utilities for MIRI software.

miri.datamodels
    Contained in MiriTE/datamodels
    Common data models for MIRI software.

miri.simulators
    Contained in MiriTE/simulators.
    MIRI simulator packages.
    
24 Jun 2010: Created.
27 Apr 2011: Add PACKAGES list
28 Nov 2013: Major simplification. Removed PACKAGES list and automatic
             test execution.
28 Sep 2016: Removed sandbox and stpipeline packages. Renamed
             miri.miritools.dataproduct to miri.datamodels.
             Renamed miri.miritools to miri.tools.
11 May 2017: Major update to new setup.py. Version number changed.
12 May 2017: MIRI software split into 4 pieces: MiriTools, MiriCalibration,
             MiriPipeline and MiriSimulators (V4.0)
19 Jun 2017: MIRI software merged with build 7.1 utilities (V5.0)
15 Jan 2018: MiriTools and MiriSimulators levels removed from package

@author: Steven Beard (UKATC)

"""
__project__ = 'MIRI Tool and Environment Software'
__author__ = 'MIRI Software Team'
__maintainer__ = __author__
__copyright__ = '2018, %s' % __author__
__version__ = '7.0'
