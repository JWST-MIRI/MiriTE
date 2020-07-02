#!/usr/bin/env python

"""
miri.tools
=========
Package miri.tools provides common utilities for the MIRI software.

Subpackages
-----------
dataproduct:
    Contains modules for the MIRI data products

measurements:
    Contains modules related to laboratory measurements and data sheets

Available modules
-----------------
filesearching:
    File searching utilities.

miriplot:
    Plotting utilities.

Scripts
-------

Data
----

:History:
29 Nov 2010: Add test section
24 Jan 2011: Added miridata
03 Feb 2011: Added mirikeyword
15 Mar 2011: Added quantum_efficiency
06 Sep 2011: miridata split into metadata and miricombdata
15 Sep 2011: Typo in miricombdata corrected
24 Oct 2011: filesearching module added.
16 Jan 2012: Test section comments added.
15 Mar 2012: Data product modules updated. Changed to V0.2dev.
02 Apr 2012: All module imports removed. The imports turned this file
             into a single point of failure, which meant that if one
             module broke all modules broke. Furthermore, one broken
             module could prevent setup.py working and stop the code
             being rebuilt. Having this file import everything also
             made tracking down the scipy/matplotlib interdependence
             bug very difficult.
06 Sep 2016: Changed to miri.tools
19 Jun 2017: Updated to use build 7.1 data models (V7.0)
05 Jan 2018: More version control information added. SVN info dropped.
10 Jun 2019: JWST build 7.3 release (v7.3)

"""
__project__ = 'MIRI Tools Software'
__author__ = 'MIRI Software Team'
__maintainer__ = 'MIRI Software Team: mirisim@roe.ac.uk'
__copyright__ = '2020, %s' % __author__
