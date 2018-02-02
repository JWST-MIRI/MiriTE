MIRI miritools release notes (:mod:`miri.datamodels`)
=====================================================

:Release: |release|
:Date: |today|

The MIRI datamodels package defines the data structures
used by MIRI data files and by MIRI calibration data
products. The MIRI data models are based on the JWST
data models package. Check the STScI portal for developer
documentation:

http://www.stsci.edu/portal

The data models support the input, management and output
of files in one of the standard JWST formats, as described
here in the user documentation:

https://jwst-docs.stsci.edu/display/JDAT/JWST+File+Names%2C+Formats%2C+and+Data+Structures

The MIRI data models require the installation of the
following packages:

jwst 0.7.8 or later;
adsf 1.2.2 or later;
astropy 3.0 or later;
yaml.

pysftp 0.2.8 or later;
paramiko 2.0.2 or later.

numpy 1.13.1 or later;
scipy 0.19.1 or later.

The following package is optional, and it will be
used when available:

matplotlib 2.0.2 or later.

These dependencies should be taken care of automatically
by the MIRICLE installation script.
