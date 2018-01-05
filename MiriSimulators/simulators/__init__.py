#!/usr/bin/env python

"""

simulators
==========
Package simulators contains a selection of MIRI simulators.

Available subpackages
---------------------
scasim:
    MIRI sensor chip assembly simulator

mirimsim:
    MIRI imager simulator

Available modules
-----------------
integrators
    PoissonIntegrator class - General purpose particle counter with
        Poisson statistics.
    ImperfectIntegrator class - As above but with added latency effects.

find_simulator_file
    Search file system for matching files.
    
Scripts
-------
run_make_filters_fits.sh - Make all filter files used by simulators
run_make_qe_fits.sh - Make all QE files used by simulators

Data
----
data/amplifiers/*.fits
    Amplifier description files.
data/detector/*.fits
    Detector description files.

:History:
16 Jan 2012: Test section comments added.
02 Apr 2012: All module imports removed.
04 May 2013: Imported PoissonIntegrator
05 Jun 2013: Moved description of top level modules from scasim to here.
05 Jan 2018: More version control information added. SVN info dropped.

"""

# If you update this file don't forget also to update defsetup.py and 
# doc/source/conf.py. The release version number is defined in three 
# places. doc/source/release_notes.rst might also need to be updated.
__project__ = 'MIRI Simulators'
__author__ = 'MIRI Software Team'
__maintainer__ = __author__
__copyright__ = '2018, %s' % __author__
__version__ = '4.0'

# Import classes that need to be accessible from the miri.simulators level.
from integrators import PoissonIntegrator, ImperfectIntegrator
