#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""

Module cosmic_ray_roperties - Defines the properties of the cosmic
rays expected when the JWST is at L2. The following sources predict
the cosmic ray flux for MIRI:

[1]: J.C. Pickel, NASA/CR 2003 212504, "Modeling Charge Collection
in Detector Arrays", June 2003

[2]: M. Ressler, MIRI DFM 308 04.02, "Analysis of the Proton Flux
on the MIRI Detector Arrays", 3 August 2009.

[3]: Massimo Robberto, JWST-STScI-001928, SM-12, "A library of
simulated cosmic ray events impacting JWST HgCdTe detectors",
2 December 2009.

[4]: Andras Gaspar, Interpixel capacitance correction of the MIRI
detector for cosmic ray hits, 2017?

[5]: George Rieke, MIRI science operations communications, 2016-2017.

:History:

25 Jul 2010: Created
02 Aug 2010: Added cosmic ray mode 'NONE' for no cosmic rays.
30 Aug 2010: Locate the configuration files using an absolute file path
             so that scasim may be run from any directory. Check for
             files contained in the JWST cr_sim_ramp_fit data directory.
27 Sep 2010: Python environment for windows verified.
06 Oct 2010: Changed to new cosmic ray library for SiAs detectors.
24 Oct 2011: Modified to use new "filesearching" module, which should
             make it easier for users to substitute their own configuration
             files (not quite so important for the cosmic ray smiulations).
25 Oct 2011: Corrected bug in energy distribution for NONE mode.
02 Sep 2013: Make more resilient against cr_ramp_fit package not being
             installed.
24 Feb 2014: Use find_cr_library function to find CR library data files
             rather than searching PYTHONPATH (which takes a long time when
             the software is installed using MIRICLE/ureka).
08 Sep 2015: Made compatible with Python 3.
09 Nov 2016: Added blurring parameter.
15 Mar 2017: Added energy map binning factor.
05 Apr 2017: Switch to the newer version of Massimo's cosmic ray library,
             for a 35 micron detector thickness and installed locally.
07 Jul 2017: Added energy scaling factor. Blurring factor adjusted to match
             more closely the power law of low energy events predicted by [1].
14 Jul 2017: Flux factor changed to 1.7 to take into account the count rate
             predicted by [2]. IPC library variant removed (since the IPC is
             convolved within SCASim). Switch back to the 470 micron cosmic
             ray library, which gives better track lengths, but scale the
             electron count by 0.2. When combined with a blurring factor of
             0.3. This combination should more closely match the effects
             predicted by George Rieke [5] (memo of 11 July 2017).
04 Aug 2017: Removed unprintable characters.
             
@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


import os
import warnings

from miri.tools.filesearching import find_file_in_path

# Get the path to the folder where the cosmic ray library files have been
# installed.
try:    
    import miri.simulators.data
    cr_sim_path = miri.simulators.data.__path__[0]
    global_cr_data_path = os.path.join(cr_sim_path, 'cosmic_rays')
except ImportError:
    warnings.warn("NOTE: miri.simulators.data package not installed.")
    cr_sim_path = ''
    global_cr_data_path = ''

def find_cr_library( fileprefix, cr_data_path=global_cr_data_path):
    """
    
    Find the folder containing files with the given named prefix, first
    in a given search path of directories and (failing that) within the
    PYTHONPATH.
    By default, this function will search for a file within the current
    directory the MIRI module  data directories and then try PYTHONPATH.
    
    This function is used when, rather than there being a single parameter
    file, there are a whole set of parameter files with a common prefix.
    
    :Parameters:
    
    fileprefix: str
        A prefix for the files to be located (for example
        ''CR_SiAs_470_SUNMIN_'', which will cause all files matching
        ''CR_SiAs_470_SUNMIN_*.*'' to be located).
    cr_data_path: str
        Top level folder where the CR library data files are installed.
        
    :Returns:
    
    fileprefixpath: str
        If matching files are found the full path of the folder containing
        the files plus the prefix is returned (e.g. ''/data/CR_SiAs_35_SUNMIN_'').
        Returns None if no files are found.

    """    
    search_path = os.pathsep + './data'
    if cr_data_path:
        search_path += os.pathsep + cr_data_path

    filename = fileprefix + '*.*'
    firsttry = find_file_in_path( filename, search_path=search_path,
                                  pathsep=os.pathsep, walkdir=True,
                                  path_only=True)
    if firsttry:
        return os.path.join(firsttry, fileprefix)
    else:
        return None


DEFAULT_CR_MODE = 'SOLAR_MIN'

# Cosmic ray flux in number/micron^2/sec incident on the MIRI detectors
# at Solar Minimum, Solar Maximum and during a Solar Flare.
# Figures from Table 1 of [3].
CR_FLUX = {}
CR_FLUX['NONE']        = 0.0
CR_FLUX['SOLAR_MIN']   = 4.8983E-8
CR_FLUX['SOLAR_MAX']   = 1.7783E-8
CR_FLUX['SOLAR_FLARE'] = 3.04683E-5
# The following parameter can be used to adjust the cosmic ray flux
# for testing and fine tuning.
# [2] measured 5 hits per second on IRAC detectors, which scales to
# 55 hits per second on the MIRI detectors (depending on shielding).
# [3] predicts 32 hits per second on the MIRI detectors at SOLAR_MIN,
# assuming some shielding.
# Multiplying the [3] predictions by 1.7 gives 54 hits per second.
CR_FLUX_MULTIPLIER = 1.7

# The following parameter scales the 21x21 pixel energy map to adjust for
# the difference between the cosmic ray trail lengths contained in the
# library and the expected lengths.
# The 35 micron library of [3] generates average track lengths of ~1.5 pixels
# The 470 micron library of [3] generates average track lengths of ~6 pixels.
# [4] predicts average track lengths of 4.3 pixels and [5] suggests an
# average of 3.5 pixels.
# The 470 micron library binned by 2 gives average track lengths of ~3.5 pixels.
CR_BINNING_FACTOR = 2  # For best results, choose (1, 2, 3 or 7)
                       # 1 gives track lengths of 1-8 pixels
                       # 2 gives track lengths of 1-4 pixels
                       # 3 gives track lengths of 1-3 pixels
                       # 7 gives track lengths of 1 pixel

#
# NOTE: If SCASim is extended to cover other instruments, choose CR
# library file by detector type. At the moment it is fixed at SiAs.
#
CR_LIBRARY_FILES = {}
CR_LIBRARY_FILES['NONE'] = None
# Library files for the NIRCAM and NIRSPEC detectors (thickness 5.5 microns).
#CR_LIBRARY_FILES['SOLAR_MIN']   = find_cr_library('CRs_MCD5.5_SUNMIN')
#CR_LIBRARY_FILES['SOLAR_MAX']   = find_cr_library('CRs_MCD5.5_SUNMAX')
#CR_LIBRARY_FILES['SOLAR_FLARE'] = find_cr_library('CRs_MCD5.5_FLARES')
# Library files for the MIRI detectors.
if cr_sim_path:
# The 470 micron libraries are chosen to give the desired track lengths.
#     CR_LIBRARY_FILES['SOLAR_MIN']   = find_cr_library('CRs_SiAs_35_SUNMIN_')
#     CR_LIBRARY_FILES['SOLAR_MAX']   = find_cr_library('CRs_SiAs_35_SUNMAX_')
#     CR_LIBRARY_FILES['SOLAR_FLARE'] = find_cr_library('CRs_SiAs_35_FLARES_')
    CR_LIBRARY_FILES['SOLAR_MIN']   = find_cr_library('CRs_SiAs_470_SUNMIN_')
    CR_LIBRARY_FILES['SOLAR_MAX']   = find_cr_library('CRs_SiAs_470_SUNMAX_')
    CR_LIBRARY_FILES['SOLAR_FLARE'] = find_cr_library('CRs_SiAs_470_FLARES_')
else:
    CR_LIBRARY_FILES['SOLAR_MIN']   = None
    CR_LIBRARY_FILES['SOLAR_MAX']   = None
    CR_LIBRARY_FILES['SOLAR_FLARE'] = None
CR_LIBRARY_FILES['MIN'] = 0
CR_LIBRARY_FILES['MAX'] = 9
# The following variants of these library files are available.
# (No variants are available at the moment.)
CR_LIBRARY_FILES['VARIANTS'] = ['']
#
# The following parameters can be used to adjust the energy or the number
# or number of electrons released by any given cosmic ray for testing and
# fine tuning.
CR_ENERGY_MULTIPLIER = 1.0
# The 470 micron library of [3] generates cosmic ray events of suitable length
# but overestimates the electron count per pixel by a factor of 10 compared
# with [1]. Scale the rate by 0.1 to correct.
CR_ELECTRON_MULTIPLIER = 0.1

CR_NUCLEONS = ('H', 'He', 'C', 'N', 'O', 'Fe', '??')

# A map of the leakage of CR flux from a detector pixel into adjacent
# pixels by capacitative coupling.
# [4] and [5] have found that adjacent pixels are affected by cosmic rays
# at a level of 3.07%. 
CR_COUPLING = ((0.0001, 0.0307, 0.0001), \
               (0.0307, 0.8772, 0.0307), \
               (0.0001, 0.0307, 0.0001))
# This is the capacitative coupling function for electrons liberated by photons
# CR_COUPLING = ((0.0,   0.015, 0.0), \
#                (0.015, 0.94,  0.015), \
#                (0.0,   0.015, 0.0))

# Blur each cosmic ray image by this sigma (in pixels) to simulate
# the additional inter-pixel scattering suggested by [4] and [5] and
# to reproduce the power law predicted by [1] more closely. A level of
# 0.3 reproduces the expectation that an average cosmic ray will pass
# through 4.3 pixels and affect 13.7 (this blurring results in a median
# track length of 4 pixels which affects 14 pixels).
CR_BLUR = 0.3

#
# If the library files are not available, this set of electron distributions
# can be used to generate random cosmic ray events independently (although
# not as realistically). The lists contain (energy,probability), where the
# probabilities are incremental and normalised to a peak of 1.0. (They are
# converted to a cumulative distribution internally by the software.)
#
CR_ELECTRONS = {}
CR_ELECTRONS['NONE'] = \
    (
     (0.0,   1.0),
     (1.0e1, 0.0),
     (1.0e2, 0.0),
     (1.0e3, 0.0),
     (1.0e4, 0.0),
     (1.0e5, 0.0),
     (1.0e6, 0.0),
     (1.0e7, 0.0),
     )
CR_ELECTRONS['SOLAR_MIN'] = \
    (
     (1.0e1, 0.0),
     (1.8e1, 0.0),
     (3.2e1, 0.0),
     (5.6e1, 0.02),
     (1.0e2, 0.05),
     (1.8e2, 0.05),
     (3.2e2, 0.05),
     (5.6e2, 0.07),
     (1.0e3, 0.08),
     (1.8e3, 0.10),
     (3.2e3, 0.20),
     (5.6e3, 0.3),
     (1.0e4, 0.5),
     (1.8e4, 1.0),
     (3.2e4, 0.6),
     (5.6e4, 0.15),
     (1.0e5, 0.08),
     (1.8e5, 0.05),
     (3.2e5, 0.03),
     (5.6e5, 0.01),
     (1.0e6, 0.01),
     (1.8e6, 0.003),
     (3.2e6, 0.001),
     (5.6e6, 0.0),
     (1.0e7, 0.0),
     (1.8e7, 0.0),
     (3.2e7, 0.0),
     (5.6e7, 0.0),
     )
CR_ELECTRONS['SOLAR_MAX'] = \
    (
     (1.0e1, 0.0),
     (1.8e1, 0.0),
     (3.2e1, 0.0),
     (5.6e1, 0.01),
     (1.0e2, 0.02),
     (1.8e2, 0.02),
     (3.2e2, 0.02),
     (5.6e2, 0.03),
     (1.0e3, 0.05),
     (1.8e3, 0.07),
     (3.2e3, 0.15),
     (5.6e3, 0.4),
     (1.0e4, 1.0),
     (1.8e4, 0.25),
     (3.2e4, 0.1),
     (5.6e4, 0.05),
     (1.0e5, 0.05),
     (1.8e5, 0.05),
     (3.2e5, 0.05),
     (5.6e5, 0.01),
     (1.0e6, 0.01),
     (1.8e6, 0.005),
     (3.2e6, 0.001),
     (5.6e6, 0.0),
     (1.0e7, 0.0),
     (1.8e7, 0.0),
     (3.2e7, 0.0),
     (5.6e7, 0.0),
     )
CR_ELECTRONS['SOLAR_FLARE'] = \
    (
     (1.0e2, 0.0),
     (1.8e2, 0.001),
     (3.2e2, 0.001),
     (5.6e2, 0.03),
     (1.0e3, 0.05),
     (1.8e3, 0.06),
     (3.2e3, 0.06),
     (5.6e3, 0.07),
     (1.0e4, 0.1),
     (1.8e4, 0.12),
     (3.2e4, 0.25),
     (5.6e4, 0.6),
     (1.0e5, 1.0),
     (1.8e5, 0.6),
     (3.2e5, 0.5),
     (5.6e5, 0.37),
     (1.0e6, 0.18),
     (1.8e6, 0.07),
     (3.2e6, 0.04),
     (5.6e6, 0.03),
     (1.0e7, 0.02),
     (1.8e7, 0.01),
     (3.2e7, 0.0),
     (5.6e7, 0.0),
     )


if __name__ == '__main__':
    print( "NOTE: The CosmicRayProperties module is supposed to be " \
        "imported by another module, not run as a main program." )
    print( "The following cosmic ray properties are defined:" )
    print( "CR_FLUX\n-------" )
    for key in CR_FLUX:
        print( "%16s = %s" % (key, CR_FLUX[key]) )
    print( "CR_LIBRARY_FILES\n----------------" )
    for key in CR_LIBRARY_FILES:
        print( "%16s = %s" % (key, CR_LIBRARY_FILES[key]) )
    print( "CR_NUCLEONS\n--------------" )
    print( CR_NUCLEONS )
    print( "CR_COUPLING\n--------------" )
    for row in CR_COUPLING:
        print( "(%.3f, %.3f, %.3f)" % row )
    print( "CR_ELECTRONS\n--------------" )
    for key in CR_ELECTRONS:
        for row in CR_ELECTRONS[key]:
            str = key
            str += ": Relative probability of %.0f electrons is %f" % row
            print( str )
    print( "Finished." )
