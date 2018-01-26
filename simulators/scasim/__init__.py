#!/usr/bin/env python

"""

simulators.scasim
=================
Subpackage simulators.scasim contains the MIRI Sensor Chip Assembly 
(SCA) simulator.
See the LICENCE file for terms and conditions.

Available modules
-----------------
sensor_chip_assembly
    SensorChipAssembly - Top level SCA simulator class

exposure_data
    LEGACY DATA MODELS.
    ExposureData class - Assembles/manages/saves/displays the data
    created by a simulation.

detector
    DetectorArray - Describes a MIRI detector

amplifier
    LEGACY CLASS - not used at the moment
    Amplifier - Describes a detector amplifier.

cosmic_ray
    CosmicRay            - Describes a cosmic ray event.
    CosmicRayEnvironment - Describes the cosmic ray environment.
    
Scripts
-------
scasim - Run an SCASIM simulation.

convert_exposure_data - Convert exposure data between FITSWRiter and STScI format

plot_exposure - Plot the contents of an exposure (or ramp) data file.

make_sca_file - Make a detector illumination file.

make_bad_pixel_mask - Define a detector bad pixel mask.

make_fringe_map - Make an artificial fringe map.

make_sca_calibration - Make an artificial calibration file.


Data
----
data/SCATestInput80x64.fits
    An example FITS file in the format readable by the SCA simulator.
    All MIRI simulators need to describe their detector illumination
    in a file of this format.

:History:
16 Jan 2012: Test section comments added.
02 Apr 2012: All module imports removed.
13 Nov 2012: Description changed to reflect new package structure.
05 Jun 2013: Moved description of top level modules to miri.simulators.
21 Nov 2017: Updated list of scripts.
05 Jan 2018: More version control information added. SVN info dropped.

"""

# If you update this file don't forget also to update defsetup.py and 
# doc/source/conf.py. The release version number is defined in three 
# places. doc/source/release_notes.rst might also need to be updated.
__project__ = 'MIRI SCA Simulator'
__author__ = 'MIRI Software Team'
__maintainer__ = 'MIRI Software Team: miri@xxx.yyy'
__copyright__ = '2018, %s' % __author__
