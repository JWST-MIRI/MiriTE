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

data_maps
    LEGACY DATA MODELS.
    DataMap class - Generic data map manager.
    Illumination map class - Loads/manages/saves/displays the description
                             of the illumination on a detector.
    DarkMap class          - Loads/manages/saves/displays dark maps.
    BadPixelMap class      - Loads/manages/saves/displays bad pixel masks.

exposure_data
    LEGACY DATA MODELS.
    ExposureData class - Assembles/manages/saves/displays the data
    created by a simulation.

detector
    DetectorArray class - Describes a MIRI detector

amplifier
    Amplifier class - Describes a detector amplifier.

cosmic_ray
    CosmicRay class      - Describes a cosmic ray event.
    CosmicRayEnvironment - Describes the cosmic ray environment.
    
Scripts
-------
scasim - Run a SCASIM simulation.

make_sca_file - Make a detector illumination file.

make_sca_calibration - Make an artificial calibration file.

make_qe_fits - Define a detector quantum efficiency measurement.

make_measurements_fits - Define a detector parameter measurement.

make_bad_pixel_mask - Define a detector bad pixel mask.

make_maps_from_dhas - Convert a DHAS format bad pixel mask to SCASim files.

Data
----
data/SCATestInput80x64.fits
    An example FITS file in the format readable by the SCA simulator.
    All MIRI simulators need to describe their detector illumination
    in a file of this format.

:History:
16 Jan 2012: Test section comments added.
02 Apr 2012: All module imports removed.
13 Nov 2012: Description changed to reflect new teams/miri package structure.
05 Jun 2013: Moved description of top level modules to miri.simulators.
"""

# If you update this file don't forget also to update defsetup.py and 
# doc/source/conf.py. The release version number is defined in three 
# places. doc/source/release_notes.rst might also need to be updated.
__author__ = 'Steven Beard, MIRI Software Team'
__maintainer__ = 'Steven Beard <steven.beard@stfc.ac.uk>'
__version__ = '4.0'
__svn_info__ = '$Id:$'


