MIRI Sensor Chip Assembly Simulator (:mod:`miri.simulators.scasim`)
===================================================================

Release Notes
~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 1

   release_notes
   licence

Introduction
~~~~~~~~~~~~
This document describes the implementation details for the MIRI Sensor Chip
Assembly (SCA) simulator software. It should be read alongside the "MIRI
Sensor Chip Assembly (SCA) Simulator User Manual" (1), which describes how
to use the software, and the "MIRI Sensor Chip Assembly (SCA) Simulator
Software Design Document" (2), which describes the design on which this
software implementation is based.

The SCA simulator software is designed to be run from the scripts - e.g.
see the description of the scasim script in "Scripts" section below. The
scripts are given basic information, such as the input and output file
names, the detector temperature and cosmic ray environment to be simulated
and the detector readout and subarray modes required. It is also possible
to call the SCA simulator directly from Python software - for example see
the description of the simulate_sca_pipeline function within the
"sensor_chip_assembly" module. More detailed information about the
properties of the SCA (such as how the dark current and read noise changes
with temperature) is contained within configuration files - see the "Data"
section for details. The "Modules" section describes the detailed contents
of the software modules.

The software depends on the tools and utilities in the miri.tools module
and the data models defined in the miri.datamodels module.

References
^^^^^^^^^^

1. MIRI Sensor Chip Assembly (SCA) Simulator User Manual, V4.0,
   4 August 2017.

2. MIRI Sensor Chip Assembly (SCA) Simulator Software Design Document, V4.0,
   4 August 2017.

3. JPL D-46944, MIRI Flight Focal Plane Module End Item
   Data Package (FPM EIDP) (edited version), A. Schneider et al.,
   Initial X1, 10 May 2010

4. The Mid-Infrared Instrument for the James Webb Space Telescope,
   I: Introduction; G.H.Rieki et al, Publications of the Astronomical
   Society of Pacific, Volume 127, Issue 953, pp. 584 (2015)
   
5. The Mid-Infrared Instrument for the James Webb Space Telescope,
   VII: The MIRI Detectors; G.H.Rieki et al., Publications of the
   Astronomical Society of Pacific, Volume 127, Issue 953, pp. 665
   (2015)

6. JPL MIRI DFM 478 04.02, MIRI FPS Exposure Time Calculations (SCE FPGA2),
   M. E. Ressler, October 2014.

7. The Mid-Infrared Instrument for the James Webb Space Telescope,
   VIII: The MIRI Focal Plane System; M. E. Ressler et al.,
   Publications of the Astronomical Society of Pacific, Volume 127,
   Issue 953, pp. 675 (2015)

8. MIRI DFM 308-04.02, Analysis of the Proton Flux on the MIRI Detector Arrays,
   M. Ressler, 3 August 2009.

9. JWST-STScI-00198, SM-12, A Library of Simulated Cosmic Ray Events Impacting
   JWST HgCdTe Detectors, M. Robberto, 9 March 2010.

10. JWST Calibration Software Online Documentation,
    https://jwst-docs.stsci.edu/display/JDAT 

11. DHAS miri_sloper Online Documentation, http://tiamat.as.arizona.edu/dhas/
    
12. MIRI CDP Calibration Data Products documentation,
    http://miri.ster.kuleuven.be/bin/view/Internal/CalDataProducts

Scripts
~~~~~~~
The scasim scripts may be found in the scasim/scripts directory.

scasim
^^^^^^
.. automodule:: scasim
   :members:

detector_latency_test
^^^^^^^^^^^^^^^^^^^^^
.. automodule:: detector_latency_test
   :members:

convert_exposure_data
^^^^^^^^^^^^^^^^^^^^^
.. automodule:: convert_exposure_data
   :members:

make_sca_file
^^^^^^^^^^^^^
.. automodule:: make_sca_file
   :members:

make_bad_pixel_mask
^^^^^^^^^^^^^^^^^^^
.. automodule:: make_bad_pixel_mask
   :members:

make_fringe_map
^^^^^^^^^^^^^^^^^^^
.. automodule:: make_fringe_map
   :members:

make_sca_calibration
^^^^^^^^^^^^^^^^^^^^
.. automodule:: make_sca_calibration
   :members:

Configuration Information
~~~~~~~~~~~~~~~~~~~~~~~~~
The parameters provided to the scripts, described in the
previous section, are used to control those aspects of a
MIRI simulation that are expected to change from observation
to observation. More detailed control of the simulations is
possible by editing the configuration files supplied with
the simulator. The files, which may be found in the
simulators/scasim directory, are::

    cosmic_ray_properties.py
    detector_properties.py
    
These files describe the cosmic ray environment,
and the properties of the detector and electronics.

To make a change to one or more of these files, copy them to
your current working directory and edit them. As long as the
new files are within your current working directory (or a
subdirectory within it), there is no need to rebuild the
SCAsim software to make the new parameters current. To revert
to the original parameters, simple delete or rename the edited
files.

SCASim will search the current directory, the MIRI data
directories and then the PYTHONPATH directories for these files,
and it will report the full path and name of each one used.
Make sure SCAsim is using the files you expect.

*NOTE: These configuration files are processed by the Python
interpreter as if they were source code. Be careful not to miss
out matching quotes or brackets. If there is a problem in
interpreting one of these files, you will see a parsing error
when SCAsim is executed. The offending file will be named.
See the "MIRI Sensor Chip Assembly (SCA) Simulator User Manual"
for details.*

Shared Modules
~~~~~~~~~~~~~~
The following modules may be found in the top level simulators
package, and are shared with other MIRI simulators.

*  find_simulator_file
*  integrators

.. module:: miri.simulators

.. toctree::
   :maxdepth: 1

   find_simulator_file
   integrators

Modules
~~~~~~~
The miri.simulators.scasim modules may be found in the simulators/scasim directory.

.. module:: miri.simulators.scasim

.. toctree::
   :maxdepth: 1

   sensor_chip_assembly
   detector
   cosmic_ray
   exposure_data

Measurement Data
~~~~~~~~~~~~~~~~
:mod:`miri.scasim` comes with configuration data describing the cosmic ray,
detector and amplifier properties, plus average dark current and
quantum efficiency measurements. Other calibration data are extracted from
the MIRI Calibration Data Products (CDP) repository. The onboard
measurement data can be provided in ASCII or FITS format files. These
files are normally stored in the scasim/data directory, but if SCASim
detects files with the same name stored in the current working directory
it will use those in preference. The ASCII files are useful for storing
measurement data cut and pasted from elsewhere, whereas the FITS files are
more useful because they can store ancilliary information (such as
measurement comments) in their header.
The 'miri.datamodels.make_measurements_fits' and 'make_qe_fits' scripts may be
used to convert the ASCII format data into FITS. The make_bad_pixel_mask
script can be used to create a tailored bad pixel mask in FITS format.
Test calibration files (including bad pixel maps and dark maps) can be made
with the make_sca_calibration script. Some examples are given below:

*NOTE: SCAsim is now configured by default to obtain bad pixel data,
dark current variation, read noise data, pixel flat-field data and
amplifier gain data from the MIRI CDP repository.*

./data/detector/dark_currentXXX.fits
   A FITS file containing a dark current variation map (for Sensor Chip
   Assembly XXX).

./data/detector/qe_measurementXXX.txt
   An ASCII file containing a detector quantum efficiency measurement
   (for Sensor Chip Assembly XXX). The equivalent FITS file has the
   same name plus a '.fits' extension.

./data/SCATestInput80x64.fits
   An example SCA-format FITS file containing an 80x64 pixel detector
   illumination image.

Dark Current vs Temperature Measurement Example (ASCII format)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: dark_currentIM.txt

Quantum Efficiency vs Wavelength Measurement Example (ASCII format)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: qe_measurementLW.txt

Bad Pixel Specification Example (ASCII format)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: bad_pixel.txt
