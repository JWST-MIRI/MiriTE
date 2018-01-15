Detector module (:mod:`miri.simulators.scasim.detector`)
========================================================

.. module:: miri.simulators.scasim.detector

Description
~~~~~~~~~~~
The detector module contains the DetectorArray class, which (together
with the amplifier module) simulates the behaviour of a MIRI sensor
chip assembly/detector combination (MIRIMAGE, MIRIUFULONG or
MIRIFUSHORT) with focal plane module (FPM) 106, 104 or 105.

*Note: Other MULTIACCUM detectors (including non-MIRI detectors) can
be added by modifying detector_properties.py*

Objects
~~~~~~~
.. autoclass:: DetectorArray
   :members:

Functions
~~~~~~~~~
None.

Configuration files
~~~~~~~~~~~~~~~~~~~
The file detector_properties.py contains configuration information.

.. literalinclude:: detector_properties.py
