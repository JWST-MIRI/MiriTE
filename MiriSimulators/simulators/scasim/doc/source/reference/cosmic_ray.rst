Cosmic ray module (:mod:`miri.simulators.scasim.cosmic_ray`)
============================================================

.. module:: miri.simulators.scasim.cosmic_ray

Description
~~~~~~~~~~~
This module contains all the classes and functions for simulating cosmic
ray events. The module contains two classes. The CosmicRayEnvironment
class simulates the cosmic ray environment and is programmed with
parameters such as the expected cosmic ray flux and energy distribution.
The simulation uses a library of cosmic ray events published Robberto,
[2]. Scale factors are applied to this simulation to match the more
recent analyses of Rieke, [5], and Gaspar, [4], which are more
consistent with Pickel [1] and Ressler [2].

If the cosmic ray library is not available, the CosmicRayEnvironment
class can also generate cosmic ray events at random. Each cosmic ray
event is described using the CosmicRay class.

Note that the frequency of cosmic ray events is determined both by the
properties of the cosmic rays, such as the expected flux, and by the
size of the detector pixels and amplifiers, which determines the target
area. The thicker pixels of the MIRI detectors produce slightly longer
cosmic ray trails than the NIRCAM and NIRSPEC detectors.

References
~~~~~~~~~~
[1]: J.C. Pickel, NASA/CR 2003 212504, "Modelling Charge Collection
in Detector Arrays", June 2003

[2]: M. Ressler, MIRI DFM 308 04.02, "Analysis of the Proton Flux
on the MIRI Detector Arrays", 3 August 2009.

[3]: Massimo Robberto, JWST-STScI-001928, SM-12, "A library of
simulated cosmic ray events impacting JWST HgCdTe detectors",
2 December 2009.

[4]: Andras Gaspar, Interpixel capacitance correction of the MIRI
detector for cosmic ray hits, 2017?

[5]: George Rieke, MIRI science operations communications, 2016-2017.

Objects
~~~~~~~
.. autoclass:: CosmicRay
   :members:

.. autoclass:: CosmicRayEnvironment
   :members:

Functions
~~~~~~~~~
.. autofunction:: load_cosmic_ray_random

.. autofunction:: load_cosmic_ray_library

.. autofunction:: plot_cosmic_ray_events

Configuration files
~~~~~~~~~~~~~~~~~~~
The file cosmic_ray_properties.py contains configuration information.

.. literalinclude:: cosmic_ray_properties.py
