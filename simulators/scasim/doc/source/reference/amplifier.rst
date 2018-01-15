Amplifier module (:mod:`miri.simulators.scasim.amplifier`)
==========================================================

.. module:: miri.simulators.scasim.amplifier

Description
~~~~~~~~~~~
This module contains the Amplifier class, which simulates the behaviour
of the readout amplifiers contained in a MIRI focal plane module.

*Note: Other amplifiers, including those used by non-MIRI MULTIACCUM
detectors can be added by modifying amplifier_properties.py*

NOTE: THE AMPLIFIER CLASS HAS BEEN REMOVED FROM THE SIMULATOR TO SUPPORT
DETECTOR-WIDE CDP CALIBRATION FILES.

Objects
~~~~~~~
.. autoclass:: Amplifier
   :members:

Functions
~~~~~~~~~
.. autofunction:: set_reference_level

.. autofunction:: random_level_change

.. autofunction:: check_readout


Configuration files
~~~~~~~~~~~~~~~~~~~
The file amplifier_properties.py contains configuration information.

.. literalinclude:: amplifier_properties.py
