MIRI Simulator Software (:mod:`miri.simulators`)
================================================

:Release: |release|
:Date: |today|

Introduction
~~~~~~~~~~~~
This document describes the implementation details for the MIRI general 
simulators package (miri.simulators).

Note that most MIRI simulators may now be found in a separate MIRISim
package.

Modules
~~~~~~~
General purpose simulator modules may be found in the simulators/lib 
directory

.. module:: miri.simulators

.. toctree::
   :maxdepth: 1

   find_simulator_file
   integrators

Subpackage - scasim
~~~~~~~~~~~~~~~~~~~
The MIRI SCA simulator.

.. toctree::
   :maxdepth: 1

   scasim/scasim
