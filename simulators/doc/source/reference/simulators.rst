MIRI Simulator Software (:mod:`miri.simulators`)
================================================

:Release: |release|
:Date: |today|

Introduction
~~~~~~~~~~~~
This document describes the implementation details for the MIRI general 
simulators package (simulators).

Modules
~~~~~~~
General purpose simulator modules may be found in the simulators/lib 
directory

.. module:: miri.simulators

.. toctree::
   :maxdepth: 1

   find_simulator_file
   integrators

Subpackage - mirimsim
~~~~~~~~~~~~~~~~~~~~~
The MIRI imager simulator. NOTE: This subpackage is empty. mirisimsim has
been moved to the MIRISim package.

.. toctree::
   :maxdepth: 2

   mirimsim/mirimsim

Subpackage - scasim
~~~~~~~~~~~~~~~~~~~
The MIRI SCA simulator.

.. toctree::
   :maxdepth: 1

   scasim/scasim
