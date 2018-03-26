MIRI simulators release notes  (:mod:`miri.simulators`)
=======================================================

:Release: |release|
:Date: |today|

General notes
-------------
The MiriTE simulators package was created historically
to contain a selection of MIRI simulators: one simulator
for each instrument mode (imager, LRS, MRS) and a
separate simulator for the MIRI detectors (SCASIm).

In practice, the MIRI instrument modes are now simulated
by a separate package called MIRiSim, so this package
contains just the detector simulator, SCASIm. MiriTE and
MIRISim are installed together to provide a fully-functional
simulator.
See http://miri.ster.kuleuven.be/bin/view/Public/MIRISim_Public

The MIRI simulators package makes use of the MIRI
tools and MIRI data models, and therefore has the same
dependencies as the miri.tools and miri.datamodels
packages.

The dependencies should be taken care of automatically
by the MIRICLE installation script.
See http://miri.ster.kuleuven.be/bin/view/Public/MirisimInstallation


Module specific release notes
-----------------------------

.. toctree::
   :maxdepth: 1

   scasim/release_notes.rst
