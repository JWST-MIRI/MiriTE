Fringe Frequencies Data (:mod:`datamodels.miri_fringe_frequencies_model`)
=========================================================================

.. module:: miri.datamodels.miri_fringe_frequencies_model

Description
~~~~~~~~~~~
This module contains the MiriMrsFringeFrequenciesModel class, which 
describes the MIRI fringe frequencies model

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriMrsFringeFrequenciesModel
   :members:

Functions
~~~~~~~~~
None

Data formats
~~~~~~~~~~~~
When stored as a FITS file, a fringe frequencies CDP contains the
following HDUs

* Primary - metadata
* FRINGE_FREQ - fringe frequencies table
