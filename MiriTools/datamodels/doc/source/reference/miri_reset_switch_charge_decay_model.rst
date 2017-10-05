RSCD Data (:mod:`datamodels.miri_reset_switch_charge_decay_model`)
==================================================================

.. module:: miri.datamodels.miri_reset_switch_charge_decay_model

Description
~~~~~~~~~~~
This module contains the MiriResetSwitchChargeDecayModel class, which describes a 
MIRI reset switch charge decay model.

MIRI data models depend on the STScI data model, found in the 
jwst.datamodels package.

Objects
~~~~~~~
.. autoclass:: MiriResetSwitchChargeDecayModel
   :members:

Functions
~~~~~~~~~
When stored as a FITS file, a RSCD CDP contains the following HDUs

* Primary - metadata
* RSCD - reset switch charge decay table
