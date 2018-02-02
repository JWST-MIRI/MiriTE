MIRI Data Models (:mod:`miri.datamodels`)
=========================================

Release Notes
~~~~~~~~~~~~~
.. toctree::
   :maxdepth: 1

   release_notes
   licence

Introduction
~~~~~~~~~~~~
This document describes the implementation details for the MIRI data
models. These data models are built on the underlying JWST data models
package. See the JWST user documentation:

https://jwst-docs.stsci.edu/display/HOM/JWST+User+Documentation+Home

Modules
~~~~~~~
General purpose miri.datamodels modules may be found in the miri/datamodels/lib 
directory

.. module:: miri.datamodels

.. toctree::
   :maxdepth: 1

   miri_model_base
   miri_measured_model
   miri_aperture_correction_model
   miri_badpixel_model
   miri_dark_reference_model
   miri_distortion_models
   miri_droop_model
   miri_flatfield_model
   miri_fluxconversion_models
   miri_fringe_frequencies_model
   miri_gain_model
   miri_ipc_model
   miri_jump_model
   miri_lastframe_model
   miri_latent_model
   miri_linearity_model
   miri_pce_model
   miri_photometric_models
   miri_pixel_saturation_model
   miri_psf_models
   miri_readnoise_model
   miri_reset_model
   miri_reset_switch_charge_decay_model
   miri_spectral_spatial_resolution_model
   miri_straylight_model
   miri_telescope_emission_model
   miri_transmission_correction_model
   miri_wavelength_correction_model
   miri_illumination_model
   miri_exposure_model
   miri_filters
   miri_measurement
   cdplib
   operations
   util
   compatibility
   dqflags
   plotting

Unit tests corresponding to these modules may be found in the 
miri/datamodels/tests directory.

Scripts
~~~~~~~
General purpose miri.datamodels scripts may be found in the miri/datamodels/scripts 
directory.

Add History Metadata (cdp_add_history)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_add_history
   :members:

Add Subarray Metadata (cdp_add_subarray)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_add_subarray
   :members:

Get documentation for CDP File (cdp_get_doc)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_get_doc
   :members:

Print Contents of CDP File (cdp_print)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_print
   :members:

Extract Subset from DARK CDP (cdp_reduce_dark)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_reduce_dark
   :members:

Verify CDP Against Requirements (cdp_verify)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: cdp_verify
   :members:

Verify Many CDP Files (multicdp_verify)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: multicdp_verify
   :members:

Find More Recent CDP (find_me_another)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: find_me_another
   :members:

Make Filter Data File (make_filters_fits)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: make_filters_fits
   :members:

Make QE Data File (make_qe_fits)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: make_qe_fits
   :members:

Make Measurement Data File (make_measurement_fits)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: make_measurements_fits
   :members:

Data
~~~~
General purpose miri.datamodels installation data may be found
in the miri/datamodels/data directory.

Example Data
^^^^^^^^^^^^
:mod:`miri.datamodels` comes with some example filter and
measurement data, provided as ASCII and FITS files (see
:mod:`miri.datamodels.miri_filters` and 
:mod:`miri.datamodels.miri_measurement`).

Configuration Data
^^^^^^^^^^^^^^^^^^
An example measurement properties file (for 
:mod:`miri.datamodels.miri_measurement`) may be found in 
miri/datamodels/lib/example_properties.py.

