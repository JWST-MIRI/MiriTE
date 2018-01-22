MIRI Data Models (:mod:`miri.datamodels`)
=========================================

:Release: |release|
:Date: |today|


LICENCE
~~~~~~~

Copyright
^^^^^^^^^
Copyright (c) 2010-2017 The JWST MIRI European Consortium software team. 
All rights reserved.

The miri data models have been developed by the MIRI EC software team 
as part of the JWST/MIRI consortium, which includes the following 
organisations: Ames Research Center, USA; Netherlands Foundation for 
Research in Astronomy; CEA Service d'Astrophysique, Saclay, France; 
Centre Spatial de Liege, Belgium; Consejo Superior de Investigacones 
Cientificas, Spain; Danish Space Research Institute; Dublin Institute 
for Advanced Studies, Ireland; EADS Astrium, Ltd., European Space 
Agency, Netherlands; UK; Institute d'Astrophysique Spatiale, France; 
Instituto Nacional de Tecnica Aerospacial, Spain; Institute of 
Astronomy, Zurich, Switzerland; Jet Propulsion Laboratory, USA; 
Laboratoire d'Astrophysique de Marseille (LAM), France; Lockheed 
Advanced Technology Center, USA; Max-Planck-Insitut fur Astronomie 
(MPIA), Heidelberg, Germany; Observatoire de Paris, France; Observatory 
of Geneva, Switzerland; Paul Scherrer Institut, Switzerland; 
Physikalishes Institut, Bern, Switzerland; Raytheon Vision Systems, USA; 
Rutherford Appleton Laboratory (RAL), UK; Space Telescope Science 
Institute, USA; Toegepast-Natuurwetenschappelijk Ondeszoek (TNOTPD), 
Netherlands; UK Astronomy Technology Centre (UKATC); University College, 
London, UK; University of Amsterdam, Netherlands; University of Arizona, 
USA; University of Cardiff, UK; University of Cologne, Germany; 
University of Groningen, Netherlands; University of Leicester, UK; 
University of Leiden, Netherlands; University of Leuven, Belgium; 
University of Stockholm, Sweden, Utah State University USA.

Terms and Conditions of Use
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This software may be used and copied free of charge only for 
non-commercial research purposes. All copies of this software must 
contain this copyright statement and disclaimer. The MIRI consortium 
must be acknowledged in any publications arising from use of this 
software. If you make modifications to this software, you must clearly 
mark the software as having been changed and you must also retain this 
copyright and disclaimer.

Where this software uses facilities developed by other members of the 
MIRI consortium (e.g. its use of the JWST data models) it is also bound
by the  licences issued with those facilities (see the LICENCE files
released with the jwst software).

Disclaimer
^^^^^^^^^^
This software is available "as is", without warranty of any kind, either 
expressed or implied, including the implied warranties of 
merchantability and fitness for a specific purpose. By using this 
software you are assuming all risks and costs. In no event is the MIRI 
EC software team or the MIRI consortium liable for any damages or losses 
that might result from the use of this software.

Introduction
~~~~~~~~~~~~
This document describes the implementation details for the MIRI data
models. These data models are built on the underlying JWST data models
package, which is described here:

http://ssb.stsci.edu/doc/jwst/jwst/datamodels/index.html

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
General purpose miri.datamodels data may be found in the miri/datamodels/data
directory.

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

