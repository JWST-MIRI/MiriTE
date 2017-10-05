#!/bin/bash -f
#
# Run the unit tests for all data models, and save the results
# to log files.
#
python test_badpixel_model.py >>test_badpixel_model.log 2>&1
python test_dark_reference_model.py >>test_dark_reference_model.log 2>&1
python test_distortion_models.py >>test_distortion_models.log 2>&1
python test_dqflags.py >>test_dqflags.log 2>&1
python test_droop_model.py >>test_droop_model.log 2>&1
python test_exposure_model.py >>test_exposure_model.log 2>&1
python test_flatfield_model.py >>test_flatfield_model.log 2>&1
python test_fluxconversion_models.py >>test_fluxconversion_models.log 2>&1
python test_fringe_frequences_model.py >>test_fringe_frequencies_model.log 2>&1
python test_gain_model.py >>test_gain_model.log 2>&1
python test_illumination_model.py >>test_illumination_model.log 2>&1
python test_ipc_model.py >>test_ipc_model.log 2>&1
python test_jump_model.py >>test_jump_model.log 2>&1
python test_lastframe_model.py >>test_lastframe_model.log 2>&1
python test_latent_model.py >>test_latent_model.log 2>&1
python test_linearity_model.py >>test_linearity_model.log 2>&1
python test_measured_model.py >>test_measured_model.log 2>&1
python test_pce_model.py >>test_pce_model.log 2>&1
python test_photometric_models.py >>test_photometric_models.log 2>&1
python test_pixel_saturation_model.py >>test_pixel_saturation_model.log 2>&1
python test_psf_models.py >>test_psf_models.log 2>&1
python test_readnoise_model.py >>test_readnoise_model.log 2>&1
python test_reset_model.py >>test_reset_model.log 2>&1
python test_reset_switch_charge_decay_model.py >>test_reset_switch_charge_decay_model.log 2>&1
python test_spectral_spatial_resolution_model.py >>test_spectral_spatial_resolution_model.log 2>&1
python test_straylight_model.py >>test_straylight_model.log 2>&1
python test_telescope_emission_model.py >>test_telescope_emission_model.log 2>&1
python test_transmission_correction_model.py >>test_transmission_correction_model.log 2>&1
python test_wavelength_correction_model.py >>test_wavelength_correction_model.log 2>&1
python test_miri_filters.py >>test_miri_filters.log 2>&1
python test_miri_measurement.py >>test_miri_measurement.log 2>&1
