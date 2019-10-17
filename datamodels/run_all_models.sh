#!/bin/bash
#
# Run the in-line tests for all data models, and save the results
# to log files.
#
python miri_aperture_correction_model.py >>miri_aperture_correction_model.log 2>&1
python miri_badpixel_model.py >>miri_badpixel_model.log 2>&1
python miri_dark_reference_model.py >>miri_dark_reference_model.log 2>&1
python miri_distortion_models.py >>miri_distortion_models.log 2>&1
python miri_droop_model.py >>miri_droop_model.log 2>&1
python miri_exposure_model.py >>miri_exposure_model.log 2>&1
python miri_flatfield_model.py >>miri_flatfield_model.log 2>&1
python miri_fluxconversion_models.py >>miri_fluxconversion_models.log 2>&1
python miri_fringe_frequencies_model.py >>miri_fringe_frequencies_model.log 2>&1
python miri_gain_model.py >>miri_gain_model.log 2>&1
python miri_illumination_model.py >>miri_illumination_model.log 2>&1
python miri_ipc_model.py >>miri_ipc_model.log 2>&1
python miri_jump_model.py >>miri_jump_model.log 2>&1
python miri_lastframe_model.py >>miri_lastframe_model.log 2>&1
python miri_latent_model.py >>miri_latent_model.log 2>&1
python miri_linearity_model.py >>miri_linearity_model.log 2>&1
python miri_measured_model.py >>miri_measured_model.log 2>&1
python miri_model_base.py >>miri_model_base.log 2>&1
python miri_pce_model.py >>miri_pce_model.log 2>&1
python miri_photometric_models.py >>miri_photometric_models.log 2>&1
python miri_pixel_saturation_model.py >>miri_pixel_saturation_model.log 2>&1
python miri_psf_models.py >>miri_psf_models.log 2>&1
python miri_readnoise_model.py >>miri_readnoise_model.log 2>&1
python miri_reset_model.py >>miri_reset_model.log 2>&1
python miri_reset_switch_charge_decay_model.py >>miri_reset_switch_charge_decay_model.log 2>&1
python miri_spectral_spatial_resolution_model.py >>miri_spectral_spatial_resolution_model.log 2>&1
python miri_straylight_model.py >>miri_straylight_model.log 2>&1
python miri_telescope_emission_model.py >>miri_telescope_emission_model.log 2>&1
python miri_transmission_correction_model.py >>miri_transmission_correction_model.log 2>&1
python miri_wavelength_correction_model.py >>miri_wavelength_correction_model.log 2>&1
python miri_filters.py >>miri_filters.log 2>&1
python miri_measurement.py >>miri_measurement.log 2>&1
python spectraldata.py >>spectraldata.log 2>&1
