Overview of MIRI data product scripts in datamodels/scripts/


The following scripts are valid for the current data models and are
installed by default.
===========================================================================

CDP verification and testing
----------------------------
cdp_print.py		- Prints the contents of a CDP file.
cdp_verify.py		- Verify a CDP is compatible with requirements.
multicdp_verify.py	- ^ Run the above script on several files.
cdp_tests.py		- Test the use of the get_cdp() function.

CDP conversion scripts
----------------------
cdp_add_history.py	- Adds HISTORY records to CDP data.
cdp_add_subarray.py	- Adds missing subarray keywords to CDP data.
multicdp_subarray.csh	- ^ Run the above script on several files.
cdp_correct_band.py	- Corrects the BAND keyword within MRS CDP data.
multicdp_band.csh	- ^ Run the above script on several files.
cdp_reduce_dark.py	- Makes a smaller version of DARK CDP data.
convert_fits_to_asdf.py	- Read data model from FITS and save to ASDF.
convert_slope_data.py	- Convert DHAS-format slope data to level 2.

CDP fetching utilities
----------------------
cdp_get_doc.py		- Finds document associated with CDP file.
find_me_another.py	- Given a CRDS file, find up to date MIRI CDP.

General examples
----------------
dqflags_examples.py	- Some examples showing the use of DQ flags.


The following  scripts depend on specific versions of CDPs or data models.
They are not currently installed or relevant but have been saved to a
"historical" folder from which they can be reused in the future.
===========================================================================

Possibly reusable conversion scripts
------------------------------------
cdp_convert.py		- Convert CDPs from one format to another
					  (currently only CDP-2 to CDP-3).
multicdp_convert.py	- ^ Run the above script on several files.
cdp_remove_junk.py	- Remove left over FITS extensions.
convert_bad_pixel_mask.py - Convert DHAS-format bad pixel mask to CDP.
convert_dark_reference.py - Convert DHAS-format DARK file to CDP.

Obsolete scripts
----------------
convert_droop.py	- Generate DROOP CDP from ASCII table.
convert_flat_field.py	- Convert old-format field-field to CDP.
convert_im_colcorr.py	- Generate COLCORR CDP from ASCII table.
convert_im_distort.py	- Convert old-format distortion map to CDP.
convert_im_fluxconv.py	- Generate FLUXCONV CDP from ASCII table.
convert_im_flux_to_phot.py - Convert MiriImagingFluxconversionModel to
                             MiriImagingPhotometricModel.
convert_linearity.py	- Convert old-format linearity file to CDP.
convert_mrs_d2c.py	- Convert old-format D2C file to CDP.
convert_mrs_straylight.py - Convert old-format straylight file to CDP.
convert_psf.py		- Convert old-format PSF file to CDP.
convert_rsrf.py		- Generate RSRF CDP from ASCII table.
convert_srf.py		- Generate SRF CDP from ASCII table.

Scripts used in the past to generate test data
----------------------------------------------
make_filters_fits.py	- Generate MiriFilter data from ASCII table.
run_make_filters_fits.sh - An example running the above script.
make_measurements_fits.py - Generate MiriMeasurement from ASCII table.
run_make_measurements_fits.sh - An example running the above script.
make_qe_fits.py		- Generate MiriQuantumEfficiency data from ASCII table.
run_make_qe_fits.sh	- An example running the above script.
