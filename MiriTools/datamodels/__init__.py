#!/usr/bin/env python

"""

datamodels
==========

Subpackage datamodels contains the MIRI data product 
definitions and related utilities.
See the LICENCE file for terms and conditions.

Available modules
-----------------
dqflags:
    DQ flags values and operations.

operations:
    Common functions for data operations.

plotting:
    Common functions for data plotting.
    
util:
    Low level global constants and utility functions.

miri_model_base:
    Generic classes describing the top level MIRI data model.

miri_measured_model:
    Generic classes for MIRI data with error and quality information.

cdplib:
    Functions for searching for and accessing MIRI Calibration
    Data products. 

cdp:
    Contains all the MIRI Calibration Data products, as listed below.

miri_aperture_correction_model:
    MIRI aperture correction data.

miri_badpixel_model:
    MIRI bad pixel data.

miri_dark_reference_model:
    MIRI dark reference data.

miri_distortion_model:
    MIRI distortion models.
    
miri_droop_model:
    MIRI droop calibration information.

miri_flatfield_model:
    MIRI pixel, fringe or sky flat-field data.

miri_fluxconversion_model:
    MIRI flux conversion model.
    
miri_fringe_frequencies_model:
    MIRI fringe frequencies data.
    
miri_gain_model:
    MIRI amplifier gain information
    
miri_ipc_model:
    MIRI ipc deconvolution model.
  
miri_jump_model:
    MIRI jump detection information.
    
miri_lastframe_model:
    MIRI last frame correction model.
    
miri_latent_model:
    MIRI latent calibration information.

miri_linearity_model:
    MIRI linearity correction coefficients.

miri_pce_model:
    MIRI Photometric Conversion Efficiency model (used by the ETC).
    
miri_photometric_models:
    MIRI photometric conversion models (used by the pipeline).

miri_pixel_saturation_model:
    MIRI pixel saturation data

miri_psf_models:
    JWST/MIRI point-spread-function data.

miri_telescope_emission_model:
    JWST/MIRI telescope emission data.

miri_transmission_correction_model:
    MIRI tranmission correction data.

miri_wavelength_calibration_model:
    MIRI wavelength calibration coefficients
    
sim:
    Contains all the simulator support data products, as listed below.
    
miri_exposure_model:
    MIRI exposure data model.
    
miri_illumination_model:
    MIRI detector illumination model
    
filters:
    Filter tools and detector quantum efficiency measurements.
    
measured_variable
    Tools for managing variable measurements.

Scripts
-------
cdp_convert.py:
    Convert a MIRI Calibration Data product from one format to another.
    
cdp_print.py:
    Display the contents of a MIRI Calibration Data product.
    
cdp_verify.py
    Verify the integrity of a MIRI Calibration Data product.
    
multicdp_convert.py
    Convert all the MIRI Calibration Data products found within a
    file path.

multicdp_verify.py
    Verify the integrity of all the MIRI Calibration Data products
    found within a file path.
    
make_filters_fits.py:
    Create FITS file containing filter data.

run_make_filters_fits.sh:
    Call make_filters_fits to generate FITS file in data/filters.

Data
----
data/filters/*.txt,*.fits:
    MIRI filters data from fm_filters.xls on the FMTestData wiki page.

data/examples/filter_example.py:
    Demonstrate usage of the filters module

:History:
13 Nov 2012: Created
14 Jan 2013: Split into dataproduct and olddataproduct
04 Jun 2013: Shortened the names of the ramp, slope and image models.
02 Oct 2013: Added dqflags and updated documentation.
11 Dec 2013: Import the open function.
02 Oct 2014: Documentation update
19 Aug 2015: Removed MiriImageModel and MiriCubeModel (not used for CDPs).
18 Aug 2016: miri.miritools.measurements merged with miri.miritools.dataproduct.
06 Sep 2016: Renamed to miri.datamodels
29 Jun 2017: Updated to use build 7.1 data models.
05 Jan 2018: More version control information added. SVN info dropped.

"""

__project__ = 'MIRI Data Model Software'
__author__ = 'MIRI Software Team'
__maintainer__ = 'Steven Beard: steven.beard@stfc.ac.uk'
__copyright__ = '2018, %s' % __author__
__version__ = '7.0'

# Import common data products
from miri.datamodels.miri_model_base import MiriDataModel
from miri.datamodels.miri_measured_model import MiriMeasuredModel, \
    MiriRampModel, MiriSlopeModel, MiriSimpleModel

# Import the calibration data products
import miri.datamodels.cdp

# Import the simulator support data products
import miri.datamodels.sim

# Import the data model open function.
from miri.datamodels.util import open

