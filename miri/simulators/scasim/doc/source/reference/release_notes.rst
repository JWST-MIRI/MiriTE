MIRI SCAsim release notes (:mod:`miri.simulators.scasim`)
=========================================================

:Release: |release|
:Date: |today|

The MIRI SCASim package simulates the behaviour of the
MIRI detectors. The package makes use of the MIRI
tools and MIRI data models, and therefore has the same
dependencies as the miri.tools and miri.datamodels
packages.

The dependencies should be taken care of automatically
by the MIRICLE installation script.
See http://miri.ster.kuleuven.be/bin/view/Public/MirisimInstallation

Full MIRISim-compatible release. The release supports the following features:

    * API via a class method instead of a global function.

    * 3 different classes can be used at the same time to keep track
      of 3 detectors.

    * Simulator parameters are obtained from scasim command line arguments
      or from the metadata contained within the illumination input files.
      Detailed detector properties are defined within configuration files
      or taken from calibration data products.

    * Parsing of intensity and wavelength arrays within an input FITS file
      of detector illumination data.
      
    * Applying a fringe pattern or flat-field pattern to the illumination
      data.
 
    * Simulation of Poisson (shot) noise.
         
    * Simulation of bad pixels, as described in a MIRI calibration 
      data product.
         
    * Simulation of the effect of quantum efficiency as a function of
      wavelength.
      
    * Simulation of dark current, which can vary with detector
      temperature, location on the detector surface, group and integration,
      and can include hot pixels; as described in a MIRI calibration 
      data product.
         
    * Simulation of read noise as a function of temperature, as described
      in a MIRI calibration data product.
                    
    * Simulation of amplifier gain, as described in a MIRI calibration
      data product.
    
    * Simulation of non-linearity using sensitivity coefficients.
           
    * Simulation of different detector readout modes, including FASTINTAVG
      and FASTGRPAVG. Can also simulate readout modes used by other JWST
      instrument detectors.
           
    * Simulation of detector subarray modes. Input data and output
      data can be in a different subarray mode.
      
    * Simulation of detector zero-point drift and latency effects using
      coefficients derived from JPL detector testing (although these
      effects only manifest themselves when a sequence of exposures is
      simulated).

    * Any of the MIRI sensor chip assemblies (MIRIMAGE, MIRIFULONG or
      MIRIFUSHORT) may be simulated, together with their focal plane
      modules (FPM S/N 106, 104 or 105). More SCAs can be added by editing
      the detector configuration file.
           
    * Simulation of cosmic ray hits on the detector pixels (using an
      STScI library of cosmic ray events, with a fall-back of random
      generation of events if the library is not available).
           
    * Writing either a DHAS FITSWriter format file or a JWST level 1b
      FITS file in either 'cube' or 'hypercube' format (with realistic
      header keywords).

The following features could be improved:

    * The non-linearity coefficients are based on an approximate inverse
      of the MIRI calibration data product and could be improved.
      
    * The simulation of detector latency and drift effects could be
      improved, in particular the effect of flux building up on the
      detector in between integrations is not yet simulated, and the
      effect in SLOW mode has not been calibrated.

    * The bias level simulated is arbitrary. Turning simulation effects
      on or off can affect the bias level, and the fact that zero-point
      drift is not simulated in SLOW mode can give the impression of a
      difference in bias between SLOW and FAST mode. The absolute DN
      zero-point may appear to be different, but the relative behaviour
      is still the same. The slopes generated are not affected.

    * The "tree ring" and first frame effects are inherited from the
      dark current calibration data product. The last frame effect is
      not currently simulated.

    * Reset switch charge decay is not currently simulated.

For details see http://miri.ster.kuleuven.be/bin/view/Public/MIRISimPublicReleases
