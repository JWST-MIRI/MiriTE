MIRI SCAsim release notes (:mod:`miri.simulators.scasim`)
=========================================================

:Release: |release|
:Date: |today|

Full MIRISim-compatible release. The release supports the following features:

    * API via a class method instead of a global function.

    * 3 different classes can be used at the same time to keep track
      of 3 detectors.

    * Simulator parameters are obtained from scasim command line arguments,
      from the parameters defined within configuration files and from
      calibration data products.

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
    
    * Simulation of non-linearity using polynomial coefficients.
           
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
      a configuration file.
           
    * Simulation of cosmic ray hits on the detector pixels (using an
      STScI library of cosmic ray events, with a fall-back of random
      generation of events if the library is not available).
           
    * Writing either a FITSWriter format file or a level 1b FITS file in
      either 'cube' or 'hypercube' format (with realistic header keywords).

The following features could be improved:
      
    * The non-linearity coefficients are assumed to be constant over the
      detector surface, rather than varying from pixel to pixel. These
      coefficients are based on an approximate inverse of the MIRI
      calibration data product and could be improved.
      
    * The simulation of detector latency and drift effects could be
      improved, in particular the effect of flux building up on the
      detector in between integrations is not yet simulated, and the
      effect in SLOW mode has not been calibrated.
      