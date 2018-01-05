#!/usr/bin/env python

"""

Module poisson_integrator - Contains the PoissonIntegrator and
ImperfectIntegrator classes.

PoissonIntegrator is a general purpose class representing an integrator,
such as a detector, which accumulates counts as a function of time. The
current count can be sampled non-destructively at any time, but the
reading obtained is subject to Poisson noise, depending on how many counts
have accumulated since the last reading. The class simulates reset,
integrate and readout functions.

ImperfectIntegrator is based on PoissonIntegrator but also simulates
imperfections, such as persistence, latency and zero-point drifts.

These general purpose classes are not specific to any simulator, although
ImperfectIntegrator implements the detector latency and zero-point drift
effects simulated by SCASim.

:History:
18 Jun 2010: Created.
28 Jun 2010: Reworked to recommended Python syntax standards.
05 Jul 2010: Removed superfluous print statement.
14 Jul 2010: Used the more efficient scipy.stats.poisson.rvs
             function instead of np.random.poisson and removed
             the limit parameter, which is no longer needed.
             Renamed from PoissionIntegrator to PoissonIntegrator
             (typo).
21 Jul 2010: Added plotting methods.
02 Aug 2010: Removed superflous blank lines and underlines from
             print output.
06 Aug 2010: Documentation formatting problems corrected.
03 Sep 2010: Added leak function.
06 Sep 2010: Added bad pixel mask.
14 Sep 2010: Added flag to turn Poisson noise on and off.
30 Sep 2010: Added a persistence factor as a class variable.
11 Oct 2010: Corrected bug in check for non-zero expected count.
15 Nov 2010: Bad pixels count removed from FITS header. Now added
             by caller.
16 Nov 2010: Added function to set the random number generator seed,
             for testing.
22 Nov 2010: Pass in the good pixel value rather than hard-coding it.
26 Jan 2010: get_header function replaced by get_poisson_noise_flag
             function.
04 Mar 2011: Documentation tweaks to resolve bad formatting.
22 Mar 2011: Plotting functions modified to use matplotlib figures
             and axes more flexibly (unfortunately the colorbar
             function doesn't work properly).
25 Mar 2011: Persistence function expanded so that a linear factor or
             polynomial can both be used.
29 Mar 2011: Check for very high counts that can no longer be stored
             as an integer.
01 Apr 2011: Documentation formatting problems corrected.
14 Jul 2011: Use of exceptions made more consistent and documented.
23 Sep 2011: Development and testing code now includes an example
             of how to demonstrate persistence in the integrator.
             Extreme levels of flux will now only cause an
             exception in the integrate() function if a saturation
             level has not been defined.
             Ability to choose whether to plot the data as a matrix
             or as an image.
29 Sep 2011: Added the ability to preflash the integerator when created,
             to allow persistence effects to be checked.
05 Oct 2011: Persistence coefficients in testing section imported from
             detector_properties.py.
20 Oct 2011: Destructor added to remove large objects.
15 Nov 2011: _MAXINT defined from sys.maxint. _MAXEXPECTED reduced to
             avoid an integer overflow from scipy.stats.poisson.
11 Jan 2012: Renamed symbol to avoid potential name clash: sum-->esum
28 Mar 2012: BUG DISCOVERED! matplotlib import crashes Python 2.7 when
             combined with the detector_properties import. Reason unknown.
             Worked around by commenting out the offending code.
02 Apr 2012: Corrected segmentation fault caused by numpy =/ operator.
             Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
21 May 2012: Plotting modified to use the miriplot utilities.
15 Nov 2012: Moved to top level simulator package.
05 Jun 2013: Allow the bad pixel mask to be cleared as well as defined.
23 Apr 2014: Changed the terminology to match common MIRI usage (e.g.
             an "integration" includes an interval between resets, not
             a fragment of that interval). Removed redundant methods.
             ImperfectIntegrator class added to separate effects
             due to other than Poisson statistics.
07 May 2014: Experimental version with a first implementation of the
             zeropoint drift and slope latency effects.
08 May 2014: Additional drift effects added using the "trapped charge"
             effect suggested by Dan and Rene.
05 Jun 2014: Simulations matching the JPL 3 test results implemented.
             Two latency effects are now implemented - a slow effect
             using an algorithm similar to Dan and Rene's, and a fast
             effect using an algorithm similar to Alistair's.
             STILL UNDER DEVELOPMENT.
06 Jun 2014: Change decay parameter from a number to timescale in seconds.
17 Jun 2014: Improved memory management. Changed verbosity definitions.
             Removed the old "preflash" latency implementation.
19 Jun 2014: Simplified the zeropoint drift function and parameterised
             the drift with time. Remodelled the persistence function so
             it operates on the signal above the zero level and therefore
             doesn't interfere with the zeropoint drift.
             Changed bad pixel mask from uint8 to uint16.
25 Jun 2014: set_linearity method added. Tests run in main section
             adjusted and documented.
01 Jul 2014: Corrected amplifier gain.
03 Sep 2015: Corrected bug in application of bad pixel mask.
08 Sep 2015: Made compatible with Python 3.
08 Oct 2015: Comparison between unequal types corrected in __str__ method.
10 Nov 2015: Added more explicit data type conversions and removed +=
             and *= operators to work around a numpy type-casting problem
             on a Mac.
04 Dec 2015: Added a logger. Renamed from poisson_integrator to integrators.
12 Feb 2016: Expected count rounded to avoid problems with extreme floating
             point values. Zero or negative values are replaced by 1.
10 Mar 2016: Bad pixel simulation moved to Detector object. set_pedestal
             function added, so the integrator can be programmed to count
             from a pre-defined zero level.
11 Mar 2016: Return of the numpy type-casting problem on a Mac. Explicitly
             typed the pedestal array to be the same type as expected_count
             and removed the += operator.
23 Mar 2016: Filter out NaNs from expected count array before they trigger
             a problem in the poisson.rvs function.
05 May 2016: Linearity renamed sensitivity to avoid confusion with the
             coefficients in the linearity CDP.
15 Jul 2016: Additional debugging.
02 Nov 2016: Test with parameters obtained from detector properties by
             default.
14 Feb 2017: Corrected bug in Poisson noise calculation. The zeropoint level
             was being included in the calculation and incorrectly boosting
             the noise. The zeropoint is now kept separate until readout.
21 Apr 2017: Restored the plot_polynomial function and added plot_corrected.
27 Jul 2017: Added an option to the wait() method to simulate a non-zero
             background flux.
26 Oct 2017: Take into account the correlation between samples when calculating
             Poisson noise.
01 Nov 2017: Limit the extrapolation of zero point drift by applying a maximum
             clock time.
14 Dec 2017: Adjust Poisson noise calculation to be based on last readout,
             not on last expected value.
04 Jan 2018: Corrected mistake which caused the zeropoint to be added to
             every frame during Poisson noise calculation. Last readout
             should be with respect to the zeropoint. Some variables
             renamed to better names.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default logger

from copy import deepcopy
import math
import sys
import scipy.stats
import numpy as np

# The maximum allowable signed integer value. This is the largest particle
# count which can be read from a Poisson integrator. The bucket size
# will be set to this value if not specified.
_MAXINT = sys.maxint

# The maximum level size of the expected count before the scipy.stats.poisson
# function overflows. It must be significantly less than _MAXINT to prevent
# an overflowing value being randomly generated.
# TODO: This is an arbitrary level set by trial and error. Can it be more exact?
_MAXEXPECTED = _MAXINT - (100L * math.sqrt(_MAXINT))

# The maximum clock time for which the slow zeropoint drift is valid.
# The detector is assumed to settle after this time has elapsed.
_MAXCLOCK = 100000.0

def linear_regression(x, y):
    """
    
    Linear regression function. Fit an intercept and a slope
    to a set of x,y coordinates
    
    """
    matrix = np.vstack( [x, np.ones_like(x)] ).T
    slope, intercept = np.linalg.lstsq(matrix,y)[0]
    return (slope, intercept)


class PoissonIntegrator(object):
    """
    
    Class PoissonIntegrator - Simulates the behaviour of a particle
    counting detector in which the actual number of particles (be they
    photons or electrons, etc...) detected varies with the expected
    number according to Poisson statistics.
    
    This class simulates a perfect integrator.
                 
    :Parameters:
    
    rows: int
        The number of detector rows. Must be at least 1.
    columns: int
        The number of detector columns. Must be at least 1.
    particle: string, optional, default="photon"
        The name of the particles being counted.
    time_unit: string, optional, default="seconds"
        The unit of time measurement.
    bucket_size: int, optional, default=None
        The maximum number of particles that may be counted.
        If None there is no limit.
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for
        example to observe what effects in a simulation are caused
        by Poisson noise).
    verbose: int, optional, default=2
        Verbosity level. Activates print statements when non-zero.
        
        * 0 - no output at all except for error messages
        * 1 - warnings and errors only
        * 2 - normal output
        * 3, 4, 5 - additional commentary
        * 6, 7, 8 - extra debugging information

    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
        
    :Raises:
    
    ValueError
        Raised if any of the initialisation parameters are out of range.
        
    :Attributes:
    
    shape: tuple of 2 ints
        The shape of the integrator array.
    particle: string
        The name of the particle/charge being counted (e.g. "photon" or
        "electron").
    time_unit: string
        The name of the time unit in which integration durations are
        measured.
    bucket_size: int
        The maximum number of paricles/charges that can be accumulated
        within any element. The element saturates when the bucket size
        is reached.
    expected_count: float array
        The particles/charges accumulated within the elements of the
        integrator at the current time.
    last_count: float array
        The expected count of the integrator when last sampled.
    last_readout: int array
        The reading obtained from the integrator when the expected_count
        was last sampled.
    flux: float array
        The last flux on which the integrator was integrated.
    nints: int
        The number of integrations executed since the PoissonIntegrator
        object was created, or since a new exposure was started.
        An integration is defined as the period between resets.
        An integration consists of one or more integration periods.
        The count may be sampled many times during an integration, and
        a graph of count against time typically looks like a ramp whose
        slope depends on the input flux.
    nperiods: int
        The number of integration periods executed since the beginning
        of the integration. An integration period is a period of
        time during which the counter is integrated on a flux. Normally,
        the counter is read out at the end of the integration flagment.
    readings: int
        The number of readings made since the beginning of the integration.
        Normally, this will be the same as the number of integration
        periods, but there is nothing to stop the counter being integrated
        more than once in between readings, or being read out more than once
        after an integration period. The latter strategy can be used to
        make a group of readings to reduce read noise.
    nperiods_at_readout: int
        The number of integration periods executed before the last
        readout.
    exposure: float
        The total exposure time for the current integration, i.e. the
        total accumulated integration time since the last reset.
    
    """

    def __init__(self, rows, columns, particle="photon", time_unit="seconds",
                 bucket_size=None, simulate_poisson_noise=True, verbose=1,
                 logger=LOGGER):
        """
        
        Constructor for class PoissonIntegrator.
        
        Parameters: See class doc string.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("integrators")
        self.verbose = int(verbose)
        if self.verbose > 3:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug("+++" + str(self.__class__.__name__) + \
                " object created with rows=" + str(rows) + \
                " columns=" + str(columns) + \
                " bucket_size=" + str(bucket_size))

        # The Poisson integrator must have a sensible size.
        if int(rows) <= 0 or int(columns) <= 0:
            strg = "The Poisson integrator must have a non-zero size."
            raise ValueError(strg)
        
        # Define the fundamental properties of the integrator: its shape,
        # units and bucket depth.
        self.shape = (int(rows), int(columns))
        self.particle = particle
        self.time_unit = time_unit
        self.bucket_size = bucket_size
        self.pedestal = None
                
        # All the internal counters are initialised to zero.
        # NOTE: The expected count is maintained in a floating point
        # array which is truncated to integer when read out.
        self.zeropoint = np.zeros(self.shape)
        self.expected_count = np.zeros(self.shape)
        self.last_count = np.zeros(self.shape)
        self.last_readout = np.zeros(self.shape, dtype=np.uint32)
        self.flux = np.zeros(self.shape)
        self.nperiods = 0    # Number of integration periods since reset.
        self.nints = 0       # Number of integrations since creation or new exposure.
        self.readings = 0    # Number of readings/groups since reset.
        self.nperiods_at_readout = 0
        self.exposure_time = 0.0  # Exposure time since reset.
        self.clock_time = 0.0     # Clock time since switch on.
        self.time_at_reset = 0.0  # Clock time at last reset
        
        # This flag controls whether Poisson noise is simulated
        self.simulate_poisson_noise = simulate_poisson_noise
        if not simulate_poisson_noise and self.verbose > 0:
            self.logger.info("NOTE: Poisson noise is turned off.")
    
    def __del__(self):
        """
        
        Destructor for class PoissonIntegrator.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        try:
            # Objects created in the constructor
            del self.flux
            del self.last_readout
            del self.last_count
            del self.expected_count
            del self.zeropoint
        except Exception:
            pass
    
    def get_title(self, description=''):
        """
        
        Get a string giving a title for the measurements.
        
        :Parameters:
        
        description: string, optional
            A description to be added to the title, if required. 
        add_comment: bool, optional
            Set to True if the comment string (if present) should
            be added to the title. The default is False.
        
        :Returns:
        
        title: string
            A title string
        
        """
        strg = "Expected %s count" % self.particle
        if description:
            strg += " %s" % description
        return strg
    
    def get_counts(self):
        """
        
        Get a copy of the array containing the number of particles counted.
        This can be used to preserve a record of a previous exposure
        before an integrator object is deleted.
        
        :Returns:
        
        count: nparray (int)
            The number of particles counted so far.

        """
        # Use a deepcopy to ensure the count is preserved even
        # when this integrator object is destroyed.
        return deepcopy(self.expected_count)

    def set_seed(self, seedvalue=None):
        """
        
        Set the seed for the numpy random number generator.
        This function can be used while testing to ensure the
        Poisson noise generated subsequently is well defined.
        
        :Parameters:
        
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator.
            If not specified, a value of None will be sent, which
            randomises the seed.
            
        """
        np.random.seed(seedvalue)
        
    def _apply_zeropoint(self, data, nresets=1):
        """
        
        Apply the zero point function to a set of data.
        The perfect integrator has a zero point of exactly zero.
        
        """
        # Data are zeroed perfectly.
        return data * 0.0

    def _apply_integrator(self, data, flux, time):
        """
        
        Apply an integrator function to the data.
        This is a perfectly linear integrator.
        
        """
        if self.verbose > 6:
            self.logger.debug("+++Perfect integrator")
        newdata = data + (flux * float(time)) # Linear sensitivity
                
        # Check that the count is not going to come anywhere near
        # overflowing (which causes problems later for the Poisson
        # noise calculation).
        if newdata.max() > _MAXEXPECTED:
            toobig = np.where(newdata > _MAXEXPECTED)
            newdata[toobig] = _MAXEXPECTED
        return newdata
 
    def set_pedestal(self, pedestal):
        """
        
        Set a pedestal level which acts as the ultimate zero point
        for the integrator. All resets and integrations start from
        this level. By default, the pedestal is zero.
        
        :Parameters:
        
        pedestal: array-like
            If provided, a pedestal level which is added to the count
            after the reset. Must be the same size and shape as the
            The array must have the same size and shape as the integrator. 
        
        """
        if pedestal is not None:
            pedestal = np.asarray(pedestal, dtype=self.expected_count.dtype)
            if pedestal.shape[0] == self.shape[0] and \
               pedestal.shape[1] == self.shape[1]:
                self.pedestal = pedestal
            else:
                strg = "Pedestal array has the wrong size (%d x %d). " % \
                    pedestal.shape
                strg += "It should be %d x %d." % self.shape
                raise TypeError(strg)
        else:
            self.pedestal = None

    def reset(self, nresets=1, new_exposure=False):
        """
        
        Reset the integrator. All internal counters are reset to the
        pedestal level (normally zero), although persistence and
        zeropoint drift might leave a residual signal.
        
        :Parameters:
        
        nresets: int (optional)
            The number of resets to be applied. (This only makes
            a difference when there are persistence effects.)
            By default there will be one reset.
        new_exposure: bool (optional)
            Set True if this reset begins a new exposure (and the
            integration count needs to be reset).
            By default this is False.
        
        """
        if self.verbose > 5:
            self.logger.debug("+++Reset " + str(nresets) + \
                              " times. New exposure=" + str(new_exposure))
        # Update the integration number
        if new_exposure:
            self.nints = 1      # Reset the integration count
        else:
            self.nints += 1     # Increment the integration count
        # Reset integration counters
        self.nperiods = 0
        self.readings = 0
        self.nperiods_at_readout = 0

        # Apply the zeropoint function, depending on the number of
        # resets requested.
        signal = self.zeropoint + self.expected_count
        self.zeropoint = self._apply_zeropoint(signal, nresets=nresets)
        if self.pedestal is not None:
            self.zeropoint = self.zeropoint + self.pedestal
        self.expected_count = self.expected_count * 0.0
        self.last_count = self.last_count * 0.0

        # Reset exposure time and last readout memory, and set the
        # clock time measured at reset.
        self.exposure_time = 0.0
        self.last_readout = self.last_readout * 0
        self.time_at_reset = self.clock_time

    def integrate(self, flux, time):
        """
        
        Integrate on a particle flux for a certain length of time.
        Calling this method adds one integration period.
        
        :Parameters:
        
        flux: array_like
            Photon flux array in particles per time unit. flux must
            either be None (indicating an integration on nothing) or be
            the same shape and size as the photon counter. The array must
            contain positive values - any negative values are treated as
            if they are zero.
        time: float
            The integration time in time units. It must not be negative.
            
        :Raises:
    
        ValueError
            Raised if the photon flux is too large.
        TypeError
            Raised if any of the parameters are of the wrong type,
            size or shape.
            
        """
        if self.verbose > 5:
            self.logger.debug("+++Integrate for " + str(time) + " " + str(self.time_unit))
        # The integration time must always be positive.
        if time < 0.0:
            raise ValueError("Integration time must be positive")
        
        # The clock time is always incremented.
        self.clock_time += float(time)
        
        # A flux array set to None implies an integration on nothing.
        # If a flux array is provided, it must be the same size and shape as
        # the internal counter array or be capable of being broadcast by
        # numpy onto this array. Anything else will raise an exception.
        #
        if flux is not None:
            # The integration period count and exposure time are increased
            # whenever there is an integration on flux.
            self.nperiods += 1
            self.exposure_time += float(time)
            
            # Make sure that the given flux is a numpy array (otherwise the
            # flux.shape attribute will not be defined and the shape test will
            # fail).
            self.flux = np.asarray(flux)
            if self.flux.shape == self.shape:      
                # The input flux array must always be positive.
                # Replace negative values with zero.
                arenegative = np.where(flux < 0.0)
                self.flux[arenegative] = 0.0
       
                # Apply the integrator function
                self.expected_count = \
                    self._apply_integrator(self.expected_count, self.flux, time)

                if self.verbose > 3:
                    strg = "Expected count after integration: "
                    strg += "min=%.2f to max=%.2f." % \
                        (self.expected_count.min(), self.expected_count.max())
                    self.logger.debug( strg )
            else:
                # Faulty flux array given - raise an exception.
                strg = "Input flux array has the wrong shape: ("
                for ii in self.flux.shape:
                    strg += "%d " % ii
                strg += ") instead of ("
                for ii in self.shape:
                    strg += "%d " % ii
                strg += ")"
                raise TypeError(strg)

    def wait(self, time, bgflux=0.0):
        """
        
        Pause for a certain length of time without integrating.
        
        :Parameters:
        
        time: float
            The wait time in time units. It must not be negative.
        bgflux: float (optional)
            A uniform background photon flux falling on the integrator
            while it is idle. This flux can affect detector persistence
            by filling charge traps. By default there is no background
            flux.
            
        :Raises:
    
        ValueError
            Raised if the time is negative.
        TypeError
            Raised if any of the parameters are of the wrong type,
            size or shape.
            
        """
        if self.verbose > 5:
            self.logger.debug("+++Wait for " + str(time) + " " + str(self.time_unit))
        # The integration time must always be positive.
        if time < 0.0:
            raise ValueError("Wait time must be positive")

        # The clock time is always incremented.
        self.clock_time += float(time)

        if bgflux > 0.001:
            # Integrate the counter on the background flux.
            # TODO: Improve this simulation.
            self.expected_count = \
                self._apply_integrator(self.expected_count, bgflux, time)
        
    def readout(self, nsamples=1):
        """
        
        Read the signal stored within the integrator and apply Poisson
        noise (if switched on).
        
        :Parameters:
        
        nsamples: int, optional, default=1
            The number of samples from the Poisson distribution.
                
        :Returns:
        
        readout_array: array_like
            The latest count as read out (with Poisson noise).
            
        """
        if self.verbose > 5:
            self.logger.debug("+++Readout " + str(self.readings+1) + \
                              " with " + str(nsamples) + " samples.")
        if int(nsamples) <= 0:
            strg = "Number of samples must be at least 1."
            raise ValueError(strg)
        
        # Simulate Poisson noise when the poisson_noise flag is True,
        if self.simulate_poisson_noise:
            # When the Poisson integrator is read out, the expected number of
            # particles is turned into an actual number of particles by making
            # one or more samples from the Poisson distribution. The Poisson
            # noise will be reduced in cases where more than one sample is made.
            # Note that scipy.stats.poisson.rvs will report an exception if
            # any of the expected counts are zero or negative, so these are
            # filtered out. The function also fails with a domain error if it
            # encounters very small negative floating point values. These are
            # removed by rounding the expected count.
            if self.verbose > 5:
                strg = "Readout %d: Expected count is %g to %g" % \
                    (self.readings+1, self.expected_count.min(), \
                     self.expected_count.max())
                self.logger.debug(strg)
            # Round extreme values.
            rounded_diff = np.around(self.expected_count-self.last_count, 4)
            # Filter out zero, negative or bad values.
            where_zero = np.where(rounded_diff <= 0)
            if where_zero:
                rounded_diff[where_zero] = 1
            where_nan = np.where(np.isnan(rounded_diff))
            if where_nan:
                rounded_diff[where_nan] = 1
            try:
                # Take a random sample from the Poisson distribution.
                read_diff = scipy.stats.poisson.rvs(rounded_diff)
                # If required, take more samples and average.
                if nsamples > 1:
                    for ii in range(1, nsamples):
                        read_diff = read_diff + \
                            scipy.stats.poisson.rvs(rounded_diff)
                    read_diff = read_diff / float(nsamples)
            except ValueError, e:
                strg = "Poisson rvs error. Difference array:\n"
                strg += str(rounded_diff)
                strg += "\n" + str(e)
                raise ValueError(strg)
            # Put the zero values back to zero.
            if where_zero:
                read_diff[where_zero] = 0
            readout_array = self.zeropoint + self.last_readout + read_diff
        else:
            readout_array = self.zeropoint + self.expected_count
        
        # The latest readout cannot be less than the previous readout or
        # (if specified) larger than the defined bucket size.
        if self.bucket_size is None:
            maxread = _MAXINT
        else:
            maxread = self.bucket_size
        readout_array = np.clip(readout_array, 0.0, maxread)
        
        # Remember the last expected count
        self.last_count = deepcopy(self.expected_count)
        
        # Increment the reading/group counter and save the last readout.
        # NOTE: The last readout is measured from the zeropoint.
        self.readings += 1
        self.last_readout = readout_array.astype(np.uint32) - self.zeropoint
        self.nperiods_at_readout = self.nperiods
        
        if self.verbose > 5:
            min = np.min( readout_array )
            max = np.max( readout_array )
            mean = np.mean( readout_array )
            self.logger.debug("   min=%.1f, max=%.1f, mean=%.1f" % (min, max, mean))
        
        # Return the readout converted to int to make a discrete count.
        return readout_array.astype(np.uint32)
 
    def __str__(self):
        """
        
        Return a string describing a PoissonIntegator object.
        
        """
        strg = "%s: %s counter of %d rows x %d columns" % \
            (self.__class__.__name__, self.particle.title(),
             self.shape[0], self.shape[1])
        if self.bucket_size is not None:
            strg += " with a bucket size of %d" % self.bucket_size
        strg += '.'
        strg += "\n    Clock time since creation/switch-on: %.1f %s" % \
            (self.clock_time, self.time_unit)
        if self.time_unit == 'seconds' or self.time_unit == 's':
            strg += " (%.2f mins)." % (self.clock_time / 60.0)
        else:
            strg += "."
        
        strg += "\n    There have been "
        if self.nperiods <= 0:
            strg += "no integration periods since the last reset."
        else:
            strg += "%d integration periods totalling %.1f %s" % \
                (self.nperiods, self.exposure_time, self.time_unit)
            elapsed_time = self.clock_time - self.time_at_reset
            if abs(elapsed_time - self.exposure_time) > 0.001:
                strg += " (%.1f %s elapsed)" % \
                    (elapsed_time, self.time_unit)
            strg += " since the last reset."
            strg += "\n    The zero point ranges from %d to %d %ss." % \
                (self.zeropoint.min(), self.zeropoint.max(), \
                 self.particle)
            strg += "\n    The expected count ranges from %d to %d %ss." % \
                (self.expected_count.min(), self.expected_count.max(), \
                 self.particle)
            if self.nperiods_at_readout <= 0:
                strg += "\n    There have been no readouts."
            else:
                strg += "\n    Reading %d (above zeropoint) after integration period %d " \
                    "obtained %d to %d %ss." % \
                    (self.readings, self.nperiods_at_readout, \
                     self.last_readout.min(), self.last_readout.max(), \
                     self.particle)
        return strg
    
    def plot(self, plotfig=None, plotaxis=None, datatype='image',
             withbar=True, title=None, description=''):
        """
        
        Plot the counter array. The subplot and show parameters
        allow more than one plot to be displayed within one plotting
        area (pyplot axis).
        
        :Parameters:
          
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        datatype: string, optional
            The type of data being displayed: 'matrix' or 'image'.
            'matrix' is better for displaying an array (with square
            pixel boundaries) and 'image' is better for displaying
            a smoothly interpolated image. The default is 'image'.
        withbar: bool, optional
            Set to False to suppress the colour bar. The default is True.
        title: string, optional
            Optional title for the plot. The default is a string
            describing the PoissonIntegrator object.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
        description: string, optional
            Additional description to be shown on the plot, if required.

        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
        # Create a new figure to contain the plot (larger than default size).
        if title is None:
            title = self.get_title(description=description)
        elif description:
            title += "\n%s" % description

        signal = self.zeropoint + self.expected_count
        mplt.plot_image(signal, plotfig=plotfig, plotaxis=plotaxis,
                        datatype=datatype, withbar=withbar,
                        xlabel=None, ylabel=None, title=title)


class ImperfectIntegrator(PoissonIntegrator):
    """
    
    Class ImperfectIntegrator - Simulates the behaviour of a particle
    counting detector in which the actual number of particles (be they
    photons or electrons, etc...) detected varies with the expected
    number according to Poisson statistics.
    
    The integrator is not perfect. Bad pixels, persistence, drift, latency and
    non-linear effects are simulated.
                 
    :Parameters:
    
    rows: int
        The number of detector rows. Must be at least 1.
    columns: int
        The number of detector columns. Must be at least 1.
    particle: string, optional, default="photon"
        The name of the particles being counted.
    time_unit: string, optional, default="seconds"
        The unit of time measurement.
    bucket_size: int, optional, default=None
        The maximum number of particles that may be counted.
        If None there is no limit.
    simulate_poisson_noise: boolean, optional, default=True
        A flag that may be used to switch off Poisson noise (for
        example to observe what effects in a simulation are caused
        by Poisson noise).
    verbose: int, optional, default=2
        Verbosity level. Activates print statements when non-zero.
        
        * 0 - no output at all except for error messages
        * 1 - warnings and errors only
        * 2 - normal output
        * 3, 4, 5 - additional commentary
        * 6, 7, 8 - extra debugging information

    logger: Logger object (optional)
        A Python logger to handle the I/O. This parameter can be used
        by a caller to direct the output to a different logger, if
        the default defined by this module is not suitable.
        
    :Raises:
    
    ValueError
        Raised if any of the initialisation parameters are out of range.

    :Attributes:
    
    persistence: float or tuple of floats
        A description of the persistence properties of the
        integrator.
        
        * If a single float is provided, this is assumed to be a
          persistence factor in the range 0.0 to 1.0, where
          0.0 means no persistence - signal is zeroed after reset.
          1.0 means 100% persistence - reset has no effect on the signal.
        
        * If a tuple of floats is provided this is taken as a
          set of polynomial coefficients describing a non-linear
          persistence effect (e.g. if the persistence is worse at
          high signals there may be a second order factor).

    sensitivity: tuple of 2 floats
        The constant and slope of the function indicating how the
        sensitivity of the integrator varies at high counts
        compared to the well depth.
        [1.0, 0.0] means the integrator is perfectly linear
        [0.0, 0.0] means the integrator is completely dead
        [1.0, -1.0] means the counter starts at 100% sensitivity and
        reduces to 0% at full well.
        [1.0, -0.5] means the counter starts at 100% sensitivity and
        reduces to 50% at full well.
        [1.5, -1.0] means the counter starts at 100% sensitivity, remains
        level at 100% sensitivity until 50% full well, then reduces to
        50% sensitivity at full well.
        The integrated sensitivity function gives the linearity function.
    zp_slow: tuple of 2 floats or None:
        A list of coefficients describing the slow change in
        zeropoint as a function of clock time. For now, only
        a linear drift is simulated. The coefficients are
        [starting_dn, driftfactor], where starting_dn is a
        constant shift in DN and driftfactor is the amount of
        drift with time in DN/s.
        Set to None to turn off the long term drift altogether.
    zp_fast: tuple of (tuple of floats) or None
        A 2-D array describing the short term zeropoint jumps
        as a function of integration number and flux.
        A tuple of floats is provided for each integration number,
        and each tuple of floats contains coefficients describing
        the zeropoint shift seen as a function of incident flux.
        For now, only a linear model is used, and each set of
        coefficients is described by 2 floats, [const, jumpfactor],
        where const is a ZP shift in DN which is always applied
        for that integration number and jumpfactor in s describes
        how the drift changes with source flux.
        Set to None to turn off zeropoint jumps altogether.
    slow_latency_params: tuple of 2 floats or None
        The gain and decay factor of the slow latency function.
        Set to None to turn off the slow latency effect altogether.
    fast_latency_parameters: tuple of 2 floats or None
        The gain and decay factor of the fast latency function.
        Set to None to turn off the fast latency effect altogether.
         
    """

    def __init__(self, rows, columns, particle="photon", time_unit="seconds",
                 bucket_size=None, simulate_poisson_noise=True, verbose=2,
                 logger=LOGGER):
        """
        
        Constructor for class ImperfectIntegrator.
        
        Parameters: See class doc string.
        
        """
        # Create the parent object.
        super(ImperfectIntegrator, self).__init__(rows, columns,
                                                  particle=particle,
                                                  time_unit=time_unit,
                                                  bucket_size=bucket_size,
                                                  simulate_poisson_noise=simulate_poisson_noise,
                                                  verbose=verbose, logger=logger)
        
        # Define quantities that make the integrator imperfect.
        # Integrator persistence in terms of fraction of count
        # remaining after a reset.
        self.persistence = 0.0
        self.persistence_poly = None # No persistence polynomial function.

        # Changes in the sensitivity as a function of fraction of full well.
        # which leads to a non-linear response.
        self.sensitivity = None # No sensitivity change - constant response.
        
        # Zero point drift as function of time and zero point
        # jumps as a function of source brightness and
        # integration number.
        self.zp_slow = None    # No slow zeropoint drift
        self.zp_fast = None    # No fast zeropoint jumps
        
        # Detector latency effects.
        # Gain and decay factor for slow and fast latent effects.
        # Begin with these effects turned off.
        self.slow_latency_params = None
        self.fast_latency_params = None

        self.last_flux = self.flux
        self.slow_latent = np.zeros(self.shape, dtype=np.uint32)
# UNCOMMENT FOR EXTRA PLOTTING
#         self.m_factor = np.zeros(self.shape, dtype=np.uint32)
#         self.s_factor = np.zeros(self.shape, dtype=np.uint32)
        
    def __del__(self):
        """
        
        Destructor for class ImperfectIntegrator.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        try:
            # Objects created in the constructor
            if self.last_flux is not None:
                del self.last_flux
            if self.slow_latent is not None:
                del self.slow_latent

            # Objects created in methods
            if self.persistence_poly is not None:
                del self.persistence_poly
             
            # Finally, tidy up the parent class
            super(ImperfectIntegrator, self).__del__()
        except Exception:
            pass
        
    def set_persistence(self, persistence):
        """
        
        Define the integrator persistence effects.
        Persistence effects are observed if a reset() does not
        completely clear the count from a previous integration.

        :Parameters:
        
        persistence: float or tuple of floats or None
            A description of the persistence properties of the
            integrator.
        
            * If a single float is provided, this is assumed to be a
              persistence factor in the range 0.0 to 1.0, where
              0.0 means no persistence - signal is zero after reset.
              1.0 means 100% persistence - reset has no effect on the signal.
        
            * If a tuple of floats is provided this is taken as a
              set of polynomial coefficients describing a non-linear
              persistence effect (e.g. if the persistence is worse at
              high signals there may be a second order factor).
              
            * Set to None to turn off persistence effects altogether
              (equivalent to setting it to 0.0).
            
        :Raises:
    
        TypeError
            Raised if the persistence parameter is invalid.
        
        """
        if self.verbose > 3:
            self.logger.debug("+++Set persistence to " + str(persistence))
        if persistence is None:
            # Turn off persistence
            self.persistence = 0.0
            self.persistence_poly = None
        elif isinstance(persistence,(float,int)):
            # A single value is provided. Assume this is a linear
            # persistence factor, which must be 0.0 to 1.0.
            if float(persistence) >= 0.0 and float(persistence) <= 1.0:
                self.persistence = float(persistence)
                self.persistence_poly = None # No polynomial function.
            else:
                strg = "Persistence factor must be between 0.0 and 1.0"
                raise ValueError(strg)
        elif isinstance(persistence,(tuple,list)):
            # The persistence factor is assumed to be polynomial.
            # All factors must be in the range 0.0 to 1.0.
            self.persistence = persistence
            # The polynomial function needs to be updated
            if self.persistence_poly is not None:
                del self.persistence_poly # Avoid memory leak.
            self.persistence_poly = np.poly1d(self.persistence)
        else:
            raise TypeError("Invalid persistence factor")

    def set_sensitivity(self, sensitivity):
        """
        
        Define the integrator sensitivity effects, which cause nonlinearity.

        :Parameters:
        
        sensitivity: tuple of 2 floats or None
            The constant and slope of the function indicating how the
            sensitivity of the integrator varies at high counts
            compared to the well depth.
            None means the integrator is perfectly linear
            [1.0, 0.0] means the integrator is perfectly linear
            [0.0, 0.0] means the integrator is completely dead
            [1.0, -1.0] means the counter starts at 100% sensitivity and
            reduces to 0% at full well.
            [1.0, -0.5] means the counter starts at 100% sensitivity and
            reduces to 50% at full well.
            [1.5, -1.0] means the counter starts at 100% sensitivity, remains
            level at 100% sensitivity until 50% full well, then reduces to
            50% sensitivity at full well.
            The integrated sensitivity function gives the linearity function.
        
        """
        if sensitivity is not None:
            self.sensitivity = sensitivity
        else:
            self.sensitivity = [1.0, 0.0]

    def set_zeropoint(self, zp_slow, zp_fast):
        """
        
        Define the integrator zeropoint effects. Zeropoint effects
        are observed if the counter does not start counting from
        zero after a reset. The zeropoint can vary with time,
        source brightness and integration number.

        :Parameters:
        
        zp_slow: tuple of 2 floats or None:
            A list of coefficients describing the slow change in
            zeropoint as a function of clock time. For now, only
            a linear drift is simulated. The coefficients are
            [starting_dn, driftfactor], where starting_dn is a
            constant shift in DN and driftfactor is the amount of
            drift with time in DN/s.
            Set to None to turn off the long term drift altogether.
        zp_fast: tuple of (tuple of floats) or None
            A 2-D array describing the short term zeropoint jumps
            as a function of integration number and flux.
            A tuple of floats is provided for each integration number,
            and each tuple of floats contains coefficients describing
            the zeropoint shift seen as a function of incident flux.
            For now, only a linear model is used, and each set of
            coefficients is described by 2 floats, [const, jumpfactor],
            where const is a ZP shift in DN which is always applied
            for that integration number and jumpfactor in s describes
            how the drift changes with source flux.
            Set to None to turn off zeropoint jumps altogether.
            
        :Raises:
    
        TypeError or AssertionError
            Raised if the zeropoint coefficients parameter is invalid.
        
        """
        if self.verbose > 3:
            self.logger.debug("+++Set zeropoint parameters to" + \
                " slow=" + str(zp_slow) + " fast=" + str(zp_fast))
        # Exit immediately if there is a problem.
        if zp_slow is not None:
            assert isinstance(zp_slow,(tuple,list))
            assert len(zp_slow) > 1
        if zp_fast is not None:
            assert isinstance(zp_fast,(tuple,list))
            assert len(zp_fast) > 1
            assert isinstance(zp_fast[0],(tuple,list))
            assert len(zp_fast[0]) > 1
        self.zp_slow = zp_slow
        self.zp_fast = zp_fast

    def set_latency(self, slow_parameters, fast_parameters):
        """
        
        Defines integrator latency effects. Latency effects are caused
        by the integrator retaining a memory of previous integrations
        and taking a finite time to react to changes. There are several
        latency effects covering different timescales:
        
        * Slow latency effects are caused by a build up of sensitivity
          during exposures, leading to each successive integration
          having a slightly greater slope than the previous integration.
          The magnitude of the effect depends on the source brightness
          and time for which each source is observed.
          The gain factor is given in (1/counts).
          The decay factor is given in seconds.
          
        * Fast latency effects are caused by rapid changes in flux.
          Instead of the flux changing instantly from one level to
          another, the integrator behaves as if the flux level is
          decaying exponentially between one level and the next.
          The gain factor is dimensionless.
          The decay factor is given in seconds.
          
        Both latency effects are described by two parameters:
        
        * The gain factor describes the magnitude of the effect.
        
        * The decay factor describes how quickly the effect decays.
        
        :Parameters:
        
        slow_parameters: tuple of 2 floats or None
            The gain and decay factor of the slow latency function.
            Set to None to turn off the slow latency effect altogether.
        fast_parameters: tuple of 2 floats or None
            The gain and decay factor of the fast latency function.
            Set to None to turn off the fast latency effect altogether.

        :Raises:
    
        TypeError or AssertionError
            Raised if either of the latency parameters are invalid.
          
        """
        if self.verbose > 3:
            self.logger.debug("+++Set latency parameters to" + \
                " slow=" + str(slow_parameters) + \
                " fast=" + str(fast_parameters))
        # Exit immediately if there is a problem.
        if slow_parameters is not None and fast_parameters is not None:
            assert isinstance(slow_parameters,(tuple,list))
            assert len(slow_parameters) > 1
        if fast_parameters is not None:
            assert isinstance(fast_parameters,(tuple,list))
            assert len(fast_parameters) > 1            
            
        self.slow_latency_params = slow_parameters
        self.fast_latency_params = fast_parameters

    def _persistence_function(self, data, coeffs, zerolevel=0.0):
        """
        
        Apply the persistence function to some data. The persistence
        function only applies to the signal above the defined zero level.
        
        """
        if coeffs is None:
            # No persistence
            newdata = data * 0.0
        elif isinstance(coeffs,(float,int)):
            # The persistence type is linear
            newdata = zerolevel + (data - zerolevel) * coeffs
        else:
            # Polynomial persistence. Use the polynomial function
            # prepared earlier or create one now.
            if self.persistence_poly is None:
                self.persistence_poly = np.poly1d(coeffs)
            newdata = zerolevel + self.persistence_poly(data - zerolevel)
        # The persistence cannot go below zero or above the bucket size.
        if self.bucket_size is None:
            maxread = _MAXINT
        else:
            maxread = self.bucket_size
        return np.clip(newdata, 0, maxread)
        
    def _zeropoint_function(self, shape, slow_coeffs, fast_coeffs, integ, flux,
                            clock_time):
        """
        
        Calculate a zeropoint function from the integration number,
        incoming flux and time since switch on.
        
        NOTE: integ starts at 1.
        
        """
        assert isinstance(shape,(list,tuple))
        if slow_coeffs is not None:
            assert isinstance(slow_coeffs,(list,tuple))
            assert len(slow_coeffs) > 1
        if fast_coeffs is not None:
            assert isinstance(fast_coeffs,(list,tuple))
            assert len(fast_coeffs) > 1
            assert isinstance(fast_coeffs[0],(list,tuple))
            assert len(fast_coeffs[0]) > 1

        # Start with a base zeropoint
        zeropoint = np.zeros(shape)
        
        # The zeropoint drifts slowly as a function of clock time in seconds,
        # up to a maximum clock time.
        if slow_coeffs is not None:
            zeropoint += slow_coeffs[0] + \
                         (slow_coeffs[1] * min(clock_time, _MAXCLOCK))
        
        # Integrations 2,3,4,5 onwards have an additional drift
        # dependent on the incoming flux.
        if (fast_coeffs is not None) and integ > 1:
            ii = min(integ - 2, (len(fast_coeffs)-1))
            zeropoint += fast_coeffs[ii][0] + (fast_coeffs[ii][1] * flux)
# Extra line left over from when the coeffs represented a 2nd order polynomial
#                             + coeffs[ii][2] * flux * flux
                    
        # The zeropoint cannot go below zero or above the bucket size.
        if self.bucket_size is None:
            maxread = _MAXINT
        else:
            maxread = self.bucket_size
        return np.clip(zeropoint, 0, maxread)
    
    def _sensitivity_function(self, slope, const, data):
        """
        
        Calculate the sensitivity of the integrator, as a function
        of current count. This function represents the ability of
        the integrator to receive new counts as a function of the
        headroom left between the current count and the bucket size.

        The change in sensitivity as a function of count headroom
        will make the integrator response non-linear.
        
        """
        if self.bucket_size is not None:
            # The sensitivity depends on the relative headroom between
            # the current count and the bucket size. Assume a linear
            # factor for now.
            sensitivity = const + (slope * data / self.bucket_size)
        else:
            # If there is no bucket size, the sensitivity is always 1.0
            sensitivity = np.ones_like(data)
        # Sensitivity must be in the range 0.0 (insensitive) to 1.0 (perfect)
        # NOTE: Actually, MIRI linearity measurements show the sensitivity
        # actually starts out around 1/0.86, so the limit could be 1.2.
        # return np.clip(sensitivity, 0.0, 1.2)
        return np.clip(sensitivity, 0.0, 1.0)

    def _latency_function(self, flux, last_flux):
        """
         
        Apply the latency functions to the incoming flux.
         
        """
        # The integrator behaves as if it has a memory of
        # the previous flux.
        # First apply the fast latency function, where a "fast_gain"
        # fraction of the difference between the previous flux and the
        # current flux is retained.
        if self.fast_latency_params is not None:
            fast_beta = self.fast_latency_params[0]
            newflux = flux + fast_beta * (last_flux - flux)
        else:
            newflux = flux
# UNCOMMENT FOR EXTRA PLOTTING
#             if np.min(flux) > 0.0:
#                 self.s_factor = 1.0 + (self.fast_latency_params[0] * (last_flux - flux) / flux)
#             else:
#                 self.s_factor = np.ones(self.shape)

        # Now apply the slow latency function, where the gain is
        # boosted by a multiplication factor dependent on the amount
        # of slow latency build-up.
        if self.slow_latency_params is not None:
            slow_beta = self.slow_latency_params[0]
            m_factor = 1.0 + slow_beta * self.slow_latent
            newflux = newflux * m_factor
#             if self.verbose > 7:
#                 self.logger.debug("+++Boosting flux by " + str(m_factor))
# UNCOMMENT FOR EXTRA PLOTTING
#             self.m_factor = m_factor
        return newflux

    def _apply_zeropoint(self, data, nresets=1):
        """
        
        Apply the zero point function to a set of data
        
        """
        # An imperfect integrator builds up a latent image based
        # on the count remaining before the reset and the slow
        # decay factor.
        if self.slow_latency_params is not None:
            slow_gamma = 1.0 - (self.exposure_time / self.slow_latency_params[1])
            self.slow_latent =  slow_gamma * self.slow_latent + data
# No evidence for an increase in slow latent effect with number of resets.
#             # Apply additional reset effects if requested.
#             # Additional resets ramp up the slow latent without decaying it,
#             # but each subsequent reset has only half the effect.
#             if nresets > 1:
#                 effect = 0.5
#                 for n in range(0,nresets-1):
#                     self.slow_latent =  self.slow_latent + data * effect
#                     effect = effect / 2.0
   
        # Remember the flux from the previous integration, smoothing
        # any sharp changes by applying the fast latent decay factor.
        if self.fast_latency_params is not None:
            fast_gamma = 1.0 - (self.exposure_time / self.fast_latency_params[1])
            self.last_flux = self.flux + \
                (self.last_flux - self.flux) * fast_gamma
        else:
            self.last_flux = self.flux

        # Calculate the current zeropoint drift, based on the clock time,
        # input flux and integration number.
        zeropoint = self._zeropoint_function(data.shape, self.zp_slow,
                                             self.zp_fast, self.nints,
                                             self.flux, self.clock_time)
        
        # An imperfect integrator will retain a fraction of the previous
        # signal (above the zeropoint), as defined by a persistence function.
        newdata = self._persistence_function(data, self.persistence,
                                             zerolevel=0.0)
#                                           zerolevel=zeropoint)
        # Apply additional resets if requested.
        if nresets > 1:
            for n in range(0,nresets-1):
                newdata = self._persistence_function(newdata, self.persistence,
                                                     zerolevel=0.0)
#                                                   zerolevel=zeropoint)
            
        newdata = zeropoint + newdata
        return newdata.astype(np.int)

    def _apply_integrator(self, data, flux, time):
        """
        
        Apply an integrator function to the data.
        This is a non-linear integrator whose sensitivity varies
        with the headroom between the existing count and the
        bucket size.
        
        """
        if self.verbose > 6:
            self.logger.debug("+++Imperfect integrator")
        # The incoming flux is affected by the integrator latency
        # function, in which a certain amount of the flux from the
        # previous integration is retained.
        newflux = self._latency_function(flux, self.last_flux)
        
        # The new count depends on the integrated flux and the sensitivity
        # of the integrator.
        if self.sensitivity is not None:
            newcounts = newflux * float(time) * \
                self._sensitivity_function(self.sensitivity[1], self.sensitivity[0],
                                           data)
        else:
            newcounts = newflux * float(time)
            
        if self.verbose > 3:
            strg = "Adding flux counts ranging from "
            strg += "min=%.2f to max=%.2f." % (newcounts.min(), newcounts.max())
            self.logger.debug( strg )
                    
        temp_data = data + newcounts
                
        # Check that the count is not going to come anywhere near
        # overflowing (which causes problems later for the Poisson
        # noise calculation). Negative counts are also not allowed.
        newdata = np.clip(temp_data, 0, _MAXEXPECTED)
        del temp_data, newcounts
        return newdata

    def anneal(self):
        """
        
        Anneal the (imperfect) integrator by leaking away the persistence.
        
        :Parameters:

        None yet
        
        """
        if self.verbose > 5:
            self.logger.debug("+++Anneal ")
           
        # Remove latents and persistence 
        self.zeropoint = np.zeros(self.shape)
        self.expected_count = np.zeros(self.shape)
        self.clock_time = 0.0     # Clock time since switch on.
            
        # Finish with a reset
        self.reset()

    def hit_by_cosmic_ray(self, energy_map, row, column):
        """
        
        Dump some cosmic ray energy onto a portion of the integrator.
        
        :Parameters:
        
        energy_map: array_like or number
            A single value or a 2-D array of cosmic ray energies in terms
            of equivalent particle count. The cosmic ray will affect the
            specified row and column plus neighbouring rows and columns,
            depending on the size of the energy map.
        row: int
            The row on which the cosmic ray energy is centred.
        column: int
            The column on which the cosmic ray energy is centred.
            
        :Raises:
    
        TypeError
            Raised if the energy map is invalid.
            
        """
        if self.verbose > 5:
            esum = 0.0
            if isinstance(energy_map,(int,float)):
                esum = energy_map
            else:
                energy_map = np.asarray(energy_map, dtype=self.expected_count.dtype)
                esum = energy_map.sum()
            self.logger.debug("Cosmic ray hit of %.0f %s around pixel (%d,%d)." % \
                (esum, self.particle, row, column))

        if isinstance(energy_map,(int,float)):
            # energy_map is a single value
            if row >= 0 and row < self.expected_count.shape[0] and \
               column >= 0 and column < self.expected_count.shape[1]:
                self.expected_count[row,column] = \
                    self.expected_count[row,column] + energy_map
            else:
                if self.verbose > 4:
                    self.logger.debug("Cosmic ray at row %d column %d ignored." % \
                        (row,column))
        else:
            # energy_map is an array. Calculate the slices affected by the
            # cosmic ray energy.
            energy_map = np.asarray(energy_map, dtype=self.expected_count.dtype)
            if energy_map.ndim != 2:
                raise TypeError("Cosmic ray energy map must be a 2-D array")
            nrows = energy_map.shape[0]
            ncolumns = energy_map.shape[1]
            rowstart = row - nrows//2
            rowend = rowstart + nrows
            colstart = column - ncolumns//2
            colend = colstart + ncolumns
        
            # Add the energy to the particle count. (If the counter
            # goes above the bucket size it will be truncated at the
            # next readout.)
            if rowstart >= 0 and rowend < self.expected_count.shape[0] and \
               colstart >= 0 and colend < self.expected_count.shape[1]:
            
                # The slicing is only simple if the energy array is entirely
                # within the particle count array.
                self.expected_count[rowstart:rowend, colstart:colend] = \
                    self.expected_count[rowstart:rowend, colstart:colend] + \
                    energy_map
            
            elif rowstart < self.expected_count.shape[0] and rowend >= 0 and \
                colstart < self.expected_count.shape[1] and colend >= 0:
            
                # The energy array is partly inside the particle count array,
                # so more complex assignment is needed. (This is a lot less
                # efficient than numpy slicing.)
                for rw in xrange(rowstart,rowend):
                    if rw >=0 and rw < self.expected_count.shape[0]:
                        for cl in xrange(colstart,colend):
                            if cl >=0 and cl < self.expected_count.shape[1]:
                                erw = rw - rowstart
                                ecl = cl - colstart
                                self.expected_count[rw,cl] = \
                                    self.expected_count[rw,cl] + \
                                    energy_map[erw,ecl]
            else:
                # Ignore a cosmic ray which is completely outside the array.
                if self.verbose > 4:
                    self.logger.debug("Cosmic ray at row %d column %d ignored." % \
                        (row,column))

    def leak(self, leakage, row, column):
        """
        
        Simulate the leakage of some charge from the integrator, as might
        happen if there was a brief short circuit or malfunction.
        
        :Parameters:
        
        leakage: int
            The number of particles leaked from the current count.
        row: int
            The row at which the leakage happens.
        column: int
            The column at which the leakage happens.
            
        """
        if self.verbose > 5:
            self.logger.debug("Leakage of %.0f %s around pixel (%d,%d)." % \
                (leakage, self.particle, row, column))
        
        # Only process leaks inside the detector boundary.
        if row >= 0 and row < self.expected_count.shape[0] and \
            column >= 0 and column < self.expected_count.shape[1]:
            # Subtract the leakage from the current count, ensuring it
            # does not go below zero.
            self.expected_count[row,column] = self.expected_count[row,column] - \
                leakage
            if self.expected_count[row,column] < 0:
                self.expected_count[row,column] = 0
            
            # Zero the last reading from this pixel, as it is no longer
            # relevant. This is one situation where a reading can be less
            # than the previous reading.
            self.last_readout[row,column] = 0
        else:
            if self.verbose > 4:
                self.logger.debug("Leakage at row %d column %d ignored." % \
                    (row,column))
 
    def __str__(self):
        """
        
        Return a string describing a PoissonIntegator object.
        
        """
        strg = super(ImperfectIntegrator, self).__str__()
        if self.persistence is not None:
            if isinstance(self.persistence, (tuple,list)):
                if len(self.persistence) > 0:
                    strg += "\n    A persistence of %s is applied to each reset." % \
                        str(self.persistence)
            elif self.persistence > 0.0:
                strg += "\n    A persistence of %f is applied to each reset." % \
                    self.persistence
        if self.zp_slow is not None:
            strg += "\n    Slow zeropoint drift: const=%f, factor=%f per %s" % \
            (self.zp_slow[0], self.zp_slow[1], self.time_unit)
        if self.zp_fast is not None:
            strg += "\n%s" % self._fast_zeropoint_str(spacer='    ')
        if self.slow_latency_params is not None:
            strg += "\n    %s" % self._slow_latency_str()
        if self.fast_latency_params is not None:
            strg += "\n    %s" % self._fast_latency_str()
#         if self.bad_pixels is not None:
#             bad_pixel_count = np.sum(self.bad_pixels != self.good_value)
#             strg += "\n    %d bad pixels are applied." % bad_pixel_count
            
        return strg
    
    def _fast_zeropoint_str(self, spacer=''):
        """
        
        Return a string describing the zeropoint coefficients
        
        """
        strg = "%sZeropoint jumps applied:" % spacer
        integ = 1
        for zpint in self.zp_fast:
            if integ % 2 > 0:
                strg += "\n%s  " % spacer
            else:
                strg += ";\t"
            strg += "Integration %d: %s" % (integ, str(zpint))
            integ += 1
        return strg

    def _slow_latency_str(self):
        """
        
        Return a string describing slow latency parameters
        
        """
        if self.slow_latency_params is not None:
            strg = "Slow latency gain=%g, timescale=%.1f %s" % \
                ((self.slow_latency_params[0], self.slow_latency_params[1],
                  self.time_unit))
        else:
            strg = 'No slow latency'
        return strg
    
    def _fast_latency_str(self):
        """
        
        Return a string describing fast latency parameters
        
        """
        if self.fast_latency_params is not None:
            strg = "Fast latency gain=%g, timescale=%.1f %s" % \
                ((self.fast_latency_params[0], self.fast_latency_params[1],
                  self.time_unit))
        else:
            strg = 'No fast latency'
        return strg

    def plot_persistence(self, plotfig=None, plotaxis=None, xmin=0,
                         xmax=None, npoints=512, description=''):
        """
        
        Plot integrator persistence as a function of counts present
        when reset to a self-contained matplotlib figure.
        
        :Parameters:
        
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xmin: int, optional, default=0
            The minimum number of counts to be plotted.
        xmax: int, optional, default=bucket size or 100000
            The maximum number of counts to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
            
        :Global variables:
        
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
            
        :Requires:
        
        matplotlib.pyplot
            
        :Raises:
     
        ImportError
            Raised if the matplotlib plotting library is not available.
            
        """
        if xmax is None:
            if self.bucket_size is not None:
                xmax = self.bucket_size
            else:
                xmax = 100000
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
        # First generate the persistence function
        counts = np.linspace(xmin, xmax, npoints)
        remaining = self._persistence_function(counts, self.persistence)
        if self.bucket_size is not None:
            toobig = np.where(remaining > self.bucket_size)
            remaining[toobig] = self.bucket_size
            
        # Create a figure to contain the plot (larger than default size).
        tstrg = "Persistence function (coeffs=%s)" % str(self.persistence)
        if description:
            tstrg += " %s" % description

        # Plot the actual gain using a blue line and gray grid.
        mplt.plot_xy(counts, remaining, plotfig=plotfig, plotaxis=plotaxis,
                     xlabel='Counts present at reset',
                     ylabel='Counts remaining after reset', title=tstrg)

    def plot_zeropoint(self, fmin=0, fmax=500, nints=5, npoints=10,
                       description=''):
        """
        
        Plot integrator zero point drift as a function of flux and
        integration number to a self-contained matplotlib figure.
        
        :Parameters:
        
        fmin: int, optional, default=0
            The minimum flux to be plotted.
        fmax: int, optional, default=500
            The maximum flux to be plotted.
        nints: int, default=5
            The number of integrations to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
            
        :Global variables:
        
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
            
        :Requires:
        
        matplotlib.pyplot
            
        :Raises:
     
        ImportError
            Raised if the matplotlib plotting library is not available.
            
        """
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
                    
        # Create a figure to contain the plot (larger than default size).
        tstrg = "Zeropoint function"
        if description:
            tstrg += " %s" % description
        plotfig = mplt.new_figure(1, figsize=(10,8), stitle=tstrg)
        plotaxis = None
        
        # Generate and plot a zeropoint function for each integration
        flux = np.linspace(fmin, fmax, npoints)
        lfmt = ['rx', 'yx', 'gx', 'bx', 'kx']
        for integ in range(1,nints+1):
            zeropoint = self._zeropoint_function(flux.shape, self.zp_slow,
                                                 self.zp_fast, integ, flux, 0.0)

            # Plot the actual gain using a blue line and gray grid.
            plotaxis = mplt.plot_xy(flux, zeropoint, plotfig=plotfig,
                                    plotaxis=plotaxis, linefmt=lfmt[integ-1],
                                    xlabel='Illumination in %ss per second' % self.particle,
                                    ylabel='Zeropoint offset in %ss' % self.particle,
                                    title='Integrations 1 to %d' % nints,
                                    showplot=False)
        mplt.show_plot()

    def plot_sensitivity(self, plotfig=None, plotaxis=None, xmin=0,
                         xmax=None, npoints=512, description=''):
        """
         
        Plot integrator sensitivity as a function of count
        to a self-contained matplotlib figure.
         
        :Parameters:
         
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xmin: int, optional, default=0
            The minimum number of counts to be plotted.
        xmax: int, optional, default=bucket size or 100000
            The maximum number of counts to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
             
        :Global variables:
         
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
             
        :Requires:
         
        matplotlib.pyplot
             
        :Raises:
      
        ImportError
            Raised if the matplotlib plotting library is not available.
             
        """
        if xmax is None:
            if self.bucket_size is not None:
                xmax = self.bucket_size
            else:
                xmax = 100000
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt

         # First generate the sensitivity function
        counts = np.linspace(xmin, xmax, npoints)
        if self.sensitivity is not None:
            tstrg = "Sensitivity function (slope=%f, const=%f)" % \
            (self.sensitivity[0], self.sensitivity[1])
            actual = self._sensitivity_function(self.sensitivity[1],
                                                self.sensitivity[0], counts)
        else:
            tstrg = "Sensitivity function (not defined)"
            actual = np.ones_like(counts)
             
        # Create a figure to contain the plot (larger than default size).
        if description:
            tstrg += " %s" % description
 
        # Plot the actual gain using a blue line and gray grid.
        mplt.plot_xy(counts, actual, plotfig=plotfig, plotaxis=plotaxis,
                     xlabel='Current count',
                     ylabel='Sensitivity', title=tstrg)

    def plot_linearity(self, plotfig=None, plotaxis=None, xmin=0,
                         xmax=None, npoints=512, description=''):
        """
         
        Plot integrator linearity as a function of count
        to a self-contained matplotlib figure.
         
        :Parameters:
         
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xmin: int, optional, default=0
            The minimum number of counts to be plotted.
        xmax: int, optional, default=bucket size or 100000
            The maximum number of counts to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
             
        :Global variables:
         
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
             
        :Requires:
         
        matplotlib.pyplot
             
        :Raises:
      
        ImportError
            Raised if the matplotlib plotting library is not available.
             
        """
        if xmax is None:
            if self.bucket_size is not None:
                xmax = self.bucket_size
            else:
                xmax = 100000
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
         # First generate the sensitivity function
        counts = np.linspace(xmin, xmax, npoints)
        if self.sensitivity is not None:
            tstrg = "Sensitivity function (slope=%f, const=%f)" % \
            (self.sensitivity[0], self.sensitivity[1])
            actual = self._sensitivity_function(self.sensitivity[1],
                                                self.sensitivity[0], counts)
        else:
            tstrg = "Sensitivity function (not defined)"
            actual = np.ones_like(counts)
        
        # The linearity function is the integral of the sensitivity
        # function, scaled by the number of points.
        linearity = np.zeros_like(actual)
        for ii in range(1,len(linearity)):
            linearity[ii] = linearity[ii-1] + actual[ii-1]
        linearity = npoints * linearity
             
        # Create a figure to contain the plot (larger than default size).
        if description:
            tstrg += " %s" % description
 
        # Plot the actual gain using a blue line and gray grid.
        mplt.plot_xy(counts, linearity, plotfig=plotfig, plotaxis=plotaxis,
                     xlabel='Input count',
                     ylabel='Expected output count', title=tstrg)

    def plot_polynomial(self, plotfig=None, plotaxis=None, xmin=0,
                         xmax=None, npoints=512, coeffs=[0.0, 1.0],
                         gain=5.5, description=''):
        """
          
        Plot a polynomial as a function of count
        to a self-contained matplotlib figure.
          
        :Parameters:
          
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xmin: int, optional, default=0
            The minimum number of counts to be plotted.
        xmax: int, optional, default=bucket size or 100000
            The maximum number of counts to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        coeffs: list of float
            Polynomial coefficients
        gain: float
            Gain in counts per DN
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
              
        :Global variables:
          
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
              
        :Requires:
          
        matplotlib.pyplot
              
        :Raises:
       
        ImportError
            Raised if the matplotlib plotting library is not available.
              
        """
        if xmax is None:
            if self.bucket_size is not None:
                xmax = self.bucket_size
            else:
                xmax = 100000
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
         # First generate the sensitivity function
        counts = np.linspace(xmin, xmax, npoints)
        counts = counts / gain
        tstrg = "Polynomial function (coeffs=%s)" % str(coeffs)
         
        # Evaluate the polynomial.
        linearity = np.polynomial.polynomial.polyval(counts, coeffs)
              
        # Create a figure to contain the plot (larger than default size).
        if description:
            tstrg += " %s" % description
  
        # Plot the actual gain using a blue line and gray grid.
        mplt.plot_xy(counts, linearity, plotfig=plotfig, plotaxis=plotaxis,
                     xlabel='Input DN',
                     ylabel='Expected output DN', title=tstrg)

    def plot_corrected(self, plotfig=None, plotaxis=None, xmin=0,
                         xmax=None, npoints=512, coeffs=[0.0, 1.0],
                         gain=5.5, description=''):
        """
         
        Generate a non-linear ramp from the sensitivity coefficients,
        correct it using the polynomial coefficients and plot the
        result to a self-contained matplotlib figure.
         
        :Parameters:
         
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
            If an axis is provided, the plot will be created within that
            axis. Providing an axis allows the caller to control exactly
            where a plot appears within a figure.
            If an axis is not provided, a new one will be created filling
            a new figure.
        xmin: int, optional, default=0
            The minimum number of counts to be plotted.
        xmax: int, optional, default=bucket size or 100000
            The maximum number of counts to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        coeffs: list of float
            Polynomial coefficients
        gain: float
            Gain in counts per DN
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
             
        :Global variables:
         
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
             
        :Requires:
         
        matplotlib.pyplot
             
        :Raises:
      
        ImportError
            Raised if the matplotlib plotting library is not available.
             
        """
        if xmax is None:
            if self.bucket_size is not None:
                xmax = self.bucket_size
            else:
                xmax = 100000
        # Import the miri.tools plotting module.
        import miri.tools.miriplot as mplt
         # First generate the sensitivity function
        counts = np.linspace(xmin, xmax, npoints)
        if self.sensitivity is not None:
            tstrg = "Sensitivity function (slope=%f, const=%f)" % \
            (self.sensitivity[0], self.sensitivity[1])
            actual = self._sensitivity_function(self.sensitivity[1],
                                                self.sensitivity[0], counts)
        else:
            tstrg = "Sensitivity function (not defined)"
            actual = np.ones_like(counts)
        counts = counts / gain
        
        # The linearity function is the integral of the sensitivity
        # function, scaled by the number of points.
        linearity = np.zeros_like(actual)
        for ii in range(1,len(linearity)):
            linearity[ii] = linearity[ii-1] + actual[ii-1]
        linearity = (npoints * linearity) / gain
        
        correction = np.zeros_like(linearity)
        for ii in range(0, npoints):
            linsquared = linearity[ii] * linearity[ii]
            correction[ii] = coeffs[0] + \
                             coeffs[1] * linearity[ii] + \
                             coeffs[2] * linsquared + \
                             coeffs[3] * linsquared * linearity[ii] + \
                             coeffs[4] * linsquared * linsquared

#         # Evaluate the polynomial.
#         correction = np.polynomial.polynomial.polyval(counts, coeffs)
#          
#         # The corrected ramp is the linearity simulation divided by
#         # the correction
#         corrected = linearity / correction
             
        # Create a figure to contain the plot (larger than default size).
        if description:
            tstrg += " %s" % description
 
        # Plot the actual gain using a blue line and gray grid.
        mplt.plot_xy(counts, correction, plotfig=plotfig, plotaxis=plotaxis,
                     xlabel='Input DN',
                     ylabel='Expected output DN', title=tstrg)

#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_poisson_integrator module.
# However, the tests here are made in verbose mode and include plotting,
# so they show what is happening in more detail.
#
if __name__ == '__main__':
    print("Testing the PoissonIntegrator and ImperfectIntegrator classes")
    # Import the miri.tools plotting module.
    import miri.tools.miriplot as mplt
    
    SIMULATE_NOISE = False

    # Adjust to determine what is plotted.
    PLOT_PERSISTENCE = False
    PLOT_ZEROPOINT   = False
    PLOT_LINEARITY   = True
    PLOT_LATENT      = False
    PLOT_RAMPS       = True
    PLOT_SLOPES      = True

    # Define the size of the integrator to be tested.
    # This should be small enough to run the tests quickly.
    ROWS = 6
    COLUMNS = 8

    # The amplifier gain is used to ensure the integrator works in electron
    # counts even though the flux levels for the JPL 3 tests are expressed
    # in DN/s. It may be defined explicitly or obtained from the amplifier
    # properties
    #
    #import miri.simulators.scasim.amplifier_properties as amplifier_properties
    #ampgain = amplifier_properties._amp493_1['GAIN'][1]
    ampgain = 0.1818      # DN/s = electrons/s * ampgain

    # The well depth determines the saturation level of the integrator
    # in electron counts. It may be defined here or obtained from the
    # detector properties.
    #
    #import miri.simulators.scasim.detector_properties as detector_properties
    #WELL_DEPTH = detector_properties._sca493['WELL_DEPTH']
    WELL_DEPTH = 250000

    # Create a perfect Poisson integrator object.
    perfect = PoissonIntegrator(ROWS, COLUMNS, verbose=3)
    print("===Status of perfect integrator after creation.")
    print( perfect )
    perfect.integrate(None, 10.0)
    
    # Check the integrator can handle strange data.
    zerodata = np.zeros( [ROWS,COLUMNS] )
    perfect.integrate(zerodata, 10.0)  
    rddata = perfect.readout(nsamples=4)
    negdata = zerodata - 1.0
    rddata = perfect.readout(nsamples=4)
    perfect.integrate(negdata, 10.0)
    # Integrate on random data 10 times.
    for rep in (0, 10):
        randdata = np.random.randn( ROWS,COLUMNS ) - 0.5
        perfect.integrate(randdata, 10.0)
        rddata = perfect.readout(nsamples=4)
    posdata = zerodata + 1.0
    perfect.integrate(posdata, 10.0)
    rddata = perfect.readout(nsamples=4)
    print("===Status of perfect integrator after some dodgy integrations.")
    print( perfect )
    del perfect
    
    # Create an imperfect Poisson integrator object.
    print("")
    integrator = ImperfectIntegrator(ROWS, COLUMNS, bucket_size=WELL_DEPTH,
                                     particle='electron', verbose=5)
    if not SIMULATE_NOISE:
        # Turn off the noise simulation, so we can concentrate on the
        # persistence and latency effects.
        integrator.simulate_poisson_noise = False
    
    # Define persistence (either explicitly or from the detector properties)
    #
    import miri.simulators.scasim.detector_properties as detector_properties
    persistence = detector_properties._sca493['PERSISTENCE']
#     persistence = [1.0e-8, 0.03, 0.0]
    integrator.set_persistence(persistence)
#     integrator.set_persistence(None)
    
    # Define the slow and fast zero point drift coefficients.
    #
    #import miri.simulators.scasim.detector_properties as detector_properties
    zpslow = detector_properties._sca493['ZP_SLOW']
    zpfast = detector_properties._sca493['ZP_FAST']
#     zpslow = [50000.0, 0.0084]
#     zpfast = [[0.0, -2.917],
#               [0.0, -2.292],
#               [0.0, -2.396],
#               [0.0, -2.408]]
    integrator.set_zeropoint(zpslow, zpfast)

    # Define a linearity/sensitivity factor
    #
    #import miri.simulators.scasim.detector_properties as detector_properties
    integrator.set_sensitivity(detector_properties._sca493['SENSITIVITY'])
#     integrator.set_sensitivity([1.3, -0.5])
    #integrator.set_sensitivity(None)
    
    # Define the slow and fast latency factors
    #
    #import miri.simulators.scasim.detector_properties as detector_properties
    (slow_gain, slow_decay) = detector_properties._sca493['LATENCY_SLOW']
    (fast_gain, fast_decay) = detector_properties._sca493['LATENCY_FAST']
#     slow_gain =  1.67e-9   # Scale factor for slow latent (electrons/s)
#     slow_decay = 136000.0  # Timescale of slow latent decay (s)
#     fast_gain = 0.002      # Scale factor for fast latent
#     fast_decay = 300.0     # Timescale of fast latent decay (s)
    integrator.set_latency([slow_gain,slow_decay], [fast_gain,fast_decay])
    
    latency_string = "%s; %s" % (integrator._slow_latency_str(),
                                 integrator._fast_latency_str())

    # Display the contents of the object. Plot if required.
    print("===Status of imperfect integrator after creation.")
    print( integrator )
    if PLOT_PERSISTENCE:
        integrator.plot_persistence(xmax=WELL_DEPTH,
            description="\nFraction of signal above zero level " + \
            "remaining after reset.")
    if PLOT_ZEROPOINT:
        integrator.plot_zeropoint(fmax=400/ampgain,
            description="(as a function of integration and flux)")
    if PLOT_LINEARITY:
        integrator.plot_sensitivity(description="\nChange in relative " + \
            "sensitivity as the count approaches full well")
        integrator.plot_linearity(description="\nIntegrated sensitivity")
#         integrator.plot_polynomial(coeffs=[0.0, 0.93435, 8.14895e-7, 2.616463e-11],
#                                    description="\nPolynomial")
        integrator.plot_polynomial(
            coeffs=[0.0, 0.86538, 4.64238e-6, -6.1809e-11, 7.2313e-16],
            description="\nPolynomial")
        integrator.plot_corrected(
            coeffs=[0.0, 0.86538, 4.64238e-6, -6.1809e-11, 7.2313e-16],
            description="\nIntegrated sensitivity / Correction")

    print("===Cosmic ray tests")
    energy_map = np.asarray(([1,2], [3,4]), dtype=np.int32)
    integrator.hit_by_cosmic_ray(energy_map, 3, 4)
    energy_map = np.asarray(([1,2], [3,4]), dtype=np.float64)
    integrator.hit_by_cosmic_ray(energy_map, 2, 3)
#     integrator.plot()

    # Test the effect of the integrator properties on the calculated slopes.
    print( "\n===Slope tests." )

    # These parameters can be matched to JPL test data.
    NEXPOSURES = 20
    NDARKEXP = 5
    NINTS = 5
    NGROUPS = 50
    FRAME_TIME = 2.785

    # Make the basic flux array all the same, to make it easier to see
    # what is going on.
    flux = np.ones( integrator.shape )
        
    # Define the flux levels (in DN/s) used in the JPL 3 tests
    dark_dn = 0.0
    very_faint_dn = 0.5
    faint_dn = 2.0
    medium_faint_dn = 5.0
    medium_dn = 20.0
    med_bright_dn = 80.0
    bright_dn = 400.0
    very_bright_dn = 1000.0
    
    # Create a variety of test arrays at different flux levels.
    # Divide by ampgain to convert DN/s to electrons/s
    darkflux = flux * dark_dn / ampgain
    faintflux = flux * faint_dn / ampgain
    mediumflux = flux * medium_dn / ampgain
    medbrightflux = flux * med_bright_dn / ampgain
    brightflux = flux * bright_dn / ampgain
    vbrightflux = flux * very_bright_dn / ampgain
    
    sequence = ''
    time = []
    ctime = []
    signal = []
    if PLOT_LATENT:
        slow_latent = []
# UNCOMMENT FOR EXTRA PLOTTING - ALSO SEE _LATENCY_FUNCTION()
#         m_factor = []
#         s_factor = []
    slope = []
    relslope = []
    firstslope = 0.0
    elapsed = 0.0
    counter = 0

    def test_exposures( name, thisflux, nexposures, nints, ngroups, frame_time ):
        """
        
        Function which encapsulates the repeated code used to test the
        simulator.
        
        """
        global integrator, sequence, time, ctime
        global signal, slow_latent, slope, relslope, elapsed, counter
        global PLOT_LATENT
        print( name )
        sequence += "%s,  " % name
        for expo in range(0,nexposures):
            for integ in range(0,nints):
                if integ == 0:
                    integrator.reset(nresets=3, new_exposure=True)
                else:
                    integrator.reset(new_exposure=False)
                cstart = counter
                for group in range(0,ngroups):
                    integrator.integrate(thisflux, frame_time)
                    array = integrator.readout(nsamples=4)
                    # Multiply by ampgain to convert electrons/s into DN/s
                    array = array * ampgain
                    elapsed += frame_time
                    time.append(elapsed)
                    signal.append(array[2][2])
                    # Needed for extra plotting
                    if PLOT_LATENT:
                        slow_latent.append(integrator.slow_latent[2][2])
# UNCOMMENT FOR EXTRA PLOTTING - ALSO SEE _LATENCY_FUNCTION()
#                         m_factor.append(integrator.m_factor[2][2])
#                         s_factor.append(integrator.s_factor[2][2])
                    del array
                    counter += 1
                m, c = linear_regression(time[cstart:],signal[cstart:])
                #ctime.append(time[cstart]/ 60.0)
                ctime.append(time[cstart])
                slope.append(m)
                if integ == 0 and expo == 0:
                    relslope.append(0.0)
                    firstslope = m
                else:
                    relslope.append(m - firstslope)

    # Simulate the JPL 3 baseline data sequence
    test_exposures('DARK', darkflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('MEDIUM', mediumflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('MED-BRIGHT', medbrightflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('DARK', darkflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    ###test_exposures('FAINT', faintflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('MEDIUM', mediumflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('DARK', darkflux, NDARKEXP, NINTS, NGROUPS, FRAME_TIME)
    ###test_exposures('FAINT', faintflux, NDARKEXP, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('FAINT', faintflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
    test_exposures('DARK', darkflux, NDARKEXP, NINTS, NGROUPS, FRAME_TIME)
#     test_exposures('V-BRIGHT', vbrightflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)
#     test_exposures('DARK', darkflux, NDARKEXP, NINTS, NGROUPS, FRAME_TIME)
#     test_exposures('MEDIUM', mediumflux, NEXPOSURES, NINTS, NGROUPS, FRAME_TIME)

    # Display the ramps, slopes and drift for a particular pixel in
    # the above sequence.
    sequence += "\n%s" % latency_string
    if PLOT_RAMPS:
        mplt.plot_xy(time, signal, linefmt='b.', linestyle='-', figsize=(20,6),
                 xlabel='Time (seconds)',
                 ylabel='Reading (DN)', title="Signal: " + sequence)
    if PLOT_LATENT:
        mplt.plot_xy(time, slow_latent, linefmt='b.', linestyle='-',
                 figsize=(20,6), xlabel='Time (seconds)',
                 ylabel='"Charge" (units)',
                 title="Psuedo trapped charge: " + sequence)
# UNCOMMENT FOR EXTRA PLOTTING - ALSO SEE _LATENCY_FUNCTION()
#         mplt.plot_xy(time, m_factor, linefmt='r.', linestyle='-',
#                  figsize=(20,6), xlabel='Time (seconds)',
#                  ylabel='M_Factor', title="Slow gain multipler: " + sequence)
#         mplt.plot_xy(time, s_factor, linefmt='r.', linestyle='-',
#                  figsize=(20,6), xlabel='Time (seconds)',
#                  ylabel='S_Factor', title="Fast gain multipler: " + sequence)
    if PLOT_SLOPES:
        mplt.plot_xy(ctime, slope, linefmt='b.', linestyle='-', figsize=(20,6),
                 xlabel='Time (seconds)',
                 ylabel='Slope (DN/second)', title="Slope: " + sequence)
        mplt.plot_xy(ctime, relslope, linefmt='b.', linestyle='-', figsize=(20,6),
                 xlabel='Time (seconds)',
                 ylabel='Subtracted Slope (DN/second)',
                 title="Slope relative to first: " + sequence)
    
    print( "\n===Cosmic ray tests." )
    # Dump in some cosmic ray energy
    energy = np.array([[200.0, 3000.0, 200.0],
                       [3000.0, 100000.0, 3000.0],
                       [200.0, 3000.0, 200.0]])
    integrator.hit_by_cosmic_ray(energy, 3, 3)
    integrator.hit_by_cosmic_ray(energy, 5, 5)

    # Make another reading in quick succession.
    # Ensure this reading is not less than the previous one.
    # NOTE: This could be verified using an assert statement.
    readout = integrator.readout()
    print( "Another reading", integrator.readings, \
           "signal (after cosmic ray):\n", readout)
   
    # Reset the integrator and check the signal has returned to zero.
    # NOTE: This could be verified using an assert statement.
    integrator.reset()
    print("\n===Status after a reset.")
    print( integrator )
    # Read out without integrating
    readout = integrator.readout()
    
    # Integrate on darkness
    integrator.integrate(None, 42.0)
    
    # Test the exception raised when the input array has the wrong shape.
#    integrator.integrate([1.0,2.0], 42.0)

    del integrator, flux, energy
    print("Test finished.")
