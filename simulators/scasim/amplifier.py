#!/usr/bin/env python

"""

Module amplifier - Contains the Amplifier class and associated functions.

The simulation is based on the parameters contained in the module
amplifier_properties.py

NOTE: This legacy class is not currently used by SCASim. It may be used
in the future to simulate a detector where the differences between the
zones read out by different amplifiers is significant. MIRI CDPs don't
distinguish between amplifiers and are applied to the whole detector.

:History:

06 Jul 2010: Created
15 Jul 2010: Modified for use with the amplifier_properties module.
20 Jul 2010: Bias levels changed to integer. maxdn attribute added.
             set_readout_mode method added and nsample used to reduce
             the readout noise.
21 Jul 2010: Added plotting methods.
06 Aug 2010: Documentation formatting problems corrected.
11 Aug 2010: Added FITS header generation.
26 Aug 2010: Fixed bug in hit_by_cosmic_ray revealed by unit tests.
01 Sep 2010: Improved attribute checking for public functions.
             Changed cosmic ray simulation algorithm. A strike on an
             amplifier will now create random glitches rather than
             increase the bias level. reset() method added.
14 Sep 2010: Added flag to turn read noise on and off.
30 Sep 2010: Added GAINTYPE.
11 Oct 2010: Added flag to turn amplifier effects on and off.
16 Nov 2010: Some FITS header content rounded to a more sensible number
             of decimal places. Gain calculation carried out in double
             precision. Added function to set the random number generator
             seed, for testing.
23 Nov 2010: Some comments updated. Noise integrity check added to test.
07 Jan 2011: Amplifier properties classified and looked up by SCA_ID
             rather than by FPM_ID.
26 Jan 2010: get_header function converted to set_metadata.
04 Feb 2011: Gain function modified to use numpy.poly1d.
24 Feb 2011: POLYCUBE gain type added. Not yet fully tested.
04 Mar 2011: Documentation tweaks to resolve bad formatting.
22 Mar 2011: Plotting functions modified to use matplotlib figures
             and axes more flexibly (unfortunately the colorbar
             function doesn't work properly).
01 Apr 2011: Documentation formatting problems corrected.
06 Apr 2011: Surplus saturation check (in which a slicing problem
             caused high flux data to appear stripey) removed.
14 Jul 2011: Use of exceptions made more consistent and documented.
14 Sep 2011: Modified to keep up with the separation of
             miri.miritools.miridata into miri.miritools.metadata and
             miri.miritools.miricombdata
20 Oct 2011: Destructor added to remove large objects.
25 Oct 2011: amplifier_properties imported using ParameterFileManager.
             This allows new parameter files to be substituted by
             putting new versions into the current working directory,
             without the need to rebuild the source code.
02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
26 Apr 2012: Brought up to date with changes in Metadata class.
21 May 2012: Plotting modified to use the miriplot utilities.
13 Nov 2012: Major restructuring of the package folder. Import statements
             updated.
14 Jan 2013: olddataproduct subpackage split from dataproduct.
21 May 2013: Removed references to the old MiriMetadata object.
04 Jun 2013: Removed use of noise_mv object. Rationalised attribute names.
13 Jun 2013: Metadata keywords modified to match new data model.
25 Oct 2013: Amplifier properties classified and looked up by DETECTOR
             rather than by SCA_ID.
10 Jun 2014: Corrected typo in metadata keywords.
17 Jun 2014: Changed verbosity definitions and some log messages.
02 Mar 2015: Reformatted some long lines of code.
08 Sep 2015: Made compatible with Python 3
04 Dec 2015: Added a logger.
30 Mar 2015: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
28 Sep 2016: miri.miritools.dataproduct renamed miri.datamodels and
             miri.miritools renamed miri.tools.
14 Feb 2017: Added more noise calculation tests.

@author: Steven Beard (UKATC)

"""
# This module is now converted to Python 3.


# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import os
import random as rn
import numpy as np

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Search for the amplifier parameters file and parse it into a
# properties dictionary. The file is searched for in 3 places:
# (a) The current directory
# (b) The directory where this Python file is being executed from
# (c) The miri.simulators.scasim installation directory.
from miri.tools.filesearching import ParameterFileManager, make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
amplifier_properties = ParameterFileManager(
                            "amplifier_properties.py",
                            search_path=search_path,
                            description="amplifier properties",
                            logger=LOGGER)


class Amplifier(object):
    """
    
    Class Amplifier - Simulates the behaviour of a detector readout amplifier.
    Each amplifier is responsible for a particular portion of the detector
    surface and has its own bias, gain, read noise and non-linearity
    properties. All amplifiers are assumed to share some common electronics
    which contributes an additional bias level which affects all the
    amplifiers.
                 
    :Parameters:
    
    number: int
        Unique amplifier number
    rowmin: int
        The minimum row of pixels for which this amplifier is responsible.
        If zero or None, the amplifier starts with the first row.
    rowmax: int
        The maximum row of pixels for which this amplifier is responsible.
        If None, the amplifier ends with the last row.
    rowstep: int
        The row increment for this amplifier.
        1 or None means every row starting with rowmin up to rowmax.
        2 means every 2nd row starting with rowmin up to rowmax.
        etc...
        The 3 row parameters form a numpy slice, rowmin:rowmax:rowstep.
    colmin: int
        The minimum column of pixels for which this amplifier is responsible.
        If zero or None, the amplifier starts with the first column.
    colmax: int
        The maximum column of pixels for which this amplifier is responsible.
        If None, the amplifier ends with the last column.
    colstep: int
        The column increment for this amplifier.
        1 or None means every column starting with colmin up to colmax.
        2 means every 2nd column starting with colmin up to colmax.
        etc...
        The 3 column parameters form a numpy slice, colmin:colmax:colstep.
    bias: int
        The amplifier bias level in electrons. A constant level added to
        all readouts.
    gain: float, tuple of floats or numpy data cube of floats
        Coefficients describing the amplifier gain. The meaning varies
        according to the gain type (see gaintype parameter). For
        POLYNOMIAL gain type, the gain tuple contains a list of
        coefficients in reverse order, i.e.
        
        * gain[-1] is a constant in DN.
        * gain[-2] is the linear gain in DN/e.
        * gain[-3] is the second order coefficient in DN/e^2
        * gain[-4] is the third order coefficient in DN/e^3
        * etc...
        
        The linear gain is expected to be the dominant term.
        For LINEAR gain type, the gain parameter is a single float
        containing just a linear gain in DN/e.
        For POLYCUBE gain type, the gain parameter is a data cube
        defining a separate set of polynomial coefficients for each
        pixel.
    read_noise: float
        The amplifier read noise in electrons.
    maxdn: float
        The maximum output that this amplifier is capable of generating
        (in DN). If None the output is unlimited.
    gaintype: string, optional, default='POLYNOMIAL'
        The type of gain function. At present only 'LINEAR',
        'POLYNOMIAL' and 'POLYCUBE' are supported.
        NOTE: POLYCUBE has not yet been fully tested.
    simulate_read_noise: boolean, optional, default=True
        A flag that may be used to switch off read noise (for example
        to observe what effects in a simulation are caused by read
        noise).
    simulate_amp_effects: boolean, optional, default=True
        A flag that may be used to switch off amplifier effects, i.e.
        the bias, gain and non-linearity introduced by each amplifier
        (for example to observe what effects in a simulation are caused
        by amplifier effects).
        *Note that when this flag is False the ratio of DNs to electrons
        is exactly 1.0*
    verbose: int, optional, default=2
        Verbosity level. Activates logging statements when non-zero.
        
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
    TypeError
        Raised of any of the initialisation parameters are of the
        wrong type.
        
    """
    # This class variable records any DC reference level which affects
    # all amplifiers. Changes in this level are detected using an extra
    # amplifier reading out only dark reference pixels.
    # Note that each amplifier also has its own bias level, which can be
    # calibrated out using the dark reference columns associated with
    # each amplifier.
    _reference_level = 0
    
    # An array of counters which may be used to keep track of which parts
    # of the detector data have had the amplifier properties applied.
    # The contents should be 0 (not applied yet) or 1 (applied). Any parts
    # which are 2 or larger indicate overlapping regions.
    _applied_array = None
    
    # The number of amplifier objects created.
    _count = 0

    def __init__(self, number, rowmin, rowmax, rowstep, colmin, colmax, colstep,
                 bias, gain, read_noise, maxdn, gaintype='POLYNOMIAL',
                 simulate_read_noise=True, simulate_amp_effects=True,
                 verbose=2, logger=LOGGER):
        """
        
        Constructor for class Amplifier.
        
        Parameters: See class doc string.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("scasim.amplifier")
        self._verbose = int(verbose)
        if verbose > 3:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug( "+++Amplifier object created with" + \
                " row slice=(" + str(rowmin) + "," + \
                                str(rowmax) + "," + \
                                str(rowstep) + ")" + \
                " column slice(=" + str(colmin) + "," + \
                                str(colmax) + "," + \
                                str(colstep) + ")" + \
                "\n\t Bias=" + str(bias) + \
                " Gain=" + str(gain) + \
                " Noise=" + str(read_noise) )
            
        # Note that all parameters (which don't have None as an option)
        # are explicitly converted into the expected data type, since
        # some of the values extracted from the properties dictionary can
        # sometimes be converted by Python into strings, which then upsets
        # formatted I/O statements.
        self.number = int(number)
        self.rowmin = rowmin
        self.rowmax = rowmax
        self.rowstep = rowstep
        self.colmin = colmin
        self.colmax = colmax
        self.colstep = colstep
       
        # This bias is associated with one particular amplifier
        # (c.f. _reference_level which affects all amplifiers).
        try:
            self.bias = int(bias)
        except (ValueError, TypeError) as e:
            strg = "Amplifier bias must be a valid number."
            strg += "\n %s" % e
            raise ValueError(strg)
        
        # Gain type must be one of the recognised options.
        self.gaintype = gaintype
        if self.gaintype == 'LINEAR':
            # For linear gain type the gain must be a single number.
            if isinstance(gain, (float,int)):
                self.gain = gain
                self._polynomial = None # No polynomial function.
            else:
                strg = "Amplifier gain must be a single number "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)
        elif self.gaintype == 'POLYNOMIAL':
            # For polynomial gain type the gain must be a tuple of numbers.
            if isinstance(gain, (tuple,list)):
                self.gain = gain
                # The polynomial function can be prepared in advance
                self._polynomial = np.poly1d(self.gain)
            else:
                strg = "Amplifier gain must be a tuple of numbers "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)
        elif self.gaintype == 'POLYCUBE':
            # A cube of polynomial coefficients. The gain must be a
            # data cube.
            self.gain = np.asarray(gain)
            if self.gain.ndim == 3:
                # The data cube is sliced in exactly the same manner as
                # the detector data.
                # NOTE: The caller must ensure this gain array is of
                # adequate size.
                # The polynomial function cannot be prepared in advance.
                self._polynomial = None
            else:
                strg = "Amplifier gain must be a data cube "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)     
        else:
            strg = "Unrecognised gain function type: %s" % gaintype
            raise ValueError(strg)

        try:
            self.read_noise = float(read_noise)
        except (ValueError, TypeError) as e:
            strg = "Read noise must be a valid number."
            strg += "\n %s" % e
            raise ValueError(strg)

        if maxdn is not None:
            try:
                self.maxdn = float(maxdn)
            except (ValueError, TypeError) as e:
                strg = "If provided, maximum DN must be a valid number."
                strg += "\n %s" % e
                raise ValueError(strg)
        else:
            self.maxdn = None
        self.nsample = 1

        # This flag controls whether read noise is simulated
        self._simulate_read_noise = simulate_read_noise
        if not simulate_read_noise and verbose > 1:
            self.logger.info( "NOTE: Amplifier %d read noise is turned off." % \
                   self._count )

        self._simulate_amp_effects = simulate_amp_effects
        if not simulate_amp_effects and verbose > 1:
            self.logger.info( "NOTE: Amplifier %d bias, gain and " \
                "non-linear gain effects are turned off." % self._count )
        
        # These variables record changes to pixel glitches and noise
        # caused by a cosmic ray hit on the amplifier.
        self._cosmic_ray_energy = 0.0
        self._cosmic_ray_glitches = 0
        self._cosmic_ray_noise = 0.0
        
        Amplifier._count += 1
        self._count = Amplifier._count

    def __del__(self):
        """
        
        Destructor for class Amplifier.
                
        """
        # Explicitly delete large objects created by this class.
        # Exceptions are ignored, since not all the objects here
        # will exist if the destructor is called as a result of
        # an exception.
        try:
            # Objects created in the constructor
            if self._polynomial is not None:
                del self._polynomial
        except Exception:
            pass

    def set_metadata(self, metadata):
        """
        
        Add amplifier keywords to the given MIRI metadata.
        
        :Parameters:
        
        metadata: dictionary-like object
            A keyword-addressable metadata object to which amplifier
            keywords should be added.
            
        :Returns:
        
        metadata: dictionary-like object
            An updated metadata object.
            
        """
        noise = 0.0
        count = 0
        num = "%1d" % self.number
        if self._simulate_read_noise:
            # Round the value to 4 decimal places
            metadata["RDNOISE" + num] = float("%.4f" % self.read_noise)
        else:
            metadata["RDNOISE" + num] = 0.0

        if self._simulate_amp_effects:
            metadata["RDBIAS" + num] =  self.bias
            metadata["GAINFN"] = self.gaintype
            
            if self.gaintype == 'LINEAR':
                metadata["GAINCF" + num] = str(self.gain)
            elif self.gaintype == 'POLYNOMIAL':
                metadata["GAINCF" + num] = str(self.gain)
            elif self.gaintype == 'POLYCUBE':
                # TODO: Write gain file name to FITS header?
                pass
            metadata["MAXDN" + num] = self.maxdn
        else:
            metadata["RDBIAS" + num] = 0.0
            metadata["GAINFN"] = 'NONE'
            metadata["GAINCF" + num] = str(1.0)
            metadata["MAXDN" + num] = 0
        
        return metadata

    def set_readout_mode(self, nsample):
        """
        
        Defines the SCA detector readout mode.
        
        :Parameters:
        
        nsample: int
            The number of (non-discarded) samples included per readout
            (which affects the read noise). Must be at least 1.
        
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
            
        """
        self.nsample = int(nsample)
        if self.nsample <= 0:
            strg = "Number of samples per readout must be at least 1."
            raise ValueError(strg)
        
    def set_seed(self, seedvalue=None):
        """
        
        Set the seed for the numpy random number generator.
        This function can be used while testing to ensure the
        randomised amplifier effects that follow are well defined.
        
        :Parameters:
        
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator.
            If not specified, a value of None will be sent, which
            randomises the seed.
            
        """
        np.random.seed(seedvalue)
        
    def _set_bias(self, bias):
        """
        
        Test function for changing the bias.
                 
        :Parameters:
    
        bias: int
            The amplifier bias level in electrons. A constant level
            added to all readouts.
        
        """
        self.bias = int(bias)

    def _set_gain(self, gain, gaintype='POLYNOMIAL'):
        """
        
        Test function for changing the gain.
                 
        :Parameters:
    
        gain: float, tuple of floats or numpy data cube of floats
            Coefficients describing the amplifier gain. The meaning varies
            according to the gain type (see gaintype parameter). For
            POLYNOMIAL gain type, the gain tuple contains a list of
            coefficients in reverse order, i.e.
        
            * gain[-1] is a constant in DN.
            * gain[-2] is the linear gain in DN/e.
            * gain[-3] is the second order coefficient in DN/e^2
            * gain[-4] is the third order coefficient in DN/e^3
            * etc...
        
            The linear gain is expected to be the dominant term.
            For LINEAR gain type the gain parameter is a single float
            containing just a linear gain in DN/e.
            For POLYCUBE gain type, the gain parameter is a data cube
            defining a separate set of polynomial coefficients for each
            pixel.
        gaintype: string, optional, default='POLYNOMIAL'
            The type of gain function. At present only 'LINEAR',
            'POLYNOMIAL' and 'POLYCUBE' are supported.
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
        TypeError
            Raised if any of the parameters are of the wrong type,
            size or shape.
        
        """
        # Gain type must be one of the recognised options.
        self.gaintype = gaintype
        # TODO: Additional GAINTYPEs?.
        if self.gaintype == 'LINEAR':
            # For linear gain type the gain must be a single number.
            if isinstance(gain, (float,int)):
                self.gain = gain
                if self._polynomial is not None:
                    del self._polynomial # Avoid memory leak.
                self._polynomial = None # No polynomial function.
            else:
                strg = "Amplifier gain must be a single number "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)
        elif self.gaintype == 'POLYNOMIAL':
            # For linear gain type the gain must be a tuple of numbers.
            if isinstance(gain, (tuple,list)):
                self.gain = gain
                # The polynomial function needs to be updated
                if self._polynomial is not None:
                    del self._polynomial # Avoid memory leak.
                self._polynomial = np.poly1d(self.gain)
            else:
                strg = "Amplifier gain must be a tuple of numbers "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)
        elif self.gaintype == 'POLYCUBE':
            # A cube of polynomial coefficients. The gain must be a
            # data cube.
            gain = np.asarray(gain)
            if gain.ndim == 3:
                # The data cube is sliced in exactly the same manner as
                # the detector data.
                rmin = self.rowmin
                rmax = self.rowmax
                rstep = self.rowstep
                cmin = self.colmin
                cmax = self.colmax
                cstep = self.colstep
                self.gain = gain[:,rmin:rmax:rstep,cmin:cmax:cstep]
                # The slice must be inside the data cube.
                if (self.gain.shape[1] > 0) and (self.gain.shape[2] > 0):
                    # The polynomial function cannot be prepared in advance
                    self._polynomial = None
                else:
                    strg = "Amplifier gain data cube is too small."
                    raise TypeError(strg)     
            else:
                strg = "Amplifier gain must be a data cube "
                strg += "for %s gain type." % self.gaintype
                raise TypeError(strg)     
        else:
            strg = "Unrecognised gain function type: %s" % gaintype
            raise ValueError(strg)

    def _set_read_noise(self, read_noise):
        """
        
        Test function for changing the read noise.
        Unlike "read_noise_off" and "read_noise_on", this function
        overwrites and loses the current read noise. It's more specific
        but less reversible.
                 
        :Parameters:
    
        noise: float
            The amplifier read noise in electrons.
        
        """
        self.read_noise = float(read_noise)
        
    def reset(self):
        """
        
        Reset all cosmic ray effects.
        
        :Parameters:
        
        None.
        
        """
        self._cosmic_ray_energy = 0.0
        self._cosmic_ray_glitches = 0
        self._cosmic_ray_noise = 0.0        
        
    def hit_by_cosmic_ray(self, energy):
        """
        
        Simulate a cosmic ray hit on the amplifier.
        
        :Parameters:
        
        energy: float
            The energy of the cosmic ray (in terms of the number of electrons
            it would have released in the detector pixels).
            
        """
        if self._verbose > 3:
            self.logger.debug( "Cosmic ray hit with energy %f." % energy )
        
        # The effect of a cosmic ray hit on an amplifier is not
        # fully known but, according to Mike Ressler, a cosmic ray
        # strike happening just as the A/D converters are processing
        # a pixel might cause a glitch in the readout of a single
        # pixel. I am also assuming that a small amount of heating
        # may occur, which might slightly increase the read noise.
        self._cosmic_ray_energy += float(energy)
        self._cosmic_ray_glitches += \
            int(energy * amplifier_properties['COSMIC_RAY_GLITCH_SENSITIVITY'])
        self._cosmic_ray_noise += \
            energy * amplifier_properties['COSMIC_RAY_NOISE_SENSITIVITY']

    def readout(self, detector_data, reset_counter=False, removeneg=True):
        """
        
        Apply the bias, gain and readout noise to the portion of the
        detector data associated with this amplifier. Other areas of
        the detector data are unaffected.
        
        *NOTE: The detector data will not be complete until the readout
        method from all the associated amplifiers has been applied to
        it. The readout counter can be used to keep track of how many
        times readout effects have been applied to different parts of
        the detector data.*
        
        :Parameters:
        
        detector_data: array_like
            The detector data on which to apply the readout.
            *NOTE: The array is converted to floating point. An integer
            array can wrap around and generate spurious large values
            if the read noise makes it go negative.*
        reset_counter: bool, optional, default=False
            Force a reset of the counter recording which areas of the
            detector have had the amplifier properties applied. Otherwise
            the counter is reset whenever the detector_data array changes
            shape.
        removeneg: bool, optional, default=True
            If True, remove negative values from the readout and replace
            them with zero.
                
        :Returns:
        
        detector_data: array_like
            Modified detector data. Same shape as the input data.
            
        """
        # Initialise the counter if it does not yet exist.
        # Otherwise reset it if the detector data changes shape
        # or a reset has been forced.
        detector_data = np.asarray(detector_data, dtype=np.float32)
        if Amplifier._applied_array is None:
            Amplifier._applied_array = np.zeros(detector_data.shape, \
                                                dtype=np.int)
        elif Amplifier._applied_array.shape != detector_data.shape:
            del Amplifier._applied_array
            Amplifier._applied_array = np.zeros(detector_data.shape, \
                                                dtype=np.int)
        elif reset_counter:
            Amplifier._applied_array *= 0
            
        # Copy the row and columns associated with this amplifier
        # to shorter variable names, to make the slicing easier.
        rmin = self.rowmin
        rmax = self.rowmax
        rstep = self.rowstep
        cmin = self.colmin
        cmax = self.colmax
        cstep = self.colstep
            
        # Take a random sample of numbers from a normal distribution
        # with zero mean and unit variance and distribute them over
        # the detector pixels.
        randArray = np.random.randn(detector_data.shape[0],
                                    detector_data.shape[1]).astype(np.double)
                
        # Determine the amount of read noise to by applied to the detector
        # data. The noise is reduced when there is more than one sample.
        # NOTE: The cosmic ray noise is not added as a root mean square
        # because it is not an independent source of noise - it's a
        # temporary boost to the same read noise.
        if self._simulate_read_noise:
            total_noise = self.read_noise + self._cosmic_ray_noise
            if self.nsample > 1:
                if self._verbose > 3:
                    self.logger.debug( \
                        "Applying readout noise %fe sampled %d times." % \
                        (total_noise, self.nsample) )
                noise = total_noise / np.sqrt(self.nsample)
            else:
                if self._verbose > 3:
                    self.logger.debug( "Applying readout noise %fe." % total_noise )
                noise = total_noise
        else:
            noise = 0.0

        # If there has been a cosmic ray strike, insert a streak of
        # bad values at a random point within randArray. If that
        # portion of the array happens to fall into the amplifier zone
        # it will greatly boost the noise and generate a glitch.
        if self._cosmic_ray_glitches > 0:
            row = int(np.random.uniform(0, detector_data.shape[0]))
            startcol = int(np.random.uniform(0, detector_data.shape[1]))
            endcol = startcol + (self._cosmic_ray_glitches*4) - 1
            if endcol >= detector_data.shape[1]:
                endcol = detector_data.shape[1] - 1
            # Set the initial glitch magnitude arbitrarily to roughly
            # the number of electrons displaced in a detector pixel.
            glitch = self._cosmic_ray_energy / (1.0 + noise)
            # Generate a streak of bad pixels along the same row.
            for col in range(startcol, endcol+1, 4):
                randArray[row, col] *= glitch
                glitch *= amplifier_properties['COSMIC_RAY_GLITCH_IPC']
        
        # Apply the electronic bias and read noise to the appropriate
        # slice of the detector data. NOTE: Although the Poisson noise
        # and read noise are not added explicitly in quadrature,
        # combining the two sets of randomly generated offsets has the
        # same effect.
        if self._simulate_amp_effects:
            detector_data[rmin:rmax:rstep, cmin:cmax:cstep] += \
                Amplifier._reference_level + \
                self.bias + \
                (randArray[rmin:rmax:rstep, cmin:cmax:cstep] * noise)
        else:
            detector_data[rmin:rmax:rstep, cmin:cmax:cstep] += \
                (randArray[rmin:rmax:rstep, cmin:cmax:cstep] * noise)
            
        del randArray
        
        # The smallest DN is zero, so the read noise shouldn't make the
        # reading go negative. If requested, replace negative values
        # with zero.
        if removeneg:
            arenegative = np.where(detector_data < 0.0)
            detector_data[arenegative] = 0.0
        
        # Apply the gain to convert electrons into DNs
        if self._simulate_amp_effects:
            detector_data[rmin:rmax:rstep, cmin:cmax:cstep] = \
                self._applygain(detector_data[rmin:rmax:rstep, \
                                              cmin:cmax:cstep])
        
        # The smallest DN is zero, so the gain shouldn't make the reading
        # go negative. If requested, replace negative values with zero.
        if removeneg:
            arenegative = np.where(detector_data < 0.0)
            detector_data[arenegative] = 0.0

        # Increment the counter.
        Amplifier._applied_array[rmin:rmax:rstep, cmin:cmax:cstep] += 1
        
        # Assume that any glitches caused by cosmic ray effects go away
        # after a readout and the excess noise gradually decays away.
        self._cosmic_ray_energy = 0.0
        self._cosmic_ray_glitches = 0
        if self._cosmic_ray_noise > 0.0:
            self._cosmic_ray_noise *= \
                amplifier_properties['COSMIC_RAY_NOISE_DECAY']
            if self._cosmic_ray_noise < 0.001:
                # An excess noise less than 0.001 electrons makes no
                # detectable difference. Stop the decay process.
                self._cosmic_ray_noise = 0.0

        return detector_data
    
    def _applygain(self, data):
        """
        
        Apply the gain function to a set of data
        
        :Parameters:
        
        data: array_like
            Array to which the gain function is to be applied.
            
        :Returns:
        
        modified_data: array_like
            Array of the same shape as data with gain function applied.
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.
            
        """
        # TODO: Additional GAINTYPEs?.
        if self.gaintype == 'LINEAR':
            # Linear gain.
            data *= self.gain
        elif self.gaintype == 'POLYNOMIAL':
            # Polynomial gain. Use the polynomial function prepared earlier
            # or create one now.
            if self._polynomial is None:
                self._polynomial = np.poly1d(self.gain)
            temp_data = self._polynomial(data)
            # NOTE: The following line only works when self._polynomial is
            # a numpy ndarray
            data = temp_data.astype(np.float)
        elif self.gaintype == 'POLYCUBE':
            # Polynomial gain with a different set of coefficients for
            # each pixel. NOTE: This calculation will fail if the gain
            # array does not have the correct shape to be broadcast
            # onto the data array.
            try:
                temp_data = self.gain[0,:,:]
                for xx in range(1, self.gain.shape[0]):
                    temp_data *= data
                    temp_data += self.gain[xx,:,:]
                data = temp_data.astype(np.float)
            except (ValueError, TypeError, IndexError) as e:
                strg = "POLYCUBE gain calculation failed.\n   %s\n" % e
                strg += "   Gain cube probably has the wrong shape: %s " % \
                    self.gain.shape.__str__()
                strg += "compared with data shape %s." % data.shape.__str__()
                raise ValueError(strg)
            
        # Clip the DN at its largest allowed value, if one has been provided.
        if self.maxdn is not None:
            aresaturated = np.where(data > self.maxdn)
            data[aresaturated] = self.maxdn

        return data
 
    def __str__(self):
        """
        
        Return a string describing an Amplifier object.
        
        """
        strg = "Readout amplifier %d out of %d " % \
            (self._count, Amplifier._count)
        if Amplifier._applied_array is not None:
            strg += "reading a detector of %d rows x %s columns" % \
                Amplifier._applied_array.shape
        strg += "\n      covering "
        strg += self._slice_str("row", \
                               self.rowmin, self.rowmax, self.rowstep)
        strg += " and "
        strg += self._slice_str("column", \
                               self.colmin, self.colmax, self.colstep)
        # TODO: Additional GAINTYPEs?.
        if self.gaintype == 'LINEAR':
            strg += ".\n      bias=%de, " \
                "gain=%f DN/e, noise=%.1fe (%d samples) " % \
                (self.bias, self.gain, self.read_noise, self.nsample)
        elif self.gaintype == 'POLYNOMIAL':
            strg += ".\n      bias=%de, " \
                "gain(DN)=%s (e), noise=%.1fe (%d samples) " % \
                (self.bias, self._poly_str(self.gain), self.read_noise,
                 self.nsample)
        elif self.gaintype == 'POLYCUBE':
            strg += ".\n      bias=%de, " \
                "gain(DN)=POLYCUBE%s, noise=%.1fe (%d samples) " % \
                (self.bias, self.gain.shape.__str__(), self.read_noise, \
                 self.nsample)

        if self.maxdn is not None:
            strg += "with max DN=%.0f " % self.maxdn
        strg += "and electronic bias=%de." % Amplifier._reference_level
        
        if Amplifier._applied_array is not None:
            strg += "\n      Readout counter ranges from %d to %d." % \
                (Amplifier._applied_array.min(), Amplifier._applied_array.max())
        return strg
    
    def _slice_str(self, description, start, end, step):
        """
        
        Returns a string describing the given slice
        
        """
        if start is None and end is None:
            if step == 1 or step is None:
                strg = "all " + description + "s"
            else:
                if step > 10 and step < 20:
                    ordinal = "th"
                elif step % 10 == 1:
                    ordinal = "st"
                elif step % 10 == 2:
                    ordinal = "nd"
                elif step % 10 == 3:
                    ordinal = "rd"
                else:
                    ordinal = "th"
                strg = "every %d%s %s" % (step, ordinal, description)
        else:            
            if start is None or start == 0:
                strg = "from the first " + description
            else:
                strg = "from %s %d" % (description, start)
            if end is None:
                strg += " to the last " + description
            else:
                strg += " to %s %d" % (description, end)
            
            if step != 1 and step is not None:
                strg += " stepping by %d" % step 
            
        return strg
    
    def _poly_str(self, coeffs):
        """
        
        Returns a string describing the polynomial coefficients
        
        """
        strg = "("
        if coeffs:
            if isinstance(coeffs, (float,int)):
                strg += "%g" % coeffs
            elif len(coeffs) == 2:
                strg += "%g x + %g" % coeffs
            else:
                for xx in range(0, len(coeffs)-1):
                    order = len(coeffs) - xx - 1
                    if order > 1:
                        strg += "%g x^%d + " % (coeffs[xx], order)
                    else:
                        strg += "%g x + " % coeffs[xx]
                strg += "%g" % coeffs[-1]
        strg += ")"
        return strg

    def plotgain(self, xmin=0, xmax=100000, npoints=512, title=None,
                  description=''):
        """
        
        Plot amplifier gain as a function of number of electrons
        to a self-contained matplotlib figure.
        
        :Parameters:
        
        xmin: int, optional, default=0
            The minimum number of electrons to be plotted.
        xmax: int, optional, default=100000
            The maximum number of electrons to be plotted.
        npoints: int, optional, default=512
            The number of points to be plotted.
        title: string, optional
            Optional title for the plot. The default is a string
            describing the Amplifier object.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
        description: string, optional
            Additional description to be shown on the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # First generate the gain function
        if self.gaintype == 'POLYCUBE':
            # TODO: 3-D gains cannot be plotted with this function
            strg = "POLYCUBE gain not supported by plotgain function"
            raise ValueError(strg)
        else:
            electrons = np.linspace(xmin, xmax, npoints)
            dns = self._applygain(electrons)

        # Create a figure to contain the plot.
        fig = mplt.new_figure(1, figsize=(10,6))
        ax = mplt.add_subplot(fig, 1, 1, 1)

        # Construct a plot title.
        if title is None:
            tstrg = "Gain function for amplifier %d" % self._count
        else:
            tstrg = title
        if description:
            tstrg += " %s" % description
   
        # Plot the actual gain using a blue line.
        ax = mplt.plot_xy(electrons, dns, plotfig=fig, plotaxis=ax,
                          xlabel='Number of electrons', ylabel='DN',
                          title=tstrg, linefmt='b-', linestyle='-')
        # If the gain is not linear, superimpose a perfect
        # linear gain using red dots.
        if self.gaintype == 'POLYNOMIAL' and len(self.gain) > 1:
            ax = mplt.plot_xy(electrons, electrons * self.gain[-2],
                              plotfig=fig, plotaxis=ax, linefmt='r-',
                              linestyle=':')                

        if self._verbose > 1:
            strg = "Plotting %s. " % tstrg
            strg += "Close the plot window to continue."
        else:
            strg = ''
        # Close and display the plot.
        mplt.show_plot(prompt=strg)

    def plot(self, plotfig=None, plotaxis=None, labels=True, withbar=True,
             description='', **kwargs):
        """
        
        Plot the amplifier readout count within the given matplotlib axis.
        This function can can be used to include a plot of this object
        in any figure.
        The plotfig method can be used to create a self-contained
        plot.

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
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        withbar: bool, optional
            Set to False to suppress the colour bar. The default is True.
        description: string, optional
            Optional description to be added to the title, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # Plot the amplifier readout count as an image with the origin at
        # the bottom.
        # NOTE: plotaxis.imshow() will plot the image but will cause the
        # colorbar() function to fail.
        if Amplifier._applied_array is not None:
            tstrg = "Count of readouts"
            if description:
                tstrg += " %s" % description
            if labels:
                xlabel = 'Columns'
                ylabel = 'Rows'
            else:
                xlabel = ''
                ylabel = ''
            mplt.plot_image(Amplifier._applied_array, plotaxis=plotaxis,
                            datatype='matrix', cmap='hot', withbar=withbar,
                            xlabel=xlabel, ylabel=ylabel, title=tstrg,
                            **kwargs)
        else:
            self.logger.info( "*** No readouts so no readout counter to plot." )


def set_reference_level(newlevel):
    """
    
    Sets a new electronic reference level level which affects all amplifiers.
        
    :Parameters:
        
    newlevel: int
        A new reference level level in electrons.
        
    """
    Amplifier._reference_level = int(newlevel)

def random_level_change(maxchange):
    """
    
    Make a random change to the electronic reference level up to the
    given limit, but don't let the level go below zero. This simulates
    random drifts that may take place within the electronics and can
    be used to verify that the pipeline software can remove such drifts
    by examining the reference pixels.
        
    :Parameters:
        
    maxchange: int
        The maximum change to be applied to the reference level in electrons.
        
    """
    change = (2.0 * maxchange * (rn.random() - 0.5)) + 0.5
    
    Amplifier._reference_level += int(change)
    if Amplifier._reference_level < 0:
        Amplifier._reference_level = 0
        
def check_readout():
    """
    
    Check that the readout counters are consistent after reading out
    all the amplifiers. All the counters should be the same, otherwise
    some parts of the detector have been missed or other parts have
    been overlapped and had the readout parameters applied more than once.
    
    :Parameters:
    
    None.
    
    :Returns:
    
    ok: bool
        True if the readout counters are consistent. False if the
        readout is unfinished or has overlapped.
        
    """
    return Amplifier._applied_array.min() == Amplifier._applied_array.max()


#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_quantum_efficiency module.
# However, the tests here are made in verbose mode and include plotting,
# so they show what is happening in more detail.
#
if __name__ == '__main__':
    print( "Testing the Amplifier class." )
    
    TEST_AMPLIFER = False
    TEST_NOISE = True
    PLOTTING = True        # Set to False to turn off plotting.
    VERBOSE = 3
    
    if TEST_AMPLIFER:
        gain0 = 0.5
        gain1 = amplifier_properties.get('_amp493_1', 'GAIN')
#        gain2 = (1e-18, 4.7e-13, -2e-7, 0.15, 1.0)
        # An example of how to set up a test 3-D gain cube.
        carr0 = np.ones((8,8)) * gain1[0]
        carr1 = np.ones((8,8)) * gain1[1]
        carr2 = np.ones((8,8)) * gain1[2]
        gain3 = np.asarray( (carr0.tolist(), carr1.tolist(), carr2.tolist()) )
      
        amp0 = Amplifier(0, None,None,None, 0,None,2, 1, gain0, 0.1, None,
                         gaintype='LINEAR', verbose=VERBOSE)
        amp1 = Amplifier(1, None,None,3, 1,None,2, 2, gain1, 0.1, 65535,
                         gaintype='POLYNOMIAL', verbose=VERBOSE)
        amp1.set_readout_mode(10)
#        amp2 = Amplifier(1, None,None,3, 1,None,2, 2, gain2, 0.1, 65535,
#                         gaintype='POLYNOMIAL', verbose=2)
        amp3 = Amplifier(1, None,None,3, 1,None,2, 2, gain3, 0.1, 65535,
                         gaintype='POLYCUBE', verbose=VERBOSE)
    
        print( amp0 )
#        amp0.plotgain()
        print( amp1 )
        if PLOTTING:
            amp1.plotgain(xmax=250000)
            amp1.plot()
#        print( amp2 )
#        if PLOTTING:
#            amp2.plotgain(xmax=250000)
#            amp2.plot()
        print( amp3 )
        
        set_reference_level( 42 )
        random_level_change( 3 )
        print( amp0 )
        print( amp1 )
    
    if TEST_NOISE:
        print( "\nVerification of the integrity of the noise calculations:" )
        import scipy.stats
        import miri.tools.miriplot as mplt
        
        # Create a 10000 element array containing a flat illumination of 100
        # electrons per second.
        SIZE = 10000
        SHAPE = [100,100]
        fluxlist = [100.0, 1000.0, 10000.0, 40000.0, 80000.0]
        for FLUX in fluxlist:
            e = FLUX + np.zeros([SIZE])
            e.shape = SHAPE
            print( "Simulated flux=%f" % FLUX)
            print( "Expected:\t mean=%f, standard deviation=%f, variance=%f" % \
                (np.mean(e), np.std(e), np.var(e)) )
            if PLOTTING:
                # Plot the ramp as an image plot.
                title = "Flux=%f: Expected count" % FLUX
                mplt.plot_image(e, title=title, xlabel='X', ylabel='Y')
            
            # Add Poisson noise to this flat array. The standard deviation should
            # be around sqrt(100) or 10.
            p = scipy.stats.poisson.rvs(e)
            print( "Poisson noise:\t mean=%f, standard deviation=%f, variance=%f" % \
                (np.mean(p), np.std(p), np.var(p)) )
            if PLOTTING:
                # Plot the ramp as an image plot.
                title = "Flux=%f: Poisson noise added" % FLUX
                mplt.plot_image(p, title=title, xlabel='X', ylabel='Y')
                
            # Add some read noise of 20 electrons by applying random adjustments
            # scattered in a normal distribution.
            RNOISE = 25.0
            rand = RNOISE * np.random.randn(SIZE)
            rand.shape = SHAPE
            r = e + rand
            print( "Read noise:\t mean=%f, standard deviation=%f, variance=%f" % \
                (np.mean(r), np.std(r), np.var(r)) )
            if PLOTTING:
                title = "Flux=%f: Read noise added" % FLUX
                mplt.plot_image(r, title=title, xlabel='X', ylabel='Y')
                
            # Now combine the Poisson noise and read noise together
            pr = p + rand
            print( "Poisson plus read noise: " + \
                "mean=%f, standard deviation=%f, variance=%f" % \
                (np.mean(pr), np.std(pr), np.var(pr)) )
            quadrature = np.sqrt(FLUX + (RNOISE * RNOISE))
            strg = "Compare the simulated noise, %f, " % np.std(pr)
            strg += "with the noise added in quadrature, %f" % quadrature
            print( strg )
            if PLOTTING:
                title = "Flux=%f: Poisson + Read noise added" % FLUX
                mplt.plot_image(pr, title=title, xlabel='X', ylabel='Y')
            del e, p, r, pr

        # Create a 40 element ramp.
        SIZE = 40
        fluxlist = [100.0, 1000.0, 10000.0, 40000.0, 80000.0]
        for FLUX in fluxlist:
            e = np.linspace(10.0, FLUX, SIZE)
            stime = 2.75
            tmin = 0.0
            tmax = SIZE * stime
            time = np.linspace(tmin, tmax, SIZE)            
            if PLOTTING:
                # Plot the ramp as an XY plot.
                title = "Flux=%f: Expected count" % FLUX
                mplt.plot_xy(time, e, linefmt='bo', xlabel='frame',
                             ylabel='electrons', title=title)
            
            
            # Add Poisson noise to this ramp array.
            p = scipy.stats.poisson.rvs(e)
            if PLOTTING:
                # Plot the ramp as an XY plot.
                title = "Flux=%f: Poisson noise added" % FLUX
                mplt.plot_xy(time, p, linefmt='bo', xlabel='frame',
                             ylabel='electrons', title=title)
                
            # Add some read noise of 20 electrons by applying random adjustments
            # scattered in a normal distribution.
            RNOISE = 20.0
            rand = RNOISE * np.random.randn(SIZE)
            r = e + rand
            if PLOTTING:
                # Plot the ramp as an XY plot.
                title = "Flux=%f: Read noise added" % FLUX
                mplt.plot_xy(time, r, linefmt='bo', xlabel='frame',
                             ylabel='electrons', title=title)
                
            # Now combine the Poisson noise and read noise together
            pr = p + rand
            if PLOTTING:
                # Plot the ramp as an XY plot.
                title = "Flux=%f: Poisson + Read noise added" % FLUX
                mplt.plot_xy(time, pr, linefmt='bo', xlabel='frame',
                             ylabel='electrons', title=title)
            del e, p, r, pr
    
    print( "\nTest finished.")
