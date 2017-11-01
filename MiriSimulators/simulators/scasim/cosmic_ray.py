#!/usr/bin/env python

"""

Module cosmic_ray - Contains the CosmicRayEnvironment and CosmicRay 
classes and associated functions.

:History:

25 Jul 2010: Created
30 Jul 2010: Used np.asarray to convert arrays passed as arguments.
06 Aug 2010: Documentation formatting problems corrected.
11 Aug 2010: Added FITS header generation.
26 Aug 2010: Corrected a bug in the loop within generate_event.
01 Sep 2010: Improved attribute checking for public functions.
13 Sep 2010: Check for the location of the test file.
06 Oct 2010: Changed to new cosmic ray library for SIA detectors.
12 Oct 2010: Corrected typo in FITS header.
16 Nov 2010: Some FITS header content rounded to a more sensible number
             of decimal places. Added function to set the random number
             generator seed, for testing.
25 Nov 2010: Don't report cosmic ray events when there are no hits
             expected.
26 Jan 2010: get_header function converted to set_metadata.
25 Mar 2011; Use matplotlib.matshow to display cosmic ray hit maps.
01 Apr 2011: Documentation formatting problems corrected.
14 Jul 2011: Use of exceptions made more consistent and documented.
14 Sep 2011: Modified to keep up with the separation of
             miritools.miridata into miritools.metadata and
             miritools.miricombdata
25 Oct 2011: cosmic_ray_properties imported using ParameterFileManager.
             This allows new parameter files to be substituted by
             putting new versions into the current working directory,
             without the need to rebuild the source code.
             Bugs corrected in exception handling of CosmicRayEnvironment
             constructor and in load_cosmic_ray_random().
11 Jan 2012: Renamed symbols to avoid potential name clash: map-->crmap
             and sum-->esum
02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
21 May 2012: Plotting modified to use the miriplot utilities.
21 May 2013: Removed references to the old MiriMetadata object.
05 Jun 2013: Rationalised attribute names.
17 Jun 2014: Changed verbosity definitions.
02 Mar 2015: Reformatted some long lines of code.
21 May 2015: Amplifier level removed and amplifier cosmic ray effects
             no longer simulated.
27 May 2015: Replaced pyfits with astropy.io.fits
08 Sep 2015: Made compatible with Python 3
04 Dec 2015: Added a logger.
01 Mar 2015: Added function which combines cosmic ray libraries into one
             large dataset. "summed" option added to histogram plots and
             plot labels clarified. Logarithmic plots made the default.
30 Mar 2015: Reduced the scope of the parameter file search so it only
             checks in 3 directories.
28 Sep 2016: miri.miritools.dataproduct renamed miri.datamodels and
             miri.miritools renamed miri.tools.
09 Nov 2016: Optionally, convolve the cosmic ray maps with the detector
             coupling function and add some blurring.
14 Mar 2017: Added the ability to rebin the cosmic ray energy maps
             when CR_BINNING_FACTOR is greater than 1.
05 Apr 2017: Switch to the newer version of Massimo's cosmic ray library,
             for a 35 micron detector thickness and installed locally.
             Corrected bug in hist_electrons plotting function.
07 Jul 2017: Added energy and electron scaling factors.
17 Jul 2017: Added rebinning function to load_cosmic_ray_libraries.
             Added flags for finer control over which tests are run.
             Statistics function added. Rebin by averaging.
31 Oct 2017: Reduced verbosity of cosmic ray event reporting.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Python logging facility
import logging
logging.basicConfig(level=logging.INFO) # Default level is informational output 
LOGGER = logging.getLogger("miri.simulators") # Get a default parent logger

import math, os
#import random as rn
import numpy as np
from scipy.signal import convolve2d
from scipy.ndimage.filters import gaussian_filter
import astropy.io.fits as pyfits

# Import the miri.tools plotting module.
import miri.tools.miriplot as mplt

# Search for the cosmic ray parameters file and parse it into a
# properties dictionary. The file is searched for in 3 places:
# (a) The current directory
# (b) The directory where this Python file is being executed from
# (c) The miri.simulators.scasim installation directory.
from miri.tools.filesearching import ParameterFileManager, make_searchpath
import miri.simulators.scasim
dir_list = ['.', os.path.dirname(__file__), miri.simulators.scasim.__path__[0]]
search_path = make_searchpath(dir_list)
cosmic_ray_properties = ParameterFileManager(
                            "cosmic_ray_properties.py",
                            search_path=search_path,
                            description="cosmic ray properties",
                            logger=LOGGER)

def _rebin(arr, new_shape, use_maximum=False, conserve_e=False):
    """
    
    Helper function which rebins a 2-D array to a new shape.
    A new array is returned with each new pixel containing the mean or
    the maximum value found in the contributing pixels.
    
    The conserve_e flag can be set to True to conserve the total electron
    count or to False to preserve the maximum electron count per pixel.
    True gives a better result when the detector pixels are larger than
    assumed by the cosmic ray library.
    False gives a better result when the detector thickness is smaller than
    assumed by the cosmic ray library.
    
    """
    if conserve_e:
        oldsum = arr.sum()
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    if use_maximum:
        newarr = arr.reshape(shape).max(-1).max(1)
    else:
        newarr = arr.reshape(shape).max(-1).mean(1)
    # If requested, conserve the total number of electrons within the image.
    if conserve_e:
        newsum = newarr.sum()
        newarr = newarr * oldsum / newsum
    return newarr     


class CosmicRay(object):
    """
    
    Class CosmicRay - Describes the properties and the effect of a
    single cosmic ray event on the MIRI SCA.
    
    :Parameters:
    
    energy: float
        The total energy of the cosmic ray event in MeV.
    target_coords: int or tuple of 2 ints
        The coordinates of the detector pixel hit by the cosmic ray
        as (row, column).
    hit_map: array_like
        A map showing the distribution of cosmic ray energy between the
        central pixel and its neighbours in electrons.
    nucleon: string, optional
        A description of the nucleon responsible for this cosmic ray.
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
    TypeError
        Raised if any of the initialisation parameters are of the
        wrong type, size or shape.

    """
    
    def __init__(self, energy, target_coords, hit_map, nucleon='', verbose=2,
                 logger=LOGGER):
        """
        
        Constructor for class CosmicRay.
        
        Parameters: See class doc string.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("scasim.cosmic_ray")
        self._verbose = int(verbose)
        if verbose > 4:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug( "+++CosmicRay object created with energy " + \
                str(energy) + " MeV hitting the detector at " + \
                str(target_coords))
            
        # Note that all parameters (which don't have None as an option)
        # are explicitly converted into the expected data type, since
        # some of the values extracted from the properties dictionary can
        # sometimes be converted by Python into strings, which then upsets
        # formatted I/O statements.
        try:
            self.energy = float(energy)
        except (ValueError, TypeError) as e:
            strg = "Cosmic ray energy must be a valid number."
            strg += "\n %s" % e
            raise ValueError(strg)
        
        if isinstance(target_coords, (tuple,list)) and len(target_coords) == 2:
            self.target_coords = target_coords
        else:
            strg = "Detector target coordinates " + \
                    "must be a tuple of 2 numbers"
            raise TypeError(strg)
        
        self.hit_map = np.asarray(hit_map)
        self.nucleon = nucleon

    def get_target_coords(self):
        """
        
        Return the target coordinates in detector pixels.

        :Returns:

        target_coords: tuple of  2 ints
            The coordinates of the detector pixel hit by the cosmic ray
            as (row, column).
        
        """
        return self.target_coords
    
    def get_energy(self):
        """
        
        Return the cosmic ray energy in Mev.

        :Returns:

        energy: float
            The total energy of the cosmic ray event in MeV.
        
        """
        return self.energy
    
    def get_electrons(self, scale=1.0):
        """
        
        Return the total number of electrons released by the cosmic ray
        scaled by the given factor.
        
        :Parameters:
        
        scale: float, optional, default=1.0
            An optional scale factor to be applied.

        :Returns:

        electrons: int
            The total number of electrons released by the cosmic ray.
            
        """
        # hit_map could be an array or a single value.
        if self.hit_map is not None:
            return (self.hit_map.sum() * float(scale))
        else:
            return 0
    
    def get_nucleon(self):
        """
        
        Return the nucleon involved.

        :Returns:

        nucleon: string
            A description of the nucleon responsible for this cosmic ray
            (which can be a null string).
        
        """
        return self.nucleon
    
    def get_hit_map(self, scale=1.0):
        """
        
        Returns the hit map scaled by the given factor.
        
        :Parameters:
        
        scale: float, optional, default=1.0
            An optional scale factor to be applied.

        :Returns:

        hit_map: array_like
            A map showing the distribution of cosmic ray energy between
            the central pixel and its neighbours in electrons.
            
        """
        if self.hit_map is not None:
            return self.hit_map * float(scale)
        else:
            return None

    def __str__(self):
        """
        
        Return a string describing an CosmicRay object.
        
        """
        strg = "Cosmic ray "
        if self.nucleon:
            strg += "of %s particle " % self.nucleon
        strg += "with energy %f MeV "  % self.energy
        
        strg += "hitting detector at (%d,%d)" % self.target_coords
        pix_affected = np.where(self.hit_map > 0.0)
        strg += "\n   and releasing %d electrons over %d pixels." % \
            (self.get_electrons(), len(self.hit_map[pix_affected]))
        if self._verbose > 4 and len(self.hit_map[pix_affected]) > 1:
            strg += "\n   Electron values: " + \
                self.hit_map[pix_affected].__str__()
        return strg

    def plot(self, plotfig=None, plotaxis=None, labels=True, withbar=False,
             title=None, description=''):
        """
        
        Plot the hit map associated with a cosmic ray event within
        the given matplotlib axis.
        This function can can be used to include a plot of this object
        in any figure.
        The plotfig method can be used to create a self-contained
        plot.

        :Parameters:
         
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.    
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        withbar: bool, optional
            Set to True to add a colour bar. The default is False.
        title: string, optional
            Optional title for the plot. The default is a string
            describing the CosmicRay object.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
        description: string, optional
            Optional description to be added to the title, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
  
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        # If the hit map is a 2-D array, plot the cosmic ray image.
        if len(self.hit_map.shape) == 2:
            # The hit map locations are with respect to the target pixel.
            irl = self.hit_map.shape[0]//2
            icl = self.hit_map.shape[1]//2
            extent = (-icl,icl,-irl,irl)
        
            if title is None:
                title = "Electrons from cosmic ray "
                if self.nucleon:
                    title += "(%s particle) " % self.nucleon
                    title += "with energy %.6g MeV"  % self.energy
            if description:
                title += "\n%s" % description
        
            # Display the cosmic ray map as an image.
            if labels:
                xlabel = 'Relative columns'
                ylabel = 'Relative rows'
            else:
                xlabel = ''
                ylabel = ''
            mplt.plot_image(self.hit_map, plotfig=plotfig, plotaxis=plotaxis,
                            xlabel=xlabel, ylabel=ylabel,
                            title=title, extent=extent)
        else:
            # TODO: What do we do with non-2D hit maps? Ignore them for now.
            pass
    

class CosmicRayEnvironment(object):
    """
    
    Class CosmicRayEnvironment - Simulates the behaviour of the cosmic ray
    radiation in the JWST environment. Galactic cosmic rays and solar
    particles are both included without being distinguished.
                 
    :Parameters:
    
    cr_flux: float
        The expected cosmic ray flux incident on the MIRI SCA in
        number of events per square micron per second.
    method: string
        The method of cosmic ray generation.
        
        * RANDOM - Generate random cosmic rays based on a probability
          distribution.
        * LIBRARY - Read cosmic ray descriptions from a prepared
          library of events
              
    energies: float array
        A vector of cosmic ray energies in MeV. This vector can be used
        in one of two ways, depending on the generation method:
        
        * RANDOM - A list of likely cosmic ray energies
        * LIBRARY - A list of cosmic ray energies read from the library.
             
    distribution: float array, optional
        A vector giving the cumulative probability distribution
        of the cosmic ray energies. Used by the RANDOM method only.
        This array must have the same number of elements as the energies
        array, and it must contain at least one non-zero value.
    images: array_like, optional
        A data cube containing a series of images of electrons released for
        each of the cosmic ray events contained in the energies list. Used by
        the LIBRARY method only.
        If not None, this array must have a slice for each element of the
        energies array.
    nucleons: int array, optional
        A list of nucleon IDs for each of the cosmic ray events contained
        in the energies list. Used by the LIBRARY method only.
        If not None, this array must be the same number of elements as the
        energies array.
    convolve_ipc: boolean (optional)
        If True, convolve all cosmic ray events with the detector
        capacitative coupling function. The default is True.
        Only relevant for the LIBRARY method.
    blur_sigma: float (optional).
        If greater than zero, apply a Gaussian blur to all cosmic ray
        events with this sigma in pixels. Defaults to zero.
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
    TypeError
        Raised if any of the initialisation parameters are of the
        wrong type, size or shape.
        
    """
    
    def __init__(self, cr_flux, method, energies, distribution=None,
                 images=None, nucleons=None, convolve_ipc=True,
                 blur_sigma=0.0, verbose=2, logger=LOGGER):
        """
        
        Constructor for class CosmicRayEnvironment.
        
        Parameters: See class doc string.
        
        """
        self.toplogger = logger
        self.logger = logger.getChild("scasim.cosmic_ray")
        if verbose > 4:
            self.logger.setLevel(logging.DEBUG)
            self.logger.debug( "+++CosmicRayEnvironment object created " + \
                "using method %s and flux %g (micron^2/sec)" % \
                (method, cr_flux))
        
        try:
            self.cr_flux = float(cr_flux)
        except (ValueError, TypeError) as e:
            strg = "Cosmic ray flux must be a valid number."
            strg += "\n %s" % e
            raise ValueError(strg)
        self.method = method
        self.energies = np.asarray(energies)
        
        if method == 'RANDOM':
            if distribution is not None:
                if len(distribution) == len(energies):
                    self.distribution = np.asarray(distribution)
                else:
                    strg = "distribution array (%d elements) " % \
                        len(distribution)
                    strg += "must be same size as energies array "
                    strg += "(%d elements)." % len(energies)
                    raise TypeError(strg)
                
                if self.distribution.max() <= 0.0:
                    strg = "The probability distribution array cannot "
                    strg += "be all zero or negative."
                    raise ValueError(strg)
            else:
                strg = "When method=RANDOM a distribution array " \
                    "must be provided"
                raise TypeError(strg)
            
            self.images = None
            self.nucleons = None
            self.coupling = None
        else:
            # Assume LIBRARY.
            self.distribution = None
            if images is not None:
                images = np.asarray(images)
                if images.ndim != 3:
                    strg = "Images array should be a data cube"
                    raise TypeError(strg)
                if images.shape[0] != len(energies):
                    strg = "Images cube (%d x %d x %d) " % images.shape
                    strg += "must have a slice for each element in the "
                    strg += "energies array (%d elements)." % len(energies)
                    raise TypeError(strg)
                self.images = np.asarray(images)
            else:
                self.images = None
            
            if nucleons is not None:
                if len(nucleons) != len(energies):
                    strg = "Nucleons array (%d elements) " % len(nucleons)
                    strg += "must be same size as energies array "
                    strg += "(%d elements)." % len(energies)
                    raise TypeError(strg)
                self.nucleons = np.asarray(nucleons)
            else:
                self.nucleons = None

            if convolve_ipc:
                self.coupling = np.array(cosmic_ray_properties['CR_COUPLING'])
            else:
                self.coupling = None
    
        self.blur_sigma = blur_sigma
        self._verbose = verbose

    def set_metadata(self, metadata):
        """
        
        Add cosmic ray keywords to the given MIRI metadata object.
        
        :Parameters:
        
        metadata: dictionary-like object
            A keyword-addressable metadata object to which cosmic
            ray environment keywords should be added.
            
        :Returns:
        
        metadata: dictionary-like object
            An updated metadata object.
            
        """
        metadata["GCR_FLUX"] = self.cr_flux
        metadata["CRMETHOD"] = self.method
        return metadata

    def set_seed(self, seedvalue=None):
        """
        
        Set the seed for the numpy random number generator.
        This function can be used while testing to ensure the set of
        RANDOM cosmic ray events that follows is well defined.
        
        :Parameters:
        
        seedvalue: int, optional, default=None
            The seed to be sent to the np.random number generator.
            If not specified, a value of None will be sent, which
            randomises the seed.
            
        """
        np.random.seed(seedvalue)

    def __str__(self):
        """
        
        Return a string describing a CosmicRayEnvironment object.
        
        """
        strg = "Cosmic ray environment modelled by %s method"  % self.method
        strg += "\n   with expected flux of %g events/micron^2/s" % \
            self.cr_flux
        strg += "\n   and energies ranging from %.1f to %.1f MeV." % \
            (self.energies.min(), self.energies.max())
        if self.images is not None:
            strg += "\n   and events ranging from %.1f to %.1f electrons." % \
                (self.images.min(), self.images.max())
        if self.coupling is not None:
            strg += "\n   An IPC coupling function is applied. "
        if self.blur_sigma > 0.0:
            strg += "\n   Events are blurred by a sigma of %.2f pixels." % \
                self.blur_sigma
        if cosmic_ray_properties['CR_BINNING_FACTOR'] > 1:
            strg += " Event pixels were binned by a factor of %d." % \
                cosmic_ray_properties['CR_BINNING_FACTOR']
        strg += "\n   There are %d simulated events." % len( self.energies )
        return strg

    def stats(self):
        """

        Return a string giving statistics for a CosmicRayEnvironment object.

        """
        strg = "Cosmic ray environment statistics\n"
        strg += "Energy array (MeV) min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
            (self.energies.min(), self.energies.max(), self.energies.mean(),
             self.energies.std(), np.median(self.energies))
        if self.images is not None:
            strg += "Images array (e) min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
                (self.images.min(), self.images.max(), self.images.mean(),
                 self.images.std(), np.median(self.images))
        
            # Generate an electron count (ignoring zero counts, which will
            # not display logarithmically).
            electrons_per_event = []
            electrons_per_pixel = []
            primary_pixels_per_event = []
            all_pixels_per_event = []
            for sl in range(0,self.images.shape[0]):
                crmap = self.images[sl,:,:]
    
                # If required, convolve the hit map with the capacitive
                # coupling and scattering functions of the detector.
                crmap = self._coupling_and_scattering(crmap)
                
                nonzero = np.where(crmap > 0)
                # Count the number of pixels affected by this event.
                if len(nonzero[0]) > 0:
                    all_pixels_per_event.append( len(nonzero[0]) )
                    if len(nonzero[0]) > 4:
                        maxe = crmap[nonzero].max()
                        mede = np.median(crmap[nonzero])
                        thresh = mede + (maxe-mede)/5.0
                        primary = np.where(crmap > thresh)
                        primary_pixels_per_event.append( len(primary[0]) )
                    else:
                        # Only one primary pixel
                        primary_pixels_per_event.append( 1 )
                    esum = 0
                    for element in crmap[nonzero]:
                        esum += element
                        electrons_per_pixel.append(element)
                    electrons_per_event.append(esum)
                
            # Calculate statistics on the cosmic ray events
            electrons_per_event = np.asarray(electrons_per_event)
            electrons_per_pixel = np.asarray(electrons_per_pixel)
            primary_pixels_per_event = np.asarray(primary_pixels_per_event)
            all_pixels_per_event = np.asarray(all_pixels_per_event)

            strg += "Electrons per event min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
                (electrons_per_event.min(), electrons_per_event.max(),
                 electrons_per_event.mean(),
                 electrons_per_event.std(), np.median(electrons_per_event))
            strg += "Electrons per pixel min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
                (electrons_per_pixel.min(), electrons_per_pixel.max(),
                 electrons_per_pixel.mean(),
                 electrons_per_pixel.std(), np.median(electrons_per_pixel))
            strg += "Primary pixels per event min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
                (primary_pixels_per_event.min(), primary_pixels_per_event.max(),
                 primary_pixels_per_event.mean(),
                 primary_pixels_per_event.std(), np.median(primary_pixels_per_event))
            strg += "All pixels per event min=%g; max=%g; mean=%g; std=%g; median=%g\n" % \
                (all_pixels_per_event.min(), all_pixels_per_event.max(),
                 all_pixels_per_event.mean(),
                 all_pixels_per_event.std(), np.median(all_pixels_per_event))
                        
        return strg
    
    def _coupling_and_scattering(self, input):
        """
        
        Helper function which applies IPC coupling and scattering to
        a 2-D array and returns a new 2-D array. Entries less than
        1.0 are replaced by 1.0.
        
        """
        output = input
        if self.coupling is not None:
            output = convolve2d(output, self.coupling)
        # If required, blur the hit map
        if self.blur_sigma > 0.0:
            output = gaussian_filter(output, self.blur_sigma)
        # Round to the nearest whole number of electrons.
        output = output.round(0)
        return output
    
    def plot_energies(self, plotfig, plotaxis, labels=True, xscale='log',
                      yscale='linear', title=''):
        """

        Add a plot of the distribution of cosmic ray energies
        to a given matplotlib axis. This function can can be used
        to include a plot of this object in any figure created by
        the caller.
        
        :Parameters:
        
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.  
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        xscale: string, optional
            X axis scaling type: 'linear' or 'log'.
            The default is 'log'.
        yscale: string, optional
            Y axis scaling type: 'linear' or 'log'.
            The default is 'linear'.
        title: string, optional
            Optional title to be shown above the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
          
        :Returns:
        
        plotaxis: matplotlib axis object
            An updated axis containing the plot.
              
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot

        """
        # Plot the energy distribution with a log scale on the X axis.
        if labels:
            xlabel = 'Energy per event (MeV)'
            ylabel = 'Relative Frequency'
        else:
            xlabel = ''
            ylabel = ''
        plotaxis = mplt.plot_xy(self.energies, self.distribution,
                                plotfig=plotfig, plotaxis=plotaxis,
                                xscale=xscale, xlabel=xlabel,
                                ylabel=ylabel, yscale=yscale, title=title)
        return plotaxis

    def hist_energies(self, plotfig, plotaxis, nbins=256, labels=True,
                      xscale='log', yscale='log', title='', ):
        """

        Add a histogram of the distribution of cosmic ray energies
        to a given matplotlib axis. This function can can be used
        to include a plot of this object in any figure created by
        the caller.
        
        :Parameters:
        
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
        nbins: int, optional, default=256
            The number of histogram bins to be plotted        
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        xscale: string, optional
            X axis scaling type: 'linear' or 'log'.
            The default is 'log'.
        yscale: string, optional
            Y axis scaling type: 'linear' or 'log'.
            The default is 'log'.
        title: string, optional
            Optional title to be shown above the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
          
        :Returns:
        
        plotaxis: matplotlib axis object
            An updated axis containing the plot.

        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot

        """
        # Plot the energies with a log scale on the X axis, displaying
        # the bins with equal width.
        if labels:
            xlabel = 'Energy per event (MeV)'
            ylabel = 'Number of events'
        else:
            xlabel = ''
            ylabel = ''
        plotaxis = mplt.plot_hist(self.energies, nbins, equalwidths=True,
                                  plotfig=plotfig, plotaxis=plotaxis,
                                  xscale=xscale, xlabel=xlabel,
                                  yscale=yscale, ylabel=ylabel, title=title)
        return plotaxis

    def hist_electrons(self, plotfig, plotaxis, nbins=256, summed=True,
                       labels=True, xscale='log', yscale='log', title=''):
        """

        Add a histogram of the distribution of released detector electrons
        to a given matplotlib axis. This function can can be used
        to include a plot of this object in any figure created by
        the caller.
        
        :Parameters:
        
        plotfig: matplotlib figure object
            Figure on which to add the plot.
        plotaxis: matplotlib axis object
            Axis on which to add the plot.
        nbins: int, optional, default=256
            The number of histogram bins to be plotted
        summed: bool, optional
            Set to True (the default) if the histogram should show the
            total number of electrons released per event.
            Set to False if the histogram should show the number of
            electrons released separately per pixel.
        labels: bool, optional
            Set to False to suppress axis labels. The default is True.
        xscale: string, optional
            X axis scaling type: 'linear' or 'log'.
            The default is 'log'.
        yscale: string, optional
            Y axis scaling type: 'linear' or 'log'.
            The default is 'log'.
        title: string, optional
            Optional title to be shown above the plot, if required.
            Note: If too much text is written on a plot it may
            overlap with other labels; especially when applied to
            subplots.
          
        :Returns:
        
        plotaxis: matplotlib axis object
            An updated axis containing the plot.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot

        """
        # Generate an electron count (ignoring zero counts, which will
        # not display logarithmically).
        electrons = []
        for sl in range(0,self.images.shape[0]):
            crmap = self.images[sl,:,:]

            # If required, convolve the hit map with the capacitive
            # coupling and scattering functions of the detector.
            crmap = self._coupling_and_scattering(crmap)
            
            nonzero = np.where(crmap > 0)
            esum = 0
            for element in crmap[nonzero]:
                if summed:
                    esum += element
                else:
                    electrons.append(element)
            if summed and esum > 0:
                electrons.append(esum)

        # Plot the electron events with a log scale on the X axis, displaying
        # the bins with equal width.
        if labels:
            if summed:
                xlabel = 'Electrons generated per event'
            else:
                xlabel = 'Electrons generated per pixel'
            if summed:
                ylabel = 'Number of events'
            else:
                ylabel = 'Number of pixels'
        else:
            xlabel = ''
            ylabel = ''
        plotaxis = mplt.plot_hist(electrons, nbins, equalwidths=True,
                                  plotfig=plotfig, plotaxis=plotaxis,
                                  xscale=xscale, xlabel=xlabel,
                                  yscale=yscale, ylabel=ylabel, title=title)
        return plotaxis

    def plot(self, nbins=256, description=''):
        """
        
        Plot the distribution of cosmic ray energies in a
        self-contained matplotlib figure and display it.
        
        :Parameters:
        
        nbins: int, optional, default=256
            The number of histogram bins to be plotted
        description: string, optional
            Additional description to be shown on the plot, if required.
            
        :Global variables:
        
        plt: matplotlib pyplot object
            The top level pyplot object, imported by this module.
            
        :Requires:
        
        miri.tools.miriplot
        matplotlib.pyplot
            
        """
        if self.method == 'LIBRARY':
            # If a cosmic ray library is being used, plot the
            # frequency of energies contained in the library.

            # Create a figure to contain the plot.
            tstrg = "Cosmic ray energy distribution from LIBRARY"
            if description:
                tstrg += " : %s" % description   
            fig = mplt.new_figure(1, figsize=(21,7), stitle=tstrg)

            # Create a subplot for the energy distribution.
            ax1 = mplt.add_subplot(fig, 1, 3, 1)
            ax1 = self.hist_energies(fig, ax1, nbins,
                        title="Cosmic ray energy distribution")

            # Create a subplot for the corresponding distribution of 
            # detector events.
            ax2 = mplt.add_subplot(fig, 1, 3, 2)
            ax2 = self.hist_electrons(fig, ax2, nbins, summed=True,
                        title="Detector event distribution (summed by event)")

            # Create a subplot for the corresponding distribution of 
            # detector events.
            ax3 = mplt.add_subplot(fig, 1, 3, 3)
            ax3 = self.hist_electrons(fig, ax3, nbins, summed=False,
                        title="Detector event distribution (for each pixel)")
                
            # Close and display the plot.
            if self._verbose > 1:
                strg =  "Plotting cosmic ray distributions. "
                strg += "Close the plot window to continue."
            else:
                strg = ''
            mplt.show_plot(prompt=strg)

        else:
            # If cosmic rays are being generated randomly
            # from a pre-defined distribution, plot that
            # distribution.                
            
            # Create a figure to contain the plot.
            tstrg = "Cosmic ray energy distribution (selected at RANDOM)"
            if description:
                tstrg += "\n%s" % description   
            fig = mplt.new_figure(1, figsize=(10,8), stitle=tstrg)

            # Plot the energy distribution, filling the whole figure.
            ax1 = mplt.add_subplot(fig, 1, 1, 1)
            ax1 = self.plot_energies(fig, ax1)

            # Close and display the plot.
            if self._verbose > 1:
                strg =  "Plotting cosmic ray distributions. "
                strg += "Close the plot window to continue."
            mplt.show_plot(prompt=strg)
    
    def generate_events(self, rows, columns, time, pixsize):
        """
        
        Generate a list of cosmic ray events which happen while a
        specified detector is integrating for a specified time.
        
        :Parameters:
        
        rows: int
            The number of detector rows.
        columns: int
            The number of detector columns.
        time: float
            The integration time in seconds.
        pixsize: float
            The detector pixel size in microns.
            
        :Returns:
        
        event_list: list of CosmicRay objects
            A list of CosmicRay objects describing the events which
            take place.
            
        """
        # Cosmic ray events will be collected into an event list.
        event_list = []

        # Generate detector events if there are detector pixels.
        if rows > 0 and columns > 0:
            # First calculate the number of detector events that are expected,
            # which is determined by the cosmic ray flux, the total detector
            # area and the integration time.
            expected_events = int(self.cr_flux * rows * columns * \
                pixsize * pixsize * time)
            
            # The actual number of events will be a number randomly
            # sampled from a Poisson distribution.
            nevents = np.random.poisson(expected_events)
            if (self._verbose > 2) and (expected_events > 0):
                self.logger.info( "Number of detector CR events " + \
                    "expected=%d and actual=%d." % (expected_events, nevents))
                
            # For each of these events choose a random target pixel and
            # select a cosmic ray energy either from the given probability
            # distribution or from the library.
            for ev in range(0,nevents):
                cosmic_ray = self.generate_event(rows, columns)
                event_list.append(cosmic_ray)
            
        return event_list
    
    def generate_event(self, rows, columns):
        """
        
        Generate a single cosmic ray event which hits the
        detector assembly at a random coordinate.
        
        :Parameters:
        
        rows: int
            The number of detector rows.
        columns: int
            The number of detector columns.
            
        :Raises:
    
        ValueError
            Raised if any of the parameters are out of range.

        :Returns:
        
        event: CosmicRay
            A cosmic ray event
            
        """
        # The cosmic ray can hit the detector anywhere over its surface
        hit_row = int(np.random.uniform(0, rows))
        hit_column = int(np.random.uniform(0, columns))
        hit_coords = (hit_row, hit_column)
            
        if self.method == 'RANDOM':
            # Select a random event from the list of possible energies
            # weighted by the given incremental probability
            # distribution.
            prob = np.random.uniform(0.0, self.distribution.sum())
            lkup = -1
            ptotal = 0.0
            while (ptotal < prob) and (lkup < len(self.energies)):
                lkup += 1
                ptotal += self.distribution[lkup]
            energy = self.energies[lkup]
            # The hit map is generated from the capacitative coupling
            hit_map = np.array(cosmic_ray_properties['CR_COUPLING']) * energy

            # If required, blur the map to simulate scattering.
            hit_map = self._coupling_and_scattering(hit_map)
            nucleon=''
        else:
            # Select a random event from the list with a uniform
            # probability distribution.
            lkup = int(np.random.uniform(0, len(self.energies)))
            energy = self.energies[lkup]
            hit_map = self.images[lkup]

            # If required, convolve the hit map with the capacitive
            # coupling and scattering functions of the detector.
            hit_map = self._coupling_and_scattering(hit_map)

            nc = self.nucleons[lkup]
            nucleon = cosmic_ray_properties.get('CR_NUCLEONS')[nc]

        cosmic_ray = CosmicRay(energy, hit_coords, hit_map=hit_map,
                               nucleon=nucleon, verbose=self._verbose,
                               logger=self.logger)
        
        return cosmic_ray        


def load_cosmic_ray_random(cosmic_ray_mode='SOLAR_MIN', verbose=2,
                           logger=LOGGER):
    """
    
    Create a CosmicRayEnvironment object based on the estimated energy
    distributions contained in cosmic_ray_properties configuration file
    This will not be as realistic as an environment based on
    the STScI cosmic ray library.
    
    :Parameters:
    
    cosmic_ray_mode: string, optional, default='SOLAR_MIN'
        The cosmic ray mode to be simulated. Available modes are:
        
        * 'NONE' - No cosmic rays.
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)
        
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
     
    :Returns:
    
    A new CosmicRayEnvironment object.
    
    """
    if (verbose > 1) and (cosmic_ray_mode != 'NONE'):
        logger.info( "Creating randomised cosmic ray environment for mode %s" % \
            cosmic_ray_mode)
    
    cr_flux = cosmic_ray_properties.get('CR_FLUX', cosmic_ray_mode) * \
        cosmic_ray_properties['CR_FLUX_MULTIPLIER']
        
    cr_e_dist = cosmic_ray_properties.get('CR_ELECTRONS', cosmic_ray_mode)
    cr_e_trans = np.transpose(cr_e_dist)
    electrons = cr_e_trans[0] * \
                cosmic_ray_properties['CR_ENERGY_MULTIPLIER']
    distribution = cr_e_trans[1]
    
    cr_env = CosmicRayEnvironment(cr_flux, 'RANDOM', energies=electrons,
                                  distribution=distribution,
                                  images=None, nucleons=None,
                                  verbose=verbose, logger=logger)
    return cr_env

def load_cosmic_ray_library(filename, cosmic_ray_mode='SOLAR_MIN',
                            convolve_ipc=True, verbose=2, logger=LOGGER):
    """
    
    Create a CosmicRayEnvironment object from the information contained
    in a STScI-format cosmic ray library file.
    
    :Parameters:
    
    filename: string
        Name of file containing cosmic ray library.
    cosmic_ray_mode: string, optional, default='SOLAR_MIN'
        The cosmic ray mode to be simulated. Available modes are:
        
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)

    convolve_ipc: boolean (optional)
        If True, convolve all cosmic ray events with the detector
        capacitative coupling function. The default is True.
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
     
    :Returns:
    
    A new CosmicRayEnvironment object.
    
    """
    cr_flux = cosmic_ray_properties.get('CR_FLUX', cosmic_ray_mode) * \
            cosmic_ray_properties['CR_FLUX_MULTIPLIER']

    if verbose > 1:
        logger.info( "Reading cosmic ray library file: \'%s\'" % filename)
    hdulist = pyfits.open(filename)
    images = hdulist[1].data * \
                cosmic_ray_properties['CR_ELECTRON_MULTIPLIER']
    nucleons = hdulist[2].data
    energies = hdulist[3].data * \
                cosmic_ray_properties['CR_ENERGY_MULTIPLIER'] 
    sigma = cosmic_ray_properties['CR_BLUR']

    # If necessary, rebin the images array down to a new size.
    binfactor = cosmic_ray_properties['CR_BINNING_FACTOR']
    if  binfactor > 1:
        newimages = []
        newshape = (images.shape[-2]//binfactor, images.shape[-1]//binfactor)
        maxshape = (newshape[0] * binfactor, newshape[1] * binfactor)
        for dim3 in range(0,images.shape[-3]):
            slice = images[dim3,:maxshape[0],:maxshape[1]]
            newslice = _rebin( slice, newshape )
            newimages.append(newslice)
        images = np.asarray(newimages)
    
    cr_env = CosmicRayEnvironment(cr_flux, 'LIBRARY', energies,
                                  distribution=None,
                                  images=images, nucleons=nucleons,
                                  convolve_ipc=convolve_ipc,
                                  blur_sigma=sigma,
                                  verbose=verbose, logger=logger)
    return cr_env

def load_cosmic_ray_libraries(cosmic_ray_mode='SOLAR_MIN', variant='',
                              convolve_ipc=True, verbose=2, logger=LOGGER):
    """
    
    Create a CosmicRayEnvironment object from the information contained
    in all known STScI-format cosmic ray library files for the specified
    cosmic ray mode.
    
    Similar to load_cosmic_ray_library, except the simulations from more
    than one library file are combined into one cosmic ray environment.
    
    :Parameters:
    
    cosmic_ray_mode: string, optional, default='SOLAR_MIN'
        The cosmic ray mode to be simulated. Available modes are:
        
        * 'SOLAR_MIN' - Solar minimum
        * 'SOLAR_MAX' - Solar maximum
        * 'SOLAR_FLARE' - Solar flare (worst case scenario)
        
    variant: string, optional, default=''
        The variant of the cosmic ray libraries to be read (if any).
        The default of '' reads the standard version.
    convolve_ipc: boolean (optional)
        If True, convolve all cosmic ray events with the detector
        capacitative coupling function. The default is True.
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
     
    :Returns:
    
    A new CosmicRayEnvironment object.
    
    """
    # First obtain a list of available files for the given cosmic ray mode.
    file_list = []
    num_min = cosmic_ray_properties.get('CR_LIBRARY_FILES', 'MIN')
    num_max = cosmic_ray_properties.get('CR_LIBRARY_FILES', 'MAX')
    for filenum in range(num_min, num_max+1):
        if variant:
            file_name = '%s0%d_%s.fits' % \
                (cosmic_ray_properties.get('CR_LIBRARY_FILES', cosmic_ray_mode), \
                 filenum, variant)
        else:
            file_name = '%s0%d.fits' % \
                (cosmic_ray_properties.get('CR_LIBRARY_FILES', cosmic_ray_mode), \
                 filenum)
        file_list.append(file_name)
    
    # The flux level depends on cosmic ray mode.
    cr_flux = cosmic_ray_properties.get('CR_FLUX', cosmic_ray_mode) * \
            cosmic_ray_properties['CR_FLUX_MULTIPLIER']

    # Read the files one by one and concatenate the contents.
    firsttime = True
    for filename in file_list:
        if verbose > 1:
            logger.info( "Reading cosmic ray library file: \'%s\'" % filename)
        hdulist = pyfits.open(filename)
        if firsttime:
            images = hdulist[1].data
            nucleons = hdulist[2].data
            energies = hdulist[3].data
            firsttime = False
        else:
            images_new = hdulist[1].data
            nucleons_new = hdulist[2].data
            energies_new = hdulist[3].data
            images = np.concatenate( (images, images_new) )
            nucleons = np.concatenate( (nucleons, nucleons_new) )
            energies = np.concatenate( (energies, energies_new) )
    images = images * cosmic_ray_properties['CR_ELECTRON_MULTIPLIER'] 
    energies = energies * cosmic_ray_properties['CR_ENERGY_MULTIPLIER'] 
    sigma = cosmic_ray_properties['CR_BLUR']

    # If necessary, rebin the images array down to a new size.
    binfactor = cosmic_ray_properties['CR_BINNING_FACTOR']
    if  binfactor > 1:
        newimages = []
        newshape = (images.shape[-2]//binfactor, images.shape[-1]//binfactor)
        maxshape = (newshape[0] * binfactor, newshape[1] * binfactor)
        for dim3 in range(0,images.shape[-3]):
            slice = images[dim3,:maxshape[0],:maxshape[1]]
            newslice = _rebin( slice, newshape )
            newimages.append(newslice)
        images = np.asarray(newimages)

    # Create a grand cosmic ray environment from the combined data.
    cr_env = CosmicRayEnvironment(cr_flux, 'LIBRARY', energies,
                                  distribution=None,
                                  images=images, nucleons=nucleons,
                                  convolve_ipc=convolve_ipc,
                                  blur_sigma=sigma,
                                  verbose=verbose, logger=logger)
    return cr_env

def plot_cosmic_ray_events(event_list, plotrows, plotcols):
    """
    
    Plot the hit maps from the first few events contain in a
    cosmic ray event list.
    
    :Parameters:
    
    event_list: list of CosmicRay objects
        A list of cosmic ray events. The first (plotrows*plotcols)
        events will be plotted.
    plotrows: int
        The number of subplot rows to be used.
    plotcols: int
        The number of subplot columns to be used.
            
    :Global variables:
        
    plt: matplotlib pyplot object
        The top level pyplot object, imported by this module.
            
    :Requires:
        
    matplotlib.pyplot
    
    :Raises:
    
    ImportError
        Raised if the matplotlib plotting library is not available.
    
    """
    # The number of events to be plotted is governed either by the
    # plotting space available or the number of events available,
    # whichever is the smaller.
    nevents = int(plotrows) * int(plotcols)
    if nevents > len(event_list):
        nevents = len(event_list)
        
    # Create a figure to contain the plot (larger than default size).
    tstrg = "Electrons generated by the first %d CR events" % nevents
    fig = mplt.new_figure(1, figsize=(10,10), stitle=tstrg)
    
    # Step through the list of and plot the first nevents.
    count = 0
    for cr in event_list:
        count += 1
        if count <= nevents:
            title = "Event %d" % count
            ax = mplt.add_subplot(fig, plotrows, plotcols, count)
            ax = cr.plot(plotfig=fig, plotaxis=ax, labels=False, title=title)
        else:
            break
    
    strg = "Plotting %s. Close the plot window to continue" % tstrg
    mplt.show_plot(prompt=strg)

#
# The following code is for development and testing only. It will run
# a few ad-hoc tests exercising the code. A more formal set of unit tests
# may be found in the scasim/tests/test_cosmic_ray module.
# However, the tests here are made in verbose mode and include more file
# I/O and plotting than the unit tests, so they show what is happening
# in more detail.
#
if __name__ == '__main__':
    print( "Testing the CosmicRayEnvironment package" )
    PLOTTING = True        # Set to False to turn off plotting.
    VERBOSE = 3
    
    TEST_RANDOM = False
    TEST_LIBRARIES = True
    TEST_LIBRARY = False
    TEST_EVENTS = True
    
    if TEST_RANDOM:
#         energies = (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000)
#         distribution = (0.05, 0.07, 0.13, 0.15, 0.16, 0.13, 0.1, 0.09, 0.07, 0.05)
#         cr_env0 = CosmicRayEnvironment(4.0e-8, 'RANDOM', energies, distribution,
#                                        verbose=VERBOSE)
        
        for env in ('SOLAR_MIN', 'SOLAR_MAX', 'SOLAR_FLARE'):
            print( "\nTesting a RANDOM distribution at %s" % env )
            cr_env1 = load_cosmic_ray_random(env, verbose=VERBOSE)
            print( cr_env1 )
            print( cr_env1.stats() )
            if PLOTTING:
                strg = "(test plot at %s)" % env
                cr_env1.plot( nbins=32, description=strg )

    if TEST_LIBRARIES:
        for crmode in ('SOLAR_MIN', 'SOLAR_MAX', 'SOLAR_FLARE'):
            for convolve_ipc in (False, True):
                print( "\nTesting combined LIBRARY " + \
                       "for %s (convolve_ipc=%s)" % \
                       (crmode, str(convolve_ipc)) )
                cr_env2 = load_cosmic_ray_libraries( crmode, variant='',
                                                     convolve_ipc=convolve_ipc,
                                                     verbose=VERBOSE )
                print( cr_env2 )
                print( cr_env2.stats() )
                if PLOTTING:
                    strg = "(multi-file test plot at %s" % crmode
                    if convolve_ipc:
                        strg += ", convolved with IPC"
                    strg += ")"
                    cr_env2.plot( nbins=1024, description=strg )
        del cr_env2

    if TEST_LIBRARY:
        for crmode in ('SOLAR_MIN', 'SOLAR_MAX', 'SOLAR_FLARE'):
            filenum = int(np.random.uniform(0,9))
            print( "\nTesting single LIBRARY for %s (number %d)" % (crmode, filenum) )
            test_file_name = '%s0%d.fits' % \
                (cosmic_ray_properties.get('CR_LIBRARY_FILES', crmode), \
                 filenum)
            try:
                cr_env2 = load_cosmic_ray_library(test_file_name, crmode,
                                                 convolve_ipc=True,
                                                 verbose=VERBOSE)
            except IOError:
                print( "***Could not find test file %s. Using RANDOM mode." % \
                       test_file_name )
                cr_env2 = load_cosmic_ray_random(crmode, verbose=VERBOSE)
     
            print( cr_env2 )
            print( cr_env2.stats() )
            if PLOTTING:
                strg = "(single file test plot at %s)" % crmode
                cr_env2.plot( nbins=256, description=strg )
        del cr_env2
    
    if TEST_EVENTS:
        crmode = 'SOLAR_MIN'
        filenum = int(np.random.uniform(0,9))
        print( "\nGenerating events for %s (number %d)" % ('SOLAR_MIN', filenum) )
        test_file_name = '%s0%d.fits' % \
            (cosmic_ray_properties.get('CR_LIBRARY_FILES', crmode), \
             filenum)
        cr_env3 = load_cosmic_ray_library(test_file_name, crmode,
                                          convolve_ipc=True,
                                          verbose=VERBOSE)
        
        event_list = cr_env3.generate_events(1024, 1024, 2.75, 25.0)
        for cr in event_list:
            print( cr )
        
        if PLOTTING:
            cr.plot(description='(last event)')
            plot_cosmic_ray_events( event_list, 5, 5)

        del event_list, cr_env3
    print( "Test finished." )
