#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""

Module miriplot - Contains general purpose plotting functions for
MIRI data.

The module is based on matplotlib, and therefore recognises the following
colours:

b = blue
g = green
r = red
c = cyan
m = magenta
y = yellow
k = black
w = white
0.5 = 50% gray
#eeefff = hex colour
(r,g,b) = RGB tuple

The following linestyles are recognised:

'-' = solid line
'--' dashed line
'-.' = dot dash line
':' = dotted line

The following markers are recognised:

. = point
, = pixel
o = circle
v = down triangle
^ = up triangle
< = left triangle
> = right triangle
s = square
p = pentagon
* = star
h = hexagon1
H = hexagon2
x = x
D = diamond
d = thin diamond 

:History:

23 Mar 2011: Created (for discussion)
24 Mar 2011: Improved with comments from Julien Morin.
11 Jan 2012: Renamed symbol to avoid potential name clash: slice-->slicenum
16 Mar 2012: Corrected Sphinx documentation typos.
02 Apr 2012: Better work around for matplotlib problem. Plotting disabled
             under Windows until the problem has been solved.
04 Apr 2012: Improvements suggested by pylint.
10 Apr 2012: matplotlib problem solved. It was caused by duplicate entries
             in the PYTHONPATH.
16 May 2012: Reworked for use by JWST/MIRI data products. "figure" and "axis"
             variants of the plotting functions merged into one function.
             Function names made easier. Some fixed parameters made
             configurable.
18 May 2012: plot_hist added, plus new-figure,add_subplot and show_plot.
22 May 2012: Clipping option added to plot_hist and plot_image.
15 Jun 2012: Removed barchar parameter from show_header().
28 Aug 2012: Added the option to draw error bars on the plots.
10 Sep 2012: showplot option added to most plotting functions.
             Setting showplot=False allows the caller to add extra
             annotation to a plot before showing it. 
11 Sep 2012: maxplots limit added to plot_xycolumn function.
19 Sep 2012: Logarithmic plotting option added to image display functions.
02 Oct 2012: Improved memory management. Copies of data not made unless
             essential.
23 Oct 2012: Removed superflous code.
13 Nov 2012: Major restructuring of the package folder. Import statements
             updated.
12 Jun 2014: Allow the Z axis label to be specified in 3-D plots.
25 Jun 2014: Close and remove figures after they have been displayed.
11 Jul 2014: Modified to use the new MiriFilter and MiriMeasurement classes
             instead of the old Filter and MeasuredVariable classes.
08 Sep 2015: Made compatible with Python 3
07 Oct 2016: Added equal_aspect parameter. Added plot_circles and plot_ellipses.
04 Jan 2017: Make the module resilient to matplotlib import problems.
             Corrected references to old miri.miritools packages.
18 May 2018: Changed deprecated logger.warn() to logger.warning(). Fixed
             incorrect log statements.

@author: Steven Beard (UKATC)

"""
# For consistency, import the same Python V3 features as the STScI data model.
from __future__ import absolute_import, unicode_literals, division, print_function

# Get a logger.
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("miri.tools.miriplot")

import math
import numpy as np
import numpy.ma as ma
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
except ImportError as e:
    plt = None
    patches = None
    _strg = "Failed to import matplotlib plotting library:"
    _strg += "\n\t%s: %s" % (e.__class__.__name__, e)
    logger.error(_strg)
#     raise ImportError(_strg)
except Exception as e:
    plt = None
    patches = None
    _strg = "Error while importing matplotlib plotting library:"
    _strg += "\n\t%s: %s" % (e.__class__.__name__, e)
    logger.error(_strg)
#     raise Exception(_strg)


#
# General purpose utility functions.
#
def test_available():
    """

    Test whether the MIRI plotting module is available for use.
    An import problem will render the module unusable.

    :Parameters:
        
    None.
        
    :Returns:
    
    available: boolean
        True if the module is ok and available.
        False if the module is not ok and cannot be used.
     
    """
    if plt is not None and patches is not None:
        return True
    else:
        return False
    
def subdivide(nplots, ratio=1.0, maxplots=None):
    """
        
    Subdivide a plotting space so that nplots can be displayed
    in a grid of subplots with the given aspect ratio.
        
    :Parameters:
        
    nplots: int
        The number of plots to be displayed. Must be at least 1.
    ratio: float, optional
        The required aspect ratio (columns/rows). Defaults to 1.0.
    maxplots: int, optional
        The maximum number of plots to be included on one figure.
        If nplots > maxplots the plots will require more than one
        figure.
        If maxplots is not provided, all the plots will be fitted
        onto one figure, no matter how many.
        
    :Returns:
    
    (subrows, subcols, nfigs): tuple of 3 ints
        The number of subplot columns and rows, plus the number of
        figures. nfigs will be 1 unless nplots > maxplots
        
    """
    nfigs = 1
    if float(ratio) <= 0.0:
        ratio = 1.0
    if maxplots is not None and nplots > maxplots:
        nfigs = int(np.ceil((nplots / float(maxplots))))
        nplots = maxplots
    if nplots > 1:
        subcols = int(np.sqrt(nplots * ratio))
        if subcols < 1:
            subcols = 1
        subrows = int(np.ceil((nplots / float(subcols))))
    else:
        subcols = 1
        subrows = 1
            
    return (subrows, subcols, nfigs)

def new_figure(page, figsize=(10,8), equal_aspect=False,
               stitle='', pyplt=plt, **kwargs):
    """
    
    Create a new figure. Matplotlib wrapper.
    This function allows a caller to create a figure and add subplots
    to it without needing to import the matplotlib library separately.

    :Parameters:
    
    page: int
        Page number for the figure.
    figsize: tuple of 2 ints, optional
        The size of the new figure. The default value is (10,8), which
        will create a 10 inch x 8 inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    stitle: str, optional
        An optional supertitle to be written at the top of the figure.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.figure    
    
    :Returns:
    
    plotfig: matplotlib figure object
        The figure object just created.
    
    """
    if pyplt is not None:
        fig = pyplt.figure(page, figsize=figsize, **kwargs)
        if equal_aspect:
            plt.gca().set_aspect('equal', adjustable='box')
        if stitle:
            fig.suptitle(stitle)
        return fig
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def add_subplot(plotfig, subrows=1, subcols=1, subno=1, **kwargs):
    """
    
    Add a new subplot to an existing figure. Matplotlib wrapper.
    This function allows a caller to create a figure and add subplots
    to it without needing to import the matplotlib library separately.

    :Parameters:
    
    plotfig: matplotlib figure object
        Figure on which to add the plot.
    subrows: int, optional
        The total number of subplot rows within the figure
    subcols: int, optional
        The total number of subplot columns within the figure
    subno: int, optional
        The running number for this particular subplot.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.add_subplot 
        
    :Returns:
       
    plotaxis: matplotlib axis object
        The axis object just created.
        
    """
    if plotfig is not None:
        return plotfig.add_subplot(subrows, subcols, subno, **kwargs)
    else:
#         logger.warning("matplotlib figure not available")
        return None

def show_plot(pyplt=plt, prompt=''):
    """
    
    Close and display the current plot. Matplotlib wrapper.
    This function allows a caller to finish and display a plot
    without needing to import the matplotlib library separately.

    :Parameters:
    
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    prompt: str, optional
        An optional string, which may be printed when the plot
        is displayed. (Program execution may halt until the plot
        window is closed.)
    
    """
    if pyplt is not None:
        if prompt:
            print( prompt )
        # Display the current plot and then clear and close the current figure.
        pyplt.show()
        pyplt.clf()
        pyplt.close()
    else:
        logger.warning("matplotlib.pyplot not available")

def close(pyplt=plt):
    """
    
    Close all open plot windows and free up resources.

    :Parameters:
    
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
        
    """
    if pyplt is not None:
        # Close all open figure windows
        pyplt.close('all')

#
# Plotting functions. Many of these functions can be used stand-alone
# without needing to set up the plotting environment with matplotlib.
# Leaving the plotaxis argument at its default setting of None will cause
# these functions to create a new figure for the plot. The functions
# can also be used within a matplotlib environment to add their plots
# to a predefined figure and axis, which can be provided through the
# plotfig and plotaxis arguments. The new_figure, add_subplot and
# show_plot wrapper functions (above) can be used to create simple
# plot layouts without needing full access to matplotlib.
#
def plot_text(text, xpos=0.05, ypos=0.95, pyplt=plt, plotfig=None,
              plotaxis=None, figsize=(10,6), equal_aspect=False, showplot=None,
              title='', horizontalalignment='left', verticalalignment='top',
              **kwargs):
    """

    Display text on a maplotlib plot. The function can add text
    to an existing matplotlib axis or create a new one.
        
    :Parameters:
    
    text: str
        The text to be added. matplotlib mathematical codes, such as
        "$\sigma_i=15" can be included to get special symbols.
    xpos: float, optional, default=0.05
        X coordinate of the start of the text, in viewport coordinates
        (0.0 to 1.0).
    ypos: float, optional, default=0.95
        Y coordinate of the start of the text, in viewport coordinates
        (0.0 to 1.0).
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    title: str, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    horizontalalignment: str, optional
        The required horizontal alignment ('left', 'right', 'center' or
        'baseline'). Defaults to 'left', which allows a paragraph of
        text to be added starting at the chosen location.
    verticalalignment: str, optional
        The required horizontal alignment ('bottom', 'top', 'center' or
        'baseline'). Defaults to 'top', which allows a paragraph of
        text to be added starting at the chosen location.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.plot
        For example fontsize, color, etc...
        See the matplotlib documentation for a full list.
          
    :Returns:
        
    plotaxis: matplotlib axis object
        The matplotlib axis containing the plot.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None:
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            # TODO: Can we remove the default tick marks created by this call?
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
    
        # Add the text and include a title if necessary.
        plotaxis.text(xpos, ypos, text, horizontalalignment=horizontalalignment,
                      verticalalignment=verticalalignment, **kwargs)
        if title:
            plotaxis.set_title(title)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            plotfig.clear()
            del plotfig, plotaxis
#             pyplt.clf()
#             pyplt.close()
            plotaxis = None
            
        # Return the axis object in which the text was added.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_xy(xdata, ydata, yerr=None, pyplt=plt, plotfig=None, plotaxis=None,
            figsize=(10,6), equal_aspect=False, showplot=None,
            xscale='linear', yscale='linear', linefmt='b-', linestyle='-',
            xlabel='', ylabel='', title='', grid=True, **kwargs):
    """

    Create an X,Y scatter plot. The function can add a scatter plot
    to an existing matplotlib axis or create a new one.
        
    :Parameters:
    
    xdata: array-like
        X axis data. This must either be a 1-D array, which is
        interpreted as a list of X coordinates, or a 2-D array, which
        is interpreted as multiple lists of X coordinates for a set
        of curves. It may be set to None, in which case the index
        numbers for the ydata array are used as X coordinates.
    ydata: array-like
        Y axis data. This must either be a 1-D array, which is
        interpreted as a list of Y coordinates, or a 2-D array, which
        is interpreted as multiple lists of Y coordinates for a set
        of curves.
    yerr: array-like, optional
        An array of Y axis errors used to draw error bars on the plot.
        If not specified, there will be no error bars.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    xscale: string, optional
        X axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    yscale: string, optional
        Y axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    linefmt: str, optional
        Matplotlib line format. Default 'b-' = blue line with no symbols.
    linestyle: str, optional.
        Matplotlib line style. Default '-' = solid line.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no label.
    ylabel: string, optional
        Optional Y axis label to be shown on the plot.
        The default is no label.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    grid: bool, optional
        If True, draw a grid on the plot. The default is True.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.plot
        For example label, linewidth, linestyle, color, marker, etc...
        See the matplotlib documentation for a full list.
          
    :Returns:
        
    plotaxis: matplotlib axis object
        The matplotlib axis containing the plot.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None:
        # Check whether matplotlib figure and axis objects have been provided.
        # If not then create new ones. By default, the plot will completely
        # fill the figure. Set a flag to remember that a new figure and/or axis
        # object have been created.
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
    
        # Ensure the Y data array is converted to a numpy ndarray object,
        # which has .ndim and .shape() method. The data can only be 1-D or 2-D.
        plotydata = np.asarray(ydata)
        if plotydata.ndim > 2:
            strg = "Data array has %d dimensions. 1 or 2 expected." % \
                plotydata.ndim
            raise TypeError(strg)
        
        # If an X axis array is not provided it defaults to the index
        # numbers of the first dimension of the Y array. Otherwise,
        # ensure the X data array is also converted to a numpy ndarray
        # object.
        if xdata is None:
            plotxdata = np.empty_like(plotydata)
            plotxdata = np.arange(0,plotydata.shape[0])
        else:
            plotxdata = np.asarray(xdata)
                    
        # Plot the data with a gray grid (if required) with the required axis
        # scaling.
        if yerr is not None:
            plotaxis.errorbar(plotxdata, plotydata, yerr=yerr, fmt=linefmt,
                              linestyle=linestyle, **kwargs)
        else:
            plotaxis.plot(plotxdata, plotydata, linefmt, linestyle=linestyle,
                          **kwargs)
        if grid:
            plotaxis.grid(color='gray', linestyle=':', linewidth=1)
        plotaxis.set_xscale(xscale)
        plotaxis.set_yscale(yscale)
    
        # Label the plot.
        if title:
            plotaxis.set_title(title)
        if xlabel:
            plotaxis.set_xlabel(xlabel)
        if ylabel:
            plotaxis.set_ylabel(ylabel)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            # Tidy up
            plotfig.clear()
            del plotfig, plotaxis
#             pyplt.clf()
#             pyplt.close()
            plotaxis = None
            
        # Return the axis object in which the plot was created.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_circles(xdata, ydata, radii, pyplt=plt, plotfig=None, plotaxis=None,
            figsize=(10,6), equal_aspect=False, showplot=None, xscale='linear', yscale='linear',
            facecolor='w', circcolor='y ', linefmt='b-', linestyle='-',
            xlabel='', ylabel='', title='', grid=True, filled=False, **kwargs):
    """

    Create an X,Y scatter plot with symbols and fixed or variable-sized circles.
    The function can add a plot to an existing matplotlib axis or create a new one.
        
    :Parameters:
    
    xdata: array-like
        X axis data. This must be a 1-D array, which is
        interpreted as a list of X coordinates.
        It may be set to None, in which case the index
        numbers for the ydata array are used as X coordinates.
    ydata: array-like
        Y axis data. This must be a 1-D array, which is
        interpreted as a list of Y coordinates.
    radii: array-like or number
        Radii of the circles. This must either be a 1-D array containing
        a radius corresponding to each ydata point, or a number corresponding
        to a fixed radius for all points.
        It may be set to None, in which case no circles are drawn.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    xscale: string, optional
        X axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    yscale: string, optional
        Y axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    facecolor, str, optional
        Matplotlib face colour. Default 'w' = white.
    circcolor, str, optional
        Matplotlib line format for drawing the circle.
        Default 'y ' = yellow with no markers.
    linefmt: str, optional
        Matplotlib line format for joining the circles.
        Default 'b-' = blue line with no symbols.
    linestyle: str, optional.
        Matplotlib line style. Default '-' = solid line.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no label.
    ylabel: string, optional
        Optional Y axis label to be shown on the plot.
        The default is no label.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    grid: bool, optional
        If True, draw a grid on the plot. The default is True.
    filled: bool, optional
        If True, draw filled circles. The default is False.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.plot
        For example label, linewidth, linestyle, color, marker, etc...
        See the matplotlib documentation for a full list.
          
    :Returns:
        
    plotaxis: matplotlib axis object
        The matplotlib axis containing the plot.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None and patches is not None:
        # Check whether matplotlib figure and axis objects have been provided.
        # If not then create new ones. By default, the plot will completely
        # fill the figure. Set a flag to remember that a new figure and/or axis
        # object have been created.
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
    
        # Ensure the Y data array is converted to a numpy ndarray object,
        # which has .ndim and .shape() method. The data can only be 1-D or 2-D.
        plotydata = np.asarray(ydata)
        if plotydata.ndim > 2:
            strg = "Data array has %d dimensions. 1 or 2 expected." % \
                plotydata.ndim
            raise TypeError(strg)
        
        # If an X axis array is not provided it defaults to the index
        # numbers of the first dimension of the Y array. Otherwise,
        # ensure the X data array is also converted to a numpy ndarray
        # object.
        if xdata is None:
            plotxdata = np.empty_like(plotydata)
            plotxdata = np.arange(0,plotydata.shape[0])
        else:
            plotxdata = np.asarray(xdata)
                    
        # Plot the circle centres with a gray grid (if required) with the
        # required axis scaling.
        plotaxis.plot(plotxdata, plotydata, linefmt, linestyle=' ')
        if grid:
            plotaxis.grid(color='gray', linestyle=':', linewidth=1)
        plotaxis.set_xscale(xscale)
        plotaxis.set_yscale(yscale)
        
        # Add the circles, if required
        if radii is not None:
            if isinstance(radii, (float,int)):
                radius = radii
                radii = np.ones_like(ydata) * radius
            for xcen, ycen, radius in zip(plotxdata, plotydata, radii):
                if filled:
                    # patches.Circle doesn't work as expected because the circles
                    # are opaque and the alpha attribute fades the edge lines as well as the face.
                    circ = patches.Circle( (xcen,ycen), radius, fill=True,
                                           facecolor=facecolor, alpha=0.5)
                    plotaxis.add_patch( circ )
                else:
                    xcirc = []
                    ycirc = []
                    for iangle in range(0,368,8):
                        rangle = math.radians(float(iangle))
                        xcirc.append(xcen + radius * math.cos(rangle))
                        ycirc.append(ycen + radius * math.sin(rangle))
                    
                    plotaxis.plot(xcirc, ycirc, circcolor, linestyle=linestyle,
                                  **kwargs)
                    # TODO: patches.CirclePolygon does not recognise linestyle
#                     circ = patches.CirclePolygon( (xcen,ycen), radius,
#                                                   resolution=30, fill=False,
#                                                   facecolor=facecolor, linestyle=linestyle )
                    plotaxis.add_patch( circ )
        
        # Label the plot.
        if title:
            plotaxis.set_title(title)
        if xlabel:
            plotaxis.set_xlabel(xlabel)
        if ylabel:
            plotaxis.set_ylabel(ylabel)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            
        # Return the axis object in which the plot was created.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot and matplotlib.patches not available")
        return None

def plot_ellipses(xdata, ydata, majors, minors, tilts, pyplt=plt, plotfig=None, plotaxis=None,
            figsize=(10,6), equal_aspect=False, showplot=None, xscale='linear', yscale='linear',
            facecolor='w', ellipsecolor='y ', linefmt='b-', linestyle='-',
            xlabel='', ylabel='', title='', grid=True, **kwargs):
    """

    Create an X,Y scatter plot with symbols and fixed or variable-sized ellipses.
    The function can add a plot to an existing matplotlib axis or create a new one.
        
    :Parameters:
    
    xdata: array-like
        X axis data. This must be a 1-D array, which is
        interpreted as a list of X coordinates.
        It may be set to None, in which case the index
        numbers for the ydata array are used as X coordinates.
    ydata: array-like
        Y axis data. This must be a 1-D array, which is
        interpreted as a list of Y coordinates.
    majors: array-like or number
        Major axes of the ellipses. This must either be a 1-D array containing
        a major axis corresponding to each ydata point, or a number corresponding
        to a fixed major axis for all points.
        It may be set to None, in which case no ellipses are drawn.
    minors: array-like or number
        Minor axes of the ellipses. This must either be a 1-D array containing
        a minor axis corresponding to each ydata point, or a number corresponding
        to a fixed minor axis for all points.
        It may be set to None, in which case no ellipses are drawn.
    tilts: array-like or number
        Tilt angles in radians.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    xscale: string, optional
        X axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    yscale: string, optional
        Y axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    facecolor, str, optional
        Matplotlib face colour. Default 'w' = white.
    circcolor, str, optional
        Matplotlib line format for drawing the circle.
        Default 'y ' = yellow with no markers.
    linefmt: str, optional
        Matplotlib line format for joining the circles.
        Default 'b-' = blue line with no symbols.
    linestyle: str, optional.
        Matplotlib line style. Default '-' = solid line.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no label.
    ylabel: string, optional
        Optional Y axis label to be shown on the plot.
        The default is no label.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    grid: bool, optional
        If True, draw a grid on the plot. The default is True.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.plot
        For example label, linewidth, linestyle, color, marker, etc...
        See the matplotlib documentation for a full list.
          
    :Returns:
        
    plotaxis: matplotlib axis object
        The matplotlib axis containing the plot.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None:
        # Check whether matplotlib figure and axis objects have been provided.
        # If not then create new ones. By default, the plot will completely
        # fill the figure. Set a flag to remember that a new figure and/or axis
        # object have been created.
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
    
        # Ensure the Y data array is converted to a numpy ndarray object,
        # which has .ndim and .shape() method. The data can only be 1-D or 2-D.
        plotydata = np.asarray(ydata)
        if plotydata.ndim > 2:
            strg = "Data array has %d dimensions. 1 or 2 expected." % \
                plotydata.ndim
            raise TypeError(strg)
        
        # If an X axis array is not provided it defaults to the index
        # numbers of the first dimension of the Y array. Otherwise,
        # ensure the X data array is also converted to a numpy ndarray
        # object.
        if xdata is None:
            plotxdata = np.empty_like(plotydata)
            plotxdata = np.arange(0,plotydata.shape[0])
        else:
            plotxdata = np.asarray(xdata)
                    
        # Plot the circle centres with a gray grid (if required) with the
        # required axis scaling.
        plotaxis.plot(plotxdata, plotydata, linefmt, linestyle=' ')
        if grid:
            plotaxis.grid(color='gray', linestyle=':', linewidth=1)
        plotaxis.set_xscale(xscale)
        plotaxis.set_yscale(yscale)
        
        # Add the ellipses, if required.
        if majors is not None and minors is not None and tilts is not None:
            if isinstance(majors, (float,int)):
                major = majors
                majors = np.ones_like(ydata) * major
            if isinstance(minors, (float,int)):
                minor = minors
                minors = np.ones_like(ydata) * minor
            if isinstance(tilts, (float,int)):
                tilt = tilts
                tilts = np.ones_like(ydata) * tilt
            for xcen, ycen, major, minor, tilt in \
                zip(plotxdata, plotydata, majors, minors, tilts):
    
                # First generate an untilted ellipse relative to the centre
                xellipse1 = []
                yellipse1 = []
                for iangle in range(0,368,8):
                    rangle = math.radians(float(iangle))
                    xee = major * math.cos(rangle) 
                    yee = minor * math.sin(rangle)
                    xellipse1.append(xee)
                    yellipse1.append(yee)
                    
                # Now rotate the ellipse by the tilt angle and add the zero point
                xellipse = []
                yellipse = []
                for xee, yee in zip(xellipse1,yellipse1):
                    xnew = xcen + (xee * math.cos(tilt) - yee * math.sin(tilt))
                    ynew = ycen + (xee * math.sin(tilt) + yee * math.cos(tilt))
                    xellipse.append(xnew)
                    yellipse.append(ynew)
                    
                plotaxis.plot(xellipse, yellipse, ellipsecolor, linestyle=linestyle,
                        **kwargs)
        
        # Label the plot.
        if title:
            plotaxis.set_title(title)
        if xlabel:
            plotaxis.set_xlabel(xlabel)
        if ylabel:
            plotaxis.set_ylabel(ylabel)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            
        # Return the axis object in which the plot was created.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_xycolumn(xdata, ylist, yerrlist=None, pyplt=plt, figsize=(8,10),
                  equal_aspect=False, maxplots=8, maxfigures=8, showplot=True,
                  xscale='linear', yscale='linear', xlabel='', ylabels=None,
                  title='', grid=True, **kwargs):
    """

    Create a figure containing a column of X,Y scatter plots with a
    shared X axis.

    This function creates a new matplotlib figure.
        
    :Parameters:
    
    xdata: array-like
        Common X axis data. This must be a 1-D array. It may be set to
        None, in which case the index numbers for the ylist arrays are
        used as X coordinates.
    ylist: list of array-like objects
        Y axis data. This must be a list of 1-D arrays.
    yerrlist: list of array-like objects, optional
        If specified, a list of Y error bar data.
        If not specified, there will be no error bars.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    maxplots: int, optional
        The maximum number of plots to be included on each figure. If
        ylist is longer than maxplots, new figures will be created to
        contain the additional plots (up to the maximum specified by
        maxfigures). The default value is 8.
    maxfigures: int, optional
        The maximum number of figures to be created by this function.
        The default value is 8.
    showplot: boolean, optional, default=True
        Controls whether the plot is finalised and shown at the end.
        If True the plot is shown. If False, the plot is never shown (useful
        if the caller needs to add something to the plot before finishing it).
    xscale: string, optional
        X axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    yscale: string, optional
        Y axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no label.
    ylabels: list of strs, optional
        List of Y axis labels to be shown on the plot.
        The default is no labels.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    grid: bool, optional
        If True, draw a grid on the plots. The default is True.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.plot
        For example label, linewidth, linestyle, color, marker, etc...
        See the matplotlib documentation for a full list.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None:
        try:
            firsty = np.asarray(ylist[0])
            if firsty.ndim > 1:
                strg = "Y axis array 1 has %d dimensions. 1 expected." % firsty.ndim
                raise TypeError(strg)
        except (IndexError, TypeError):
                strg = "ylist must be a list of arrays, not %s." % \
                    ylist.__class__.__name__
                raise TypeError(strg)
    
        # If an X axis array is not provided it defaults to the index
        # numbers of the first dimension of the Y array. Otherwise,
        # ensure the X data array is also converted to a numpy ndarray
        # object.
        if xdata is None:
            plotxdata = np.empty_like(firsty)
            plotxdata = np.arange(0,firsty.shape[0])
        else:
            plotxdata = np.asarray(xdata)
        if plotxdata.ndim > 1:
            strg = "X axis array has %d dimensions. 1 expected." % plotxdata.ndim
            raise TypeError(strg)
    
        # Create the figure and add the main title if required.
        page = 1
        fig = pyplt.figure(page, figsize=figsize)
        if equal_aspect:
            plt.gca().set_aspect('equal', adjustable='box')
        if title:
            fig.suptitle(title)
    
        # Local function to encapsulate the new page check code within the
        # loop below.
        def _dopage(pyplt, thispage, thissubno, thisfig, maxplots, maxfigures,
                    figsize, equal_aspect, title):
            if thissubno >= maxplots:
                thissubno = 0
                thispage += 1
                if thispage <= maxfigures:
                    thisfig = pyplt.figure(thispage, figsize=figsize)
                    if equal_aspect:
                        plt.gca().set_aspect('equal', adjustable='box')
                    if title:
                        thisfig.suptitle(title)
            return (thispage, thissubno, thisfig)
    
        # Local function to encapsulate the plotting code within the loop below.
        def _doplot(subno, plotno, subcols, thisplt, thisfig, xscale, yscale,
                    xlabel, ylabels, plotxdata, ydata, yerr, **kwargs):
            # Ensure the plot data is converted to a numpy ndarray and
            # has the dimensionality expected.
            plotydata = np.asarray(ydata)
            if plotydata.ndim > 1:
                strg = "Y axis array %d has %d dimensions. 1 expected." % \
                    (subno+1 % plotydata.ndim)
                raise TypeError(strg)
    
            ax = thisfig.add_subplot(subcols, 1, subno+1)
            if ylabels is not None and len(ylabels)>plotno:
                ystrg = ylabels[plotno]
            else:
                ystrg = ''
            plot_xy(plotxdata, plotydata, yerr=yerr,
                    pyplt=thisplt, plotfig=thisfig, plotaxis=ax,
                    xscale=xscale, yscale=yscale, xlabel=xlabel,
                    ylabel=ystrg, title='', grid=True, showplot=False,
                    **kwargs)        
            return
    
        # Step through each data array in the list:
        subno = 0
        plotno = 0
        subcols = len(ylist)
        if subcols > maxplots:
            subcols = maxplots
        if yerrlist is None:
            for ydata in ylist:
                (page, subno, fig) = _dopage(pyplt, page, subno, fig, maxplots,
                                             maxfigures, figsize, equal_aspect, title)
                if page > maxfigures:
                    break
                _doplot(subno, plotno, subcols, plt, fig, xscale, yscale, xlabel,
                        ylabels, plotxdata, ydata, None, **kwargs)
                subno += 1
                plotno += 1
        else:
            for ydata,yerr in zip(ylist,yerrlist):
                (page, subno, fig) = _dopage(pyplt, page, subno, fig, maxplots,
                                             maxfigures, figsize, equal_aspect, title)
                if page > maxfigures:
                    break
                _doplot(subno, plotno, subcols, plt, fig, xscale, yscale, xlabel,
                        ylabels, plotxdata, ydata, yerr, **kwargs)
                subno += 1
                plotno += 1
    
        if showplot:
            # Show the plot and tidy up.
            pyplt.show()
            fig.clear()
            del fig
            pyplt.clf()
            pyplt.close()
    else:
        logger.warning("matplotlib.pyplot not available")

def plot_hist(data, bins=30, equalwidths=True, pyplt=plt, plotfig=None,
              plotaxis=None, figsize=(10,6), equal_aspect=False, showplot=None,
              clip=None, xscale='linear', yscale='linear', xlabel='', ylabel='',
              title='', **kwargs):
    """

    Create an X,Y scatter plot. The function can add a histogram plot
    to an existing matplotlib axis or create a new one.
        
    :Parameters:
    
    data: array-like
        Array of data to be histogrammed.
    bins: int or list of numbers
        Either: An int giving the required number of histogram bins
        Or: A list of values giving the boundaries of the histogram
        bins.
    equalwidths: bool, optional
        If nbins is an int and equalwidths is True, the bin boundaries
        will be calculated to give equal width bins on the plot.
        This is especially useful when the xscale is not linear.
        The default is True.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,6), which will create a 10 inch x 6
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    clip: tuple of 2 floats, optional
        The lower and upper clipping limits expressed as a fraction of
        the data range, such that (0.0, 1.0) displays the full range of
        the data. If not specified, the full range of data is displayed.
    xscale: string, optional
        X axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
        If xscale is set to 'log', bins is an int and equalwidths is True,
        the histogram bins will be recalculated so their boundaries are
        distributed logarithmically.
    yscale: string, optional
        Y axis scaling type: 'linear' or 'log'.
        The default is 'linear'.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no label.
    ylabel: string, optional
        Optional Y axis label to be shown on the plot.
        The default is no label.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    grid: bool, optional
        If True, draw a grid on the plot. The default is True.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.hist
        For example histtype, facecolor, edgecolor, etc...
        See the matplotlib documentation for a full list.
          
    :Returns:
        
    plotaxis: matplotlib axis object
        The matplotlib axis containing the plot.
            
    :Requires:
        
    matplotlib.pyplot

    """
    if pyplt is not None:
        # Ensure the Y data array is converted to a numpy ndarray object,
        # which has .ndim and .shape() method. The data can only be 1-D or 2-D.
        plotdata = np.asarray(data)
    
        # If the X scale is logarithmic, and the bin boundaries have not
        # been provided, recalculate the bin boundaries so they are
        # logarithmic (so the plot has even-sized bins).
        # This is only possible if the data has a sensible range.
        if isinstance(bins, (int,float)) and xscale == 'log' and equalwidths:
            dmin = plotdata.min()
            dmax = plotdata.max()
            # Negative values cannot be scaled logarithmically. In this case
            # we leave Matplotlib to work out what to do (it will give a valid
            # plot but without the equal width bins).
            if dmin > 0.0 and dmax > 0.0:
                lowlog = math.log10(dmin)
                highlog = math.log10(dmax)
                nbins = int(bins)
                bins = np.logspace(lowlog, highlog, nbins)
            else:
                logger.warning( "plot_hist: Equal width bins not possible - " + \
                    "Negative X values and xscale=%s" % xscale )
    
        # If specified, the clipping parameters are expected to be in the range
        # 0-1 and must be in a sensible order.
        if clip is not None:
            if not isinstance(clip, (tuple,list)) or len(clip) < 2:
                strg = "Clipping parameter must be a tuple of 2 values."
                raise TypeError(strg)
            for cl in clip:
                if cl < 0.0 or cl > 1.0:
                    strg = "Clipping parameters %f .. %f " % clip
                    strg += "are expected to be in the range 0.0-1.0."
                    raise ValueError(strg)
            if clip[1] <= clip[0]:
                strg = "Clipping parameters %f .. %f " % clip
                strg += "are in the wrong order."
                raise ValueError(strg)
    
        # Check whether matplotlib figure and axis objects have been provided.
        # If not then create new ones. By default, the plot will completely
        # fill the figure. Set a flag to remember that a new figure and/or axis
        # object have been created.
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
    
        # If clipping parameters have been specified, calculate the limits.
        if clip is not None:
            lower = plotdata.min() + (plotdata.max()-plotdata.min()) * clip[0]
            upper = plotdata.min() + (plotdata.max()-plotdata.min()) * clip[1]
            hrange = (lower,upper)
        else:
            hrange = None
    
        # Histogram the data with the required axis scaling.
        plotaxis.hist(plotdata, bins, range=hrange, **kwargs)
        plotaxis.set_xscale(xscale)
        plotaxis.set_yscale(yscale)
    
        # Label the plot.
        if title:
            plotaxis.set_title(title)
        if xlabel:
            plotaxis.set_xlabel(xlabel)
        if ylabel:
            plotaxis.set_ylabel(ylabel)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            # Tidy up
            plotfig.clear()
            del plotfig, plotaxis
#             pyplt.clf()
#             pyplt.close()
            plotaxis = None
            
        # Return the axis object in which the plot was created.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_image(data, pyplt=plt, plotfig=None, plotaxis=None, figsize=(10,10),
              equal_aspect=False, showplot=None, datatype='matrix',
              datascale='linear', maxplots=16, maxfigures=8, cmap='hot',
              withbar=False, clip=None, xlabel=None, ylabel=None, zlabel=None,
              title='', **kwargs):
    """
    
    Plot N-D data as an image or series of images.

    1-D or 2-D data are plotted in a single matplotlib axis. The
    function can add the plot to an existing axis object or create a
    new one.
    
    3-D data are plotted as series of 2-D slices - each contained in
    a separate matplotlib axis. Large numbers of slices can be displayed
    in multiple matplotlib figures. Since more than one figure and axis
    are involved, there is no point in trying to plot 3-D data within
    an existing axis.

    Currently 1-D, 2-D or 3-D data are supported.

    :Parameters:

    data: array-like
        The 1-D, 2-D or 3-D data to be plotted.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
        Only sensible for 1-D or 2-D data.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is provided, the plot will be created within that
        axis. Providing an axis allows the caller to control exactly
        where a plot appears within a figure.
        If an axis is not provided, a new one will be created filling
        the current figure.
        Only sensible for 1-D or 2-D data.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,10), which will create a 10 inch x 10
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    datatype: string, optional
        The type of data being displayed: 'matrix' or 'image'. 'matrix'
        is better for displaying an array (with square pixel boundaries)
        and 'image' is better for displaying a smoothly interpolated
        image. The default is 'matrix'.
    datascale: string, optional
        How the data are scaled before being displayed:
        If 'linear', the data are not changed before display.
        If 'log', the logarithmic of the data is displayed.
        If 'squared', the square of the data is displayed.
        If 'sqrt', the square root of the data is displayed.
        The default is 'linear'.
    maxplots: int, optional
        The maximum number of plots to be included on each figure. If
        3-D data contains more slices than maxplots, new figures will
        be created to contain the additional plots (up to the maximum
        specified by maxfigures). The default value is 16.
    maxfigures: int, optional
        The maximum number of figures to be created by this function.
        The default value is 8.
    cmap: string, optional
        Name of the matplotlib colour map to be used. Defaults to 'hot'.
    withbar: bool, optional
        Set to True to add a colour bar. The default is False.
    clip: tuple of 2 floats, optional
        The lower and upper clipping limits expressed as a fraction of
        the data range, such that (0.0, 1.0) displays the full range of
        the data. If not specified, the full range of data is displayed.
        NOTE: If datascale is specified, the clipping is applied to the
        scaled data.
    xlabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is no 'Columns' for a single image or no label for
        multiple images.
    ylabel: string, optional
        Optional X axis label to be shown on the plot.
        The default is 'Rows' for a single image or no label for
        multiple images.
    zlabel: string, optional
        Optional Z axis label to be shown on the plot. Only valid for
        3-D images.  The default is 'Slice'.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.imshow
        For example aspect, interpolation, norm, filternorm, etc...
        See the matplotlib documentation for a full list.

    :Requires:
        
    matplotlib.pyplot
    
    """
    if pyplt is not None:
        # Ensure the data array is converted to a numpy ndarray object,
        # which has an .ndim attribute.
        plotdata = np.asarray(data)
    
        # What happens next depends on the dimensionality of the data.
        if plotdata.ndim <= 2:
            # Check whether matplotlib figure and axis objects have been provided.
            # If not then create new ones. By default, the plot will completely
            # fill the figure. Set a flag to remember that a new figure and/or axis
            # object have been created.
            newaxis = False
            if plotfig is None:
                plotfig = pyplt.figure(1, figsize=figsize)
                if equal_aspect:
                    plt.gca().set_aspect('equal', adjustable='box')
                newaxis = True
            if plotaxis is None:
                plotaxis = plotfig.add_subplot(1,1,1)
                newaxis = True
            if newaxis and (showplot is None):
                showplot = True
                
            # The default labels for a single plot are 'Columns' and 'Rows'.
            if xlabel is None:
                xlabel = 'Columns'
            if ylabel is None:
                ylabel = 'Rows' 
            # A 1-D or 2-D image is plotted within a single subplot of a single
            # figure.
            plot_image2D(plotdata, pyplt=pyplt, plotfig=plotfig, plotaxis=plotaxis,
                         datatype=datatype, datascale=datascale, cmap=cmap,
                         withbar=withbar, clip=clip, xlabel=xlabel, ylabel=ylabel,
                         title=title, **kwargs)
            # If this function has created a new figure or axis object, or if
            # showplot has been explicitly set to True, finish the figure by
            # displaying it. (Otherwise the caller is responsible for doing this).
            if showplot:
                pyplt.show()
                # Tidy up
                plotfig.clear()
                del plotfig, plotaxis
#                 pyplt.clf()
#                 pyplt.close()

        elif plotdata.ndim == 3:
            # For 3-D data there is no point in the caller providing a
            # figure or axis object.
            if plotaxis is not None:
                logger.warning( "plot_image: The existing matplotlib axis is ignored " + \
                    "for 3-D data. The plot will start from a new figure." )
    
            # For a multiple plot, the default is to have no X,Y labels
            # and a Z label of 'Slice'.
            if xlabel is None:
                xlabel = ''
            if ylabel is None:
                ylabel = '' 
            if zlabel is None:
                zlabel = 'Slice'
            
            # A 3-D image can be plotted within multiple subplots of multiple
            # figures.
            plot_image3D(plotdata, pyplt=pyplt, plotfig=plotfig, figsize=figsize,
                         equal_aspect=equal_aspect, datatype=datatype,
                         datascale=datascale, maxplots=maxplots,
                         maxfigures=maxfigures, cmap=cmap, withbar=withbar,
                         clip=clip, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                         title=title, **kwargs)
            
        else:
            strg = "Data array has %d dimensions. 1, 2 or 3 supported." % \
                plotdata.ndim
            raise TypeError(strg)
    else:
        logger.warning("matplotlib.pyplot not available")

def plot_images(datalist, pyplt=plt, figsize=(10,10), equal_aspect=False,
                showplot=True, datatype='matrix', datascales=None, cmap='hot',
                withbar=False, clip=None, xlabel='', ylabel='', title='',
                subtitles=None, **kwargs):
    """
    
    Plot multiple sets of 1-D or 2-D data as a series of images on
    the same page. This function creates a new matplotlib figure.

    Currently 1-D or 2-D data are supported.

    :Parameters:

    datalist: list of array-like objects
        A list of data arrays to be plotted on the same page.
        The function will attempt to plot all the data, no matter
        how large the list. Large lists may result in unreadable plots.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the figure.
        The default value is (10,10), which will create a 10 inch x 10
        inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional, default=True
        Controls whether the plot is finalised and shown at the end.
        If True the plot is shown. If False, the plot is never shown (useful
        if the caller needs to add something to the plot before finishing it).
    datatype: string, optional
        The type of data being displayed: 'matrix' or 'image'. 'matrix'
        is better for displaying an array (with square pixel boundaries)
        and 'image' is better for displaying a smoothly interpolated
        image. The default is 'matrix'.
    datascales: list of strings, optional
        How each data array is scaled before being displayed:
        If 'linear', the data are not changed before display.
        If 'log', the logarithmic of the data is displayed.
        If 'squared', the square of the data is displayed.
        If 'sqrt', the square root of the data is displayed.
        The default is 'linear' for all arrays.
        If a list is supplied, it should be the same length as the
        datalist (missing elements are assumed 'linear').
    cmap: string, optional
        Name of the matplotlib colour map to be used. Defaults to 'hot'.
    withbar: bool, optional
        Set to True to add a colour bar. The default is False.
    clip: tuple of 2 floats, optional
        The lower and upper clipping limits expressed as a fraction of
        the data range, such that (0.0, 1.0) displays the full range of
        the data. If not specified, the full range of data is displayed.
        NOTE: If datascales is specified, the clipping is applied to the
        scaled data.
    xlabel: string, optional
        Optional X axis label to be shown on all the plots.
        The default is 'Columns' for a single image or no label for
        multiple images.
    ylabel: string, optional
        Optional X axis label to be shown on all the plots.
        The default is 'Rows' for a single image or no label for
        multiple images.
    title: str, optional
        Optional main title to be shown at the top of the page, if required.
        The default is no title.
    subtitles: list of str, optional
        An optional list of subtitles to be shown at the top of each
        individual plot. If a list is supplied, it should be the same
        length as the datalist (missing elements are assumed blank).
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.imshow
        For example aspect, interpolation, norm, filternorm, etc...
        See the matplotlib documentation for a full list.

    :Requires:
        
    matplotlib.pyplot
    
    """
    if pyplt is not None:
        # Subdivide the figure into a grid of subplots, one for each
        # data array supplied.
        (subrows, subcols, nfigs) = subdivide(len(datalist))
        if nfigs > 1:
            strg = "plot_images: The full plot would require %d " % nfigs
            strg += "figures, but the maximum is 1."
            strg += "Some data will not be plotted."
            logger.warning( strg )
        
        # The default labels for a single image are 'Columns' and 'Rows'.
        if len(datalist) < 2:
            if not xlabel:
                xlabel = 'Columns'
            if not ylabel:
                ylabel = 'Rows' 
    
        # Create the figure and add the main title if required.
        fig = pyplt.figure(1, figsize=figsize)
        if equal_aspect:
            plt.gca().set_aspect('equal', adjustable='box')
        if title:
            fig.suptitle(title)
    
        # Step through each data array in the list:
        subno = 0
        for data in datalist:
            if (datascales is not None) and (len(datascales) > subno):
                datascale = datascales[subno]
            else:
                datascale = 'linear'
            # Ensure the plot data is converted to a numpy ndarray and
            # has the dimensionality expected.
            plotdata = np.asarray(data)
            if plotdata.ndim > 2:
                strg = "Data array has %d dimensions. 1 or 2 expected." % \
                    plotdata.ndim
                raise TypeError(strg)
    
            ax = fig.add_subplot(subrows, subcols , subno+1)
            if (subtitles is not None) and (len(subtitles) > subno):
                strg = subtitles[subno]
            else:
                strg = ''
            plot_image2D(plotdata, pyplt=pyplt, plotfig=fig, plotaxis=ax,
                         datatype=datatype, datascale=datascale, cmap=cmap,
                         withbar=withbar, clip=clip, xlabel=xlabel, ylabel=ylabel,
                         title=strg, **kwargs)
            subno += 1
    
        if showplot:
            pyplt.show()
            # Tidy up
            fig.clear()
            del fig
#             pyplt.clf()
#             pyplt.close()
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_image2D(data, pyplt=plt, plotfig=None, plotaxis=None, figsize=(10,10),
                 equal_aspect=False, showplot=None, datatype='matrix',
                 datascale='linear', cmap='hot', withbar=False, clip=None,
                 xlabel='Columns', ylabel='Rows', title='', **kwargs):
    """
        
    Plot a 1-D or 2-D image within the given matplotlib axis.
    The function can add an image plot to an existing matplotlib
    axis or create a new one..

    :Parameters:

    data: array-like
        The 2-D data to be plotted. 1-D data is also accepted,
        although it will be plotted as a narrow line of colour.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot.
        If a figure is not provided, a new one will be created.
    plotaxis: matplotlib axis object, optional
        Axis on which to add the plot.
        If an axis is not provided, a new one will be created,
        filling the current figure.
    figsize: tuple of 2 ints, optional
        If a new matplotlib figure is created, the size of the
        figure. The default value is (10,10), which will create
        a 10 inch x 10 inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional
        Controls whether the plot is finalised and shown at the end.
        If True the plot is always shown.
        If False, the plot is never shown.
        If not specified, the plot is only shown when a new axis
        has been created.
    datatype: string, optional
        The type of data being displayed: 'matrix' or 'image'.
        'matrix' is better for displaying an array (with square
        pixel boundaries) and 'image' is better for displaying
        a smoothly interpolated image. The default is 'matrix'.
    datascale: string, optional
        How the data are scaled before being displayed:
        If 'linear', the data are not changed before display.
        If 'log', the logarithmic of the data is displayed.
        If 'squared', the square of the data is displayed.
        If 'sqrt', the square root of the data is displayed.
        The default is 'linear'.
    cmap: string, optional
        Name of the matplotlib colour map to be used. Defaults to
        'hot'.
    withbar: bool, optional
        Set to True to add a colour bar. The default is False.
    clip: tuple of 2 floats, optional
        The lower and upper clipping limits expressed as a fraction of
        the data range, such that (0.0, 1.0) displays the full range of
        the data. If not specified, the full range of data is displayed.
        NOTE: If datascale is specified, the clipping is applied to the
        scaled data.
    xlabel: string, optional
        Optional X axis label to be shown on the plot. The default
        is 'Columns'.
    ylabel: string, optional
        Optional X axis label to be shown on the plot. The default
        is 'Rows'.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.imshow
        For example aspect, interpolation, norm, filternorm, etc...
        See the matplotlib documentation for a full list.

    :Requires:
        
    matplotlib.pyplot
            
    """
    if pyplt is not None:
        # Ensure the plot data is converted to a numpy ndarray and
        # has the dimensionality expected.
        plotdata = np.asarray(data)
        if plotdata.ndim > 2:
            strg = "Data array has %d dimensions. 1 or 2 expected." % plotdata.ndim
            raise TypeError(strg)
        elif plotdata.ndim == 1:
            # A 1-D image can be converted to a 2-D image with a single row.
            plotdata.shape = (1, plotdata.shape[0])
            # A 1-D image is best displayed with 2 extra line breaks
            # at the end of the title (otherwise matplotlib overlaps
            # the labels).
            if title:
                title += "\n\n"
                
        # Scale the data if necessary. Log or sqrt scaling will cause negative
        # values to be filled with the minimum.
        if datascale == 'log':
            logdata = ma.log10(plotdata)
            smin = ma.min(logdata)
            del plotdata
            plotdata = logdata.filled(smin)
        elif datascale == 'square':
            sqdata = plotdata * plotdata
            del plotdata
            plotdata = sqdata
        elif datascale == 'sqrt':
            sqrtdata = ma.sqrt(plotdata)
            smin = ma.min(sqrtdata)
            del plotdata
            plotdata = sqrtdata.filled(smin)
    
        # If specified, the clipping parameters are expected to be in the range
        # 0-1 and must be in a sensible order.
        if clip is not None:
            if not isinstance(clip, (tuple,list)) or len(clip) < 2:
                strg = "Clipping parameter must be a tuple of 2 values."
                raise TypeError(strg)
            for cl in clip:
                if cl < 0.0 or cl > 1.0:
                    strg = "Clipping parameters %f .. %f " % clip
                    strg += "are expected to be in the range 0.0-1.0."
                    raise ValueError(strg)
            if clip[1] <= clip[0]:
                strg = "Clipping parameters %f .. %f " % clip
                strg += "are in the wrong order."
                raise ValueError(strg)
    
        # Check whether matplotlib figure and axis objects have been provided.
        # If not then create new ones. By default, the plot will completely
        # fill the figure. Set a flag to remember that a new figure and/or axis
        # object have been created.
        newaxis = False
        if plotfig is None:
            plotfig = pyplt.figure(1, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
            newaxis = True
        if plotaxis is None:
            plotaxis = plotfig.add_subplot(1,1,1)
            newaxis = True
        if newaxis and (showplot is None):
            showplot = True
        
        # Display the data as an image with the origin at
        # the bottom.
        # NOTE: ax.imshow() will plot the image but will cause the
        # colorbar() function to fail.
        if datatype == 'matrix':
            imgplot = pyplt.matshow(plotdata, fignum=0, origin='lower', **kwargs)
        else:  
            imgplot = pyplt.imshow(plotdata, origin='lower', **kwargs)
    
        # If clipping parameters have been specified, calculate the limits.
        if clip is not None:
            lower = plotdata.min() + (plotdata.max()-plotdata.min()) * clip[0]
            upper = plotdata.min() + (plotdata.max()-plotdata.min()) * clip[1]
            imgplot.set_clim(lower, upper)
    
        imgplot.set_cmap(cmap)
        if withbar:
            # NOTE: imgplot.colorbar() does not work!
            pyplt.colorbar()
            
        # Label the plot
        if title:
            plotaxis.set_title(title)
        if xlabel:
            plotaxis.set_xlabel(xlabel)
        if ylabel:
            plotaxis.set_ylabel(ylabel)
    
        # If this function has created a new figure or axis object, or if
        # showplot has been explicitly set to True, finish the figure by
        # displaying it. (Otherwise the caller is responsible for doing this).
        if showplot:
            pyplt.show()
            # Tidy up
            plotfig.clear()
            del plotfig, plotaxis
#             pyplt.clf()
#             pyplt.close()
            plotaxis = None
            
        # Delete the temporary plot data
        del plotdata
            
        # Return the axis object in which the plot was created.
        return plotaxis
    else:
        logger.warning("matplotlib.pyplot not available")
        return None

def plot_image3D(data, pyplt=plt, plotfig=None, figsize=(10,10),
                 equal_aspect=False, showplot=True,
                 datatype='matrix', datascale='linear', maxplots=16,
                 maxfigures=8, cmap='hot', withbar=False, clip=None,
                 xlabel='', ylabel='', zlabel='Slice', title='', **kwargs):
    """
        
    Plot a 3-D image as a series of 2-D slices within one or more
    matplotlib figures and axes.
    
    This function requires or creates a new matplotlib figure.

    :Parameters:

    data: array-like
        The 3-D data to be plotted.
    pyplt: matplotlib pyplot object, optional
        A top level matplotlib pyplot object.
        Defaults to the pyplot object imported by the miriplot module.
    plotfig: matplotlib figure object, optional
        Figure on which to add the plot. The figure should be
        newly created.
        If a figure is not provided, a new one will be created.
    figsize: tuple of 2 ints, optional
        When a new matplotlib figure is created, the size of the
        figure. The default value is (10,10), which will create
        a 10 inch x 10 inch window on the display.
    equal_aspect: boolean, optional
        If True, modifies the figure size to keep the X/Y aspect
        ratio constant. Setting this to True may change the size requested
        in the figsize parameter. The default is False.
    showplot: boolean, optional, default=True
        Controls whether the plot is finalised and shown at the end.
        If True the plot is shown. If False, the plot is never shown (useful
        if the caller needs to add something to the plot before finishing it).
    datatype: string, optional
        The type of data being displayed: 'matrix' or 'image'.
        'matrix' is better for displaying an array (with square
        pixel boundaries) and 'image' is better for displaying
        a smoothly interpolated image. The default is 'matrix'.
    datascale: string, optional
        How the data are scaled before being displayed:
        If 'linear', the data are not changed before display.
        If 'log', the logarithmic of the data is displayed.
        If 'squared', the square of the data is displayed.
        If 'sqrt', the square root of the data is displayed.
        The default is 'linear'. The same scaling is applied to all
        slices within a 3-D array.
    cmap: string, optional
        Name of the matplotlib colour map to be used. Defaults to
        'hot'.
    withbar: bool, optional
        Set to True to add a colour bar. The default is False.
    clip: tuple of 2 floats, optional
        The lower and upper clipping limits expressed as a fraction of
        the data range, such that (0.0, 1.0) displays the full range of
        the data. If not specified, the full range of data is displayed.
        NOTE: If datascale is specified, the clipping is applied to the
        scaled data.
    xlabel: string, optional
        Optional X axis label to be shown on the plot. The default
        is no label.
    ylabel: string, optional
        Optional X axis label to be shown on the plot. The default
        is no label.
    zlabel: string, optional
        Optional Z axis label to be shown in the title for each slice.
        The default is 'Slice'.
    title: string, optional
        Optional title to be shown above the plot, if required.
        The default is no title.
    \*\*kwargs:
        All other keyword arguments will be passed to matplotlib.imshow
        For example aspect, interpolation, norm, filternorm, etc...
        See the matplotlib documentation for a full list.

    :Requires:
        
    matplotlib.pyplot
            
    """
    if pyplt is not None:
        # Ensure the plot data is converted to a numpy ndarray and
        # has the dimensionality expected.
        plotdata = np.asarray(data)
        if plotdata.ndim != 3:
            strg = "Data array has %d dimensions. 3 expected." % plotdata.ndim
            raise TypeError(strg)
    
        # 3-D data is plotted in a series of 2-D slices, each displayed
        # in its own subplot, up to a maximum of maxplots per figure.
        page = 0
        # If a figure has been provided, start with it, otherwise
        # create a new one.
        if plotfig is None:
            fig = pyplt.figure(page, figsize=figsize)
            if equal_aspect:
                plt.gca().set_aspect('equal', adjustable='box')
        else:
            fig = plotfig
        if title:
            fig.suptitle(title)
            
        # Subdivide each figure into a grid of subplots.
        nslices = plotdata.shape[0]        
        (subrows, subcols, nfigs) = subdivide(nslices, maxplots=maxplots)
        if nfigs > maxfigures:
            strg = "plot_image3D: The full plot would require %d " % nfigs
            strg += "figures, but the maximum is %d. " % maxfigures
            strg += "Some data will not be plotted."
            logger.warning( strg )
            
        # Step through each slice of the 3-D data and display each one
        # in a separate subplot.
        for slicenum in range(0, nslices):
            if zlabel:
                strg = "%s %d" % (zlabel,(slicenum+1))
            else:
                strg = "Slice %d" % (slicenum+1)
            subno = slicenum % maxplots
            if subno == 0:
                # New figure needed
                page += 1
                if page > maxfigures:
                    # Maximum number of figures exceeded. Break out of the for loop.
                    break
                elif page > 1:
                    fig = pyplt.figure(page, figsize=figsize)
                    if equal_aspect:
                        plt.gca().set_aspect('equal', adjustable='box')
                    if title:
                        fig.suptitle(title)
            ax = fig.add_subplot(subrows, subcols , subno+1)
            plotslice = plotdata[slicenum,:,:]
            plot_image2D(plotslice, pyplt=pyplt, plotfig=fig, plotaxis=ax,
                         datatype=datatype, datascale=datascale, cmap=cmap,
                         withbar=withbar, clip=clip, xlabel=xlabel, ylabel=ylabel,
                         title=strg, **kwargs)
        if showplot:
            pyplt.show()
            # Tidy up
            fig.clear()
            del fig
#             pyplt.clf()
#             pyplt.close()
    else:
        logger.warning("matplotlib.pyplot not available")

#
# Other display utilities can be added here. Here are some examples...
#

def plot( *args ):
    """
    
    A general purpose plotting function. Create a plot displaying the
    contents of all the objects listed. This function works in a similar
    manner to the Python print command, in that all objects provided as
    arguments are combined to make a single plot. The last object in the
    list can be used to give the overall plot a title. For example:
    
    plot( obj1, obj2, obj3, 'Overall title of the plot')
    
    Each object is displayed within a subplot contained within a
    single figure. The subplots are arranged in a pattern as close
    as possible to a square, so for example a list of 8 objects would
    be plotted in a 3x3 pattern with one subplot blank.
    
    Two kinds of object may be plotted:
    
    1) For simple list, tuple and array-like objects, the data
    contained within those objects is plotted as an X,Y plot.
    Each plot is labelled with the class name of the object.
    
    2) For general objects, the .plot() method is used to create the
    plot, labelled with a title created by the .get_title() method.
    NOTE: The first argument of the plot() method is assumed to be
    a matplotlib axis.
    
    Primitive objects and objects without a .plot() method are ignored.
        
    :Parameters:
        
    args: object
        An object whose contents are to be plotted. Several may be provided.
        If the last argument is a string it is interpreted as the overall
        title for the plot.  

    :Global variables:
        
    plt: matplotlib pyplot object
        The top level pyplot object, imported by the miriplot module.
            
    :Requires:
        
    matplotlib.pyplot
            
    """
    if test_available():
        largs = list(args)
        # Create the top level figure and add an overall title (if one has
        # been specified).
        # NOTE: If only one string is provided this is assumed to be something
        # to be plotted, rather than a title for a null plot.
        fig = plt.figure(figsize=(10,10))
    #     plt.gca().set_aspect('equal', adjustable='box')
        if len(largs) > 1:
            lastarg = largs[-1]
            if isinstance(lastarg, str):
                # Add the title and remove it from the arguments
                fig.suptitle(lastarg)
                #largs.remove(lastarg) # Causes problem with JWST data models
                largs = largs[:-1]
            
        # Arrange subplots into a square grid filling only one figure.
        naxes = len(largs)
        (pltrows, pltcols, nfigs) = subdivide( naxes )
        if nfigs > 1:
            strg = "plot: The full plot would require %d " % nfigs
            strg += "figures, but the maximum is 1. "
            strg += "Some data will not be plotted."
            logger.warning( strg )
        
        # The maximum title length depends on the number of columns.
        tlen = 80 // pltrows
        
        # The font size depends on the number of rows and columns.
        if naxes < 2:
            fontsize = 14
        elif naxes < 4:
            fontsize = 12
        elif naxes < 8:
            fontsize = 10
        else:
            fontsize = 8
                
        # Plot each object within its own subplot. Skip any object
        # that doesn't have a .plot() method.
        count = 1
        for obj in largs:
                
            # Create a new subplot/
            ax = fig.add_subplot(pltrows, pltcols, count)
            
            # Check the data type.
            if isinstance(obj, str):
                
                # Text is simply written within an otherwise empty axis.
                ax.set_title( obj.__class__.__name__ )
                ax = plot_text(obj, plotfig=fig, plotaxis=ax, fontsize=fontsize)
            elif isinstance(obj, (float,int,bool)):
                
                # A scalar object is converted to text and centred within
                # an otherwise empty axis.
                ax.set_title( obj.__class__.__name__ )
                ax = plot_text(str(obj), xpos=0.5, ypos=0.5,
                               horizontalalignment='center',
                               verticalalignment='center',
                               fontsize=(fontsize+4),
                               plotfig=fig, plotaxis=ax)
            else:
                if isinstance(obj, (list,tuple,np.ndarray)):
                    # Array-like objects are labelled with their class name,
                    # unless they have a .get_title() method which returns
                    # a title. The kind of plot generated depends on the
                    # dimensionality of the data.
                    if hasattr(obj, 'get_title'):
                        strg = obj.get_title()
                        slines = strg.split("\n")
                        tstrg = slines[0]
                    else:
                        tstrg = obj.__class__.__name__
                    ax.set_title( tstrg[:tlen] )
                    obj2 = np.asarray(obj)  # Just in case of a list or tuple
                    if obj2.ndim > 1:
                        # Image plot
                        ax = plot_image(obj, plotfig=fig, plotaxis=ax)
                    else:
                        # X,Y plot
                        ax = plot_xy(None, obj, plotfig=fig, plotaxis=ax)
                    del obj2
                else:
                    # For miscellaneous objects, look for a .get_title()
                    # method to obtain a plot title and a .plot() method
                    # to plot them.
                    if hasattr(obj, 'get_title'):
                        strg = obj.get_title()
                        slines = strg.split("\n")
                        tstrg = slines[0]
                    else:
                        tstrg = obj.__class__.__name__
                    ax.set_title( tstrg[:tlen] )
                    if hasattr(obj, 'plot'):
                        # Aattempt to invoke the plot method with the
                        # matplotlib axis as first argument.
                        ax = obj.plot(ax)
                    else:
                        # Ignore objects without a .plot() method.
                        pass
            count += 1     
        plt.show()
        del largs
    else:
        logger.warning("matplotlib plotting utilities are not available")


#
# A minimal test is run when this file is run as a main program.
# For a more substantial test see miritools/tests/test_miriplot.py.
#
if __name__ == '__main__':
    import os
    print( "Testing the miriplot module\n" )

    print( "Plain text..." )
    plot_text("Hello, this is some text.\nThis is another line of text.",
              fontsize=18, color='blue', title="Plot with nothing but text")
 
    print( "XY plots..." )
    a = np.linspace(0, 100, 64)
    err = np.ones_like(a) * 1000.0
    err2 = np.ones_like(a) / 10.0
    b = a * a
    c = np.sin( a/10.0 )
    d = np.cos( a/20.0 )
    plot_xy( a, b, xlabel='Value', ylabel='Value squared',
             title='Test plot of Y = X^2', linewidth=4)
    plot_xy( a, b, xlabel='Value', ylabel='Value squared', yscale='log',
             title='Test plot of Y = X^2 (logarithmic scaling)', marker='o')
    plot_xy( a, b, yerr=err, xlabel='Value', ylabel='Value squared',
             title='Test plot of Y = X^2 with error bars', linewidth=4)
    plot_xy( a, b, yerr=err, xlabel='Value', ylabel='Value squared',
             yscale='log',
             title='Test plot of Y = X^2 with error bars (logarithmic scaling)',
             marker='o')
    close()
 
    ylist = [b, c, d]
    plot_xycolumn( a, ylist, yerrlist=None, xlabel='Value',
                   ylabels=('Square', 'Sine', 'Cosine'),
                   title='Test plot of XY columns')
    yerrlist = [err, err2, err2]
    plot_xycolumn( a, ylist, yerrlist=yerrlist, xlabel='Value',
                   ylabels=('Square', 'Sine', 'Cosine'),
                   title='Test plot of XY columns with error bars')
    close()
 
    print( "Histogram plots..." )
    plot_hist( b, bins=30, xlabel='Value squared', ylabel='Frequency',
               title='Test histogram of Y = X^2' )
    plot_hist( b+0.1, bins=30, xlabel='Value squared', ylabel='Frequency',
               xscale='log', title='Test histogram of Y = X^2 in log space' )
    close()
 
    print( "Plots of 1-D image..." )
    plot_xy( None, a, xlabel='Columns', ylabel='Value',
             title='Test plot of 1-D image as XY graph')
    plot_image( a, cmap='hot', withbar=True, xlabel='Columns',
                ylabel='Row', title='Test plot of 1-D image')
 
    imagelist = [a, b]
    plot_images( imagelist, title='Test plot of multiple images',
                 subtitles=('Image 1', 'Image 2'))
    close()
         
    print( "Plots of 2-D image..." )
    a.shape = (8, 8)
    plot_xy( None, a, xlabel='Columns', ylabel='Value',
             title='Test plot of 2-D image as XY graph')
    plot_image(a, cmap='hot', withbar=True, xlabel='Columns', ylabel='Rows',
               title='Test plot of 2-D image', interpolation='nearest')    
    plot_image(a, datascale='log', cmap='hot', withbar=True, xlabel='Columns',
               ylabel='Rows',
               title='Test plot of 2-D image (logarithmic scaling)',
               interpolation='nearest')
    plot_image( a, datascale='square', cmap='hot', withbar=True,
                xlabel='Columns', ylabel='Rows',
                title='Test plot of 2-D image (square scaling)',
                interpolation='nearest')
    plot_image( a, datascale='sqrt', cmap='hot', withbar=True,
                xlabel='Columns', ylabel='Rows',
                title='Test plot of 2-D image (square root scaling)',
                interpolation='nearest')
    close()
 
    print( "Plots of 3-D image..." )
    a.shape = (4, 4, 4)
    plot_image( a, cmap='hot', withbar=True, xlabel='Columns', ylabel='Rows',
                title='Test plot of 3-D image')
    plot_image( a, datascale='log', cmap='hot', withbar=True, xlabel='Columns',
                ylabel='Rows',
                title='Test plot of 3-D image (logarithmic scaling)')
    close()
 
    b = np.tile(a, 10)
    b.shape = (40, 4, 4)
    plot_image( b, datatype='image', cmap='hot', withbar=False, xlabel='',
                maxfigures=2, ylabel='', title='Test plot of large 3-D image')
    b.shape = (40, 16)
    close()

#     print( "Ellipse..." )
#     xcen = 0.0
#     ycen = 0.0
#     major = 5.0
#     minor = 3.0
#     for tilt in range(0,366,30):
#         orient = math.radians(float(tilt))
#         plot_ellipses([xcen], [ycen], [major], [minor], [orient],
#                       title="Tilt angle %.2f degrees" % tilt,
#                       equal_aspect=True)
 
    print( "Testing subdivide..." )
    for maxplt in (None, 20):
        for rratio in (1.0, 1.5, 1.61, 0.75, 0.5, 0.0):
            for nplt in (100, 64, 20, 10, 6, 3, 2, 1, 0):
                (rows, cols, nfigs) = subdivide( nplt, ratio=rratio,
                                          maxplots=maxplt )
                print( "maxplots=", maxplt, ":", nplt, \
                    "plots with ratio", rratio, \
                    "divided into rows=", rows, "x cols=", cols, \
                    "with", nfigs, "figures." )

    #This test requires objects that have .plot() methods.
    print( "Miscellaneous objects..." )
    a = np.linspace(0, 100, 64)
    b = (1,2,3,4,3,2,1)
    c = 42
    
    from miri.datamodels.miri_filters import MiriFilter
    wavelength = np.linspace(2.0, 20.0, 180)
    # Crude top hat filter
    transmission = 0.5*np.ones((180,))
    transmission[:50] = 0.0
    transmission[-50:] = 0.0
    filter_table = []
    for (wav,trans) in zip(wavelength, transmission):
        filter_table.append( [wav,trans] )

    f1 = MiriFilter( filter_table=filter_table, filter_name='FND',
                     filter_type='Blocking' )

    #from miri.datamodels.miri_measurement import MiriMeasurement, \
    #    ascii_to_measurement
    #import miri.datamodels.data
    #datapath = miri.datamodels.data.__path__[0]
    #test_file_name = os.path.join(datapath,"example_measurement.txt")
    #mv = ascii_to_measurement(test_file_name,
    #                          title="Test of MeasuredVariable class",
    #                          name='Dark current', unit='electrons',
    #                          vname='Temperature', vunit='K',
    #                          interptype='LINEAR')
 
    # Plot a collection of objects. Those without .plot() methods are ignored.
    # The last argument is taken as the plot title.
    plot( a, b, 'Two arrays only' )
    plot( c, 'Text', 'Two scalars only' )
    plot( f1, 'MiriFilter object')
    plot( str(f1.get_meta_str()), 'Metadata only' )
    plot( a, b, c, 'Text', f1.values, str(f1.get_meta_str()),
          'Two arrays, two scalars and metadata' )
    #plot( a, b, c, 'Text', mv.values, str(mv.get_meta_str()),
    #      'Two arrays, two scalars and a measurement' )
    close()

    del a, b, c
    print( "Test finished." )
