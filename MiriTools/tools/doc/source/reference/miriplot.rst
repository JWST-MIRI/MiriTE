Plotting module (:mod:`miri.tools.miriplot`)
============================================

.. module:: miri.tools.miriplot

Description
~~~~~~~~~~~
This module contains general purpose functions for plotting MIRI data. 
In fact, the functions can be used for plotting any numpy array data 
which matches the function parameters described. It is designed to be 
used as a toolbox by other MIRI software. The "plot" function may be run 
directly from a Python shell.

Note: The plotting functions use matplotlib. It is important to use the
close() function to tidy up and remove old figures after making a series
of plots. Failure to tidy up regularly may result in old figures popping
back into view on the next call to showplot().

Functions
~~~~~~~~~
.. autofunction:: subdivide
.. autofunction:: new_figure
.. autofunction:: add_subplot
.. autofunction:: show_plot
.. autofunction:: close
.. autofunction:: plot_text
.. autofunction:: plot_xy
.. autofunction:: plot_xycolumn
.. autofunction:: plot_hist
.. autofunction:: plot_image
.. autofunction:: plot_images
.. autofunction:: plot_image2D
.. autofunction:: plot_image3D
.. autofunction:: plot
