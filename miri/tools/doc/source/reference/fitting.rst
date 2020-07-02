Fitting tools (:mod:`miri.tools.fitting`)
=========================================================

.. module:: miri.tools.fitting

Description
~~~~~~~~~~~
This module contains general purpose fitting functions for
MIRI data, including:

gaussPlusPoly - Generate a Gaussian plus a polynomial

gaussian - Generate a Gaussian

nonLinFit - Nonlinear least squares fit using the Levenberg-Marquardt method
      in the scipy.optimize.curve_fit method
      this should be possible with any function in the format of
      ``ydata = f(xdata, *params)``
      
LinFit - Class for linear fitting


Objects
~~~~~~~
.. autoclass:: LinFit
   :members:

Functions
~~~~~~~~~
.. autofunction:: gaussPlusPoly
.. autofunction:: gaussian
.. autofunction:: nonLinFit
