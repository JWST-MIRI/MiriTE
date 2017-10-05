Poisson Integrators Module (:mod:`miri.simulators.integrators`)
===============================================================

.. module:: miri.simulators.integrators

Description
~~~~~~~~~~~
This module simulates the behaviour of a detector made of electron-counting
pixels, when the incoming flux is small. It can be reused in any simulation
where a particle integrator obeying Poisson statistics is needed.

The imperfect integrator adds cosmic rays, zero-point and latency effects.

Objects
~~~~~~~
.. autoclass:: PoissonIntegrator
   :members:
   
.. autoclass:: ImperfectIntegrator
   :members:

Functions
~~~~~~~~~
None.
