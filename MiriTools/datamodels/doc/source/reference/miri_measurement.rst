Measurement Utilities Module (:mod:`miri.measurements.miri_measurement`)
========================================================================

.. module:: miri.datamodels.miri_measurement

Description
~~~~~~~~~~~
This module contains the MiriMeasurement class, which is used to 
describe any set of measurements which vary in a monotonic manner with 
respect to one variable parameter. For example, dark current and 
read noise are expected to vary with detector temperature. This class 
stores the measurements as a lookup table from which the property at any
parameter setting can be predicted by interpolating the measurements

*NOTE: It is assumed that calibration measurements will be made over the 
full range of operational conditions. The MiriMeasurement class will 
not give good predictions outside the range of measurements because it 
does not extrapolate.*

Objects
~~~~~~~
.. autoclass:: MiriMeasurement
   :members:

Functions
~~~~~~~~~
.. autofunction:: ascii_to_measurement

Configuration files
~~~~~~~~~~~~~~~~~~~
MIRI measured variables are defined in the files 
data/detector/dark_currentFPMSNxxx.txt and 
data/amplifier/read_noiseFPMSNxxx.txt where "SNxxx" is the ID of a 
particular focal plane module.

ASCII files contain data in the following format. Lines starting with a 
'#' character are ignored. The file contains two or more columns:

  * Column #0: The environmental parameter value (e.g. temperature).
  * Column #1: The corresponding variable measurement.
  * Column #2: Another corresponding measurement of the same variable.
  * etc...
  
Columns 2 onwards may be used to store separate measurements of the same 
parameter. For example, one file can be used to store the read noise 
measurements as a function of temperature for each amplifier within the 
same focal plane module.

FITS files are stored in the format specified by the JWST data model.
The measurement data may be found in a binary table HDU called 'MEASUREMENT'.

