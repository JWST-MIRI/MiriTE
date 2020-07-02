Filter Utilities Module (:mod:`datamodels.miri_filters`)
========================================================

.. module:: miri.datamodels.miri_filters

Description
~~~~~~~~~~~
The miri_filters module is used to read, store and apply a filter,
which is described as transmission function which varies with wavelength.
The  transmission can vary from 0.0 (no flux transmitted) to 1.0 (flux 
unaffected).

A filter can also be used to describe the quantum efficiency of a 
detector, where the transmission function in this case represents the QE 
as a function of wavelength.

Objects
~~~~~~~
.. autoclass:: MiriFilter
   :members:

.. autoclass:: MiriBandPassFilter
   :members:

.. autoclass:: MiriQuantumEfficiency
   :members:

Functions
~~~~~~~~~
.. autofunction:: ascii_to_filter

Configuration files
~~~~~~~~~~~~~~~~~~~
A Filter object can be stored either within a FITS file or an ASCII 
table. A FITS file contains all the information about a Filter object 
whereas a ASCII file only contains the transmission as a function of 
wavelength.

A quantum efficiency measurement is defined in a similar file format, 
except the transmission column contains QE measurements.

Naming convention
-----------------
MIRI filter measurements are defined in the files 
data/filters/IM-VV-FXXXX.txt where "VV" is a version number and "FXXXX" 
is the name of the filter (e.g. "F2100W").

MIRI quantum efficiency measurements are defined in the files 
data/detector/qe_measurementXXX.txt where "XXX" is the focal plane 
module ID (e.g. "495").

ASCII file format
-----------------
ASCII files contain data in the following format. Lines starting with a 
'#' character are ignored. There is no header, and the file can contain 
two columns (where there is a measurement at a single temperature):

* Column #1: The wavelength in microns.
* Column #2: The quantum efficiency measurement at this wavelength.

or three columns (where there are measurements at room temperature and 
cryogenic temperature):

* Column #1: wavelength (micron) %.2f    
* Column #2: transmission (%) %.5f
* Column #3: transmission (%) %.5f 

.. note:: Column #2 is not used by :mod:`miri.datamodels.filters`, it corresponds to room temperature measurements. Column #3 corresponds to cryogenic temperature measurement in MIRI data and is actually used by :mod:`miri.datamodels.filters`.

FITS File format
----------------
When stored as a FITS file, a MIRI filter model contains the following HDUs

* Primary - metadata (including FILTERTY=filter type,
  WAVECENT=central wavelength and FWHM=FWHM wavelength)
* FILTER - filter transmission table
