Sensor chip assembly module (:mod:`miri.simulators.scasim.sensor_chip_assembly`)
================================================================================

.. module:: miri.simulators.scasim.sensor_chip_assembly

Description
~~~~~~~~~~~
The sensor_chip_assembly module contains the SensorChipAssembly class (the
top level class of the MIRI Sensor Chip Assembly simulator) plus associated
functions. The top level simulate_sca function runs the SCA simulator;
reading a detector illumination input file and writing an output file.
The simulate_sca_fromdata function provides a way of running the SCA
simulator directly from data, without the need to write a detector
illumination file. The overall class hierarchy is::

    SensorChipAssembly
        CosmicRayEnvironment
            CosmicRay
        MiriIlluminationModel
        MiriBadPixelMaskModel
        MiriExposureModel or ExposureData
           MiriMeasuredModel
        DetectorArray
            ImperfectIntegrator
                PoissonIntegrator
            MiriMeasurement
        ParameterFileManager
        MiriQuantumEfficiency
            MiriFilter

The MiriMeasurement, MiriQuantumEfficiency and ParameterFileManager classes are
obtained from the miri.datamodels package. The software also uses the following
configuration files to describe the properties of the MIRI Sensor Chip
Assembly and its environment::

    cosmic_ray_properties.py
    detector_properties.py

Objects
~~~~~~~
.. autoclass:: SensorChipAssembly
   :members:

Functions
~~~~~~~~~
Global functions simulate_sca, simulate_sca_list and simulate_pipeline 
are deprecated. Global function simulate_sca_fromdata has been removed.
