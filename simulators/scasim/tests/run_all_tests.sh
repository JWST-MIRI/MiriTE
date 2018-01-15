#!/bin/bash -f
#
# Run the unit tests for SCASim, and save the results
# to log files.
#
python test_amplifier.py >>test_amplifier.log 2>&1
python test_cosmic_ray.py >>test_cosmic_ray.log 2>&1
python test_detector.py >>test_detector.log 2>&1
python test_exposure_data.py >>test_exposure_data.log 2>&1
python test_sensor_chip_assembly.py >>test_sensor_chip_assembly.log 2>&1
