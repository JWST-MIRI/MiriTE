#!/bin/bash
cd ./source/reference
ln -fs ../../../../doc/source/reference/find_simulator_file.rst .
ln -fs ../../../../doc/source/reference/integrators.rst .
#
ln -fs ../../../../data/detector/bad_pixel.txt .
ln -fs ../../../amplifier_properties.py .
ln -fs ../../../cosmic_ray_properties.py .
ln -fs ../../../detector_properties.py .
#
ln -fs ../../../../data/detector/dark_currentIM.txt .
ln -fs ../../../../data/detector/dark_currentLW.txt .
ln -fs ../../../../data/detector/dark_currentSW.txt .
#
ln -fs ../../../../data/detector/qe_measurementIM.txt .
ln -fs ../../../../data/detector/qe_measurementLW.txt .
ln -fs ../../../../data/detector/qe_measurementSW.txt .
#
ln -fs ../../../../data/amplifiers/read_noiseIM.txt .
ln -fs ../../../../data/amplifiers/read_noiseLW.txt .
ln -fs ../../../../data/amplifiers/read_noiseSW.txt .
cd ../..
