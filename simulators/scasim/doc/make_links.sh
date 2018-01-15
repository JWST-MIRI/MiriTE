#!/bin/bash
cd ./source/reference
ln -fs ../../../../doc/source/reference/find_simulator_file.rst find_simulator_file.rst
ln -fs ../../../../doc/source/reference/integrators.rst integrators.rst
#
ln -fs ../../../../data/detector/bad_pixel.txt bad_pixel.txt
ln -fs ../../../amplifer_properties.py amplifier_properties.py
ln -fs ../../../cosmic_ray_properties.py cosmic_ray_properties.py
ln -fs ../../../detector_properties.py detector_properties.py
#
ln -fs ../../../../data/detector/dark_currentIM.txt dark_currentIM.txt
ln -fs ../../../../data/detector/dark_currentLW.txt dark_currentLW.txt
ln -fs ../../../../data/detector/dark_currentSW.txt dark_currentSW.txt
#
ln -fs ../../../../data/detector/qe_measurementIM.txt qe_measurementIM.txt
ln -fs ../../../../data/detector/qe_measurementLW.txt qe_measurementLW.txt
ln -fs ../../../../data/detector/qe_measurementSW.txt qe_measurementSW.txt
#
ln -fs ../../../../data/amplifiers/read_noiseIM.txt read_noiseIM.txt
ln -fs ../../../../data/amplifiers/read_noiseLW.txt read_noiseLW.txt
ln -fs ../../../../data/amplifiers/read_noiseSW.txt read_noiseSW.txt
cd ../..
