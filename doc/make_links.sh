#!/bin/bash

cd ./source

mirimod=$(echo datamodels tools)
for module in $mirimod; do
        echo linking miri/$module
        cd ./reference
        ln -fs ../../../$module/doc/source/reference ./$module  &&\
        cd ..
done

# Simulators needs to be handled separately because of scasim
echo linking miri/simulators
mkdir -p ./reference/simulators
cd ./reference/simulators
ln -fs ../../../../simulators/doc/source/reference/* .
mkdir -p ./scasim
cd ./scasim
# This link might already exist
ln -fs ../../../../../simulators/scasim/doc/source/reference/* .
cd ../../..


cd ..
