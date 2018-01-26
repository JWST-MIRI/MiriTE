#!/bin/bash

cd ./source

mirimod=$(echo datamodels tools)
for module in $mirimod; do
        echo linking miri/$module
        ln -fs ../../$module/doc/source/release_notes.rst ./release_notes_$module.rst  &&\
        cd ./reference
        ln -fs ../../../$module/doc/source/reference ./$module  &&\
        cd ..
done

# Simulators needs to be handled separately because of scasim
echo linking miri/simulators
ln -fs ../../simulators/doc/source/release_notes.rst ./release_notes_simulators.rst
ln -fs ../../simulators/scasim/doc/source/release_notes.rst ./release_notes_scasim.rst
mkdir -p ./reference/simulators
cd ./reference/simulators
ln -fs ../../../../simulators/doc/source/reference/* .
mkdir -p ./scasim
cd ./scasim
ln -fs ../../../../../simulators/scasim/doc/source/reference/* .
cd ../../..


cd ..
