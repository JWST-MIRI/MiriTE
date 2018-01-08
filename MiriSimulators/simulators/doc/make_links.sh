#!/bin/bash

cd ./source
mirimod=$(echo mirimsim scasim)
for module in $mirimod; do
	echo linking miri/MiriSimulators/simulators/$module
	ln -fs ../../$module/doc/source/release_notes.rst ./release_notes_$module.rst  &&\
	cd ./reference
	ln -fs ../../../$module/doc/source/reference ./$module  &&\
	cd ..
done
cd ..
