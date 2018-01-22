#!/bin/bash

cd ./source
# There is now only one module: scasim
mirimod=$(echo scasim)
for module in $mirimod; do
	echo linking miri/simulators/$module
	ln -fs ../../$module/doc/source/release_notes.rst ./release_notes_$module.rst  &&\
	cd ./reference
	ln -fs ../../../$module/doc/source/reference ./$module  &&\
	cd ..
done
cd ..
