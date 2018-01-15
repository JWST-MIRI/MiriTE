#!/bin/bash

cd ./source
#mirimod=$(echo tools datamodels simulators)
#for module in $mirimod; do
#	echo linking miri/$module
#	ln -fs ../../$module/doc/source/release_notes.rst ./release_notes_$module.rst  &&\
#	cd ./reference
#	ln -fs ../../../$module/doc/source/reference ./$module  &&\
#	cd ..
#done
#
# tools
echo linking miri/tools
ln -fs ../../tools/doc/source/release_notes.rst ./release_notes_tools.rst
mkdir -p reference
cd ./reference
ln -fs ../../../tools/doc/source/reference ./tools
cd ..
# datamodels
echo linking miri/datamodels
ln -fs ../../datamodels/doc/source/release_notes.rst ./release_notes_datamodels.rst
mkdir -p reference
cd ./reference
ln -fs ../../../datamodels/doc/source/reference ./datamodels
cd ..
# simulators
echo linking miri/simulators
ln -fs ../../simulators/doc/source/release_notes.rst ./release_notes_simulators.rst
mkdir -p reference
cd ./reference
ln -fs ../../../simulators/doc/source/reference ./simulators
cd ..

cd ..
