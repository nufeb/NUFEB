#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory
currentDir=$PWD
cd hdf5

./configure --enable-parallel --enable-shared 
make
make install

echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$currentDir/hdf5/hdf5/lib/" >> ~/.bashrc


