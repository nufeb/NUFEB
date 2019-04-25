#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD
cd vtk
mkdir vtk-build
cd vtk-build
mkdir vtk-8.0
intallpath=$PWD/vtk-8.0
echo $intallpath

cmake \
 -DCMAKE_INSTALL_PREFIX:PATH=$intallpath ../

make -j4
make install

echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$currentDir/vtk/vtk-build/vtk-8.0/lib" >> ~/.bashrc
