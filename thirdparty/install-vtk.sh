#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

mkdir vtk-build
cd vtk-build
mkdir vtk-6.3
intallpath=$PWD/vtk-6.3
echo $intallpath

cmake \
 -DCMAKE_INSTALL_PREFIX:PATH=$intallpath \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTK_USE_SYSTEM_ZLIB:BOOL=ON \
 -DBUILD_SHARED_LIBS:BOOL=OFF  ../vtk

make -j4
make install
