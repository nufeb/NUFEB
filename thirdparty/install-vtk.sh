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
 -DBUILD_TESTING=OFF \
 -DCMAKE_INSTALL_PREFIX:PATH=$intallpath ../

make -j4
make install

version=`uname`
# set LD path according to different versions

if grep -q $intallpath ~/.bashrc; then
  echo -n
else
  echo "Writing path to .bashrc"
  if [ $version == "Linux" ]; then
    echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$intallpath/lib" >> ~/.bashrc
  elif [ $version == "Darwin" ]; then
    echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$intallpath/lib" >> ~/.bashrc
  fi
fi

exit 0
