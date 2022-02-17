#!/bin/bash

set -euo pipefail

cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD
cd hdf5 || exit 1

./configure --enable-parallel --enable-shared 
make
make install

version=`uname`
# set LD path according to different versions
intallpath=$currentDir/hdf5/hdf5

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

