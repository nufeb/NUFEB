#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory
currentDir=$PWD
cd hdf5

./configure --enable-parallel --enable-shared 
make
make install

version=`uname`
# set LD path according to different versions
if [ $version == "Linux" ] 
then
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$currentDir/hdf5/hdf5/lib/" >> ~/.bashrc
elif [ $version == "Darwin" ] 
then
echo "export DYLD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$currentDir/hdf5/hdf5/lib/" >> ~/.bashrc
fi

