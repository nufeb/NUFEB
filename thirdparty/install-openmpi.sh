#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD
tar xvzf openmpi-3.0.6.tar.gz

cd openmpi-3.0.6
mkdir ompi-build
intallpath=$PWD/ompi-build
echo $intallpath

./configure --prefix=$intallpath

make all install

version=`uname`
# set DYLD_LIBRARY_PATH and PATH according to different versions
if [ $version == "Linux" ] 
then
echo "export LD_LIBRARY_PATH=$currentDir/openmpi-3.0.6/ompi-build/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export PATH=$currentDir/openmpi-3.0.6/ompi-build/bin:\$PATH" >> ~/.bashrc
elif [ $version == "Darwin" ] 
then
echo "export DYLD_LIBRARY_PATH=$currentDir/openmpi-3.0.6/ompi-build/lib:\$DYLD_LIBRARY_PATH" >> ~/.bashrc
echo "export PATH=$currentDir/openmpi-3.0.6/ompi-build/bin:\$PATH" >> ~/.bashrc
fi


