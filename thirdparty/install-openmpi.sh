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
# set LD_LIBRARY_PATH and PATH according to different versions

if grep -q $intallpath ~/.bashrc; then
  echo -n
else
  echo "Writing path to .bashrc"
  echo "export PATH=$intallpath/bin:\$PATH" >> ~/.bashrc

  if [ $version == "Linux" ]; then
    echo "export LD_LIBRARY_PATH=$intallpath/lib:\$LD_LIBRARY_PATH" >> ~/.bashrc
  elif [ $version == "Darwin" ]; then
    echo "export DYLD_LIBRARY_PATH=$intallpath/lib:\$DYLD_LIBRARY_PATH" >> ~/.bashrc
  fi
fi

exit 1


