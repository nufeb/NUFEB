#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD
tar xvzf openmpi-1.10.2.tar.gz

cd openmpi-1.10.2
mkdir ompi-build
intallpath=$PWD/ompi-build
echo $intallpath

./configure --prefix=$intallpath

make all install
