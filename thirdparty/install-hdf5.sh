#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

cd hdf5

./configure --enable-parallel --enable-shared 
make
make install



