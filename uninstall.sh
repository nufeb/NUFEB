#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

echo "Uninstalling NUFEB.."
currentDir=$PWD

cd $currentDir/lammps/src
make no-all
make clean-all
