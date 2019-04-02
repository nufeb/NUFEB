#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Uninstalling nufebFoam.."
nufebfoamDir=$PWD


cd $nufebfoamDir/src
rm Make/options

wclean dragModels
wclean 

