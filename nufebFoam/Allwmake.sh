#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing nufebFoam (for mac/linux).."
nufebfoamDir=$PWD

# Determine if the directory of LAMMPS exists or not.
# If not, look for LAMMPS in the default directory.

cd ..
lammpsDir=$PWD/lammps
nufebDir=$PWD

#### Copy package and lib files to LAMMPS directory #####
cd ..
echo "Copying packages to LAMMPS.."
cp -rf $nufebDir/src/* $lammpsDir/src/
cp -rf $nufebDir/lib/* $lammpsDir/lib/
cp -rf $nufebfoamDir/MAKE/* $lammpsDir/src/MAKE/


echo "Directory of LAMMPS is: " $lammpsDir

cd $lammpsDir/src/

# Make packages
make yes-GRANULAR
make yes-USER-NUFEB
make yes-COLLOID
make yes-USER-CFDDEM

version=`uname`
# Use different options according to different versions
if [ $version == "Linux" ]
then
    echo "The version you choose is openmpi version"
    make -j4 nufebfoam mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_nufebfoam.so .
    cd $nufebfoamDir/src/Make/
    cp options-ubuntu-openmpi options
elif [ $version == "Darwin" ]
then
    echo "The version you choose is mac version"
    make -j4 make nufebfoam mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_nufebfoammac.so .
    cd $nufebfoamDir/src/Make/
    cp options-mac-openmpi options
else
    echo "Sorry, we haven't got the required version."
fi

cd $nufebfoamDir/src/
wmake libso dragModels
wmake libso chPressureGrad
wmake libso lammpsFoamTurbulenceModels
wmake 

