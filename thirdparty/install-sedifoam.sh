#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing lammpsFoam (for mac/linux).."
sedifoamDir=$PWD/sedifoam
cd ..
nufebDir=$PWD
lammpsDir=$PWD/lammps
lammpsSRC=$lammpsDir/src

echo "Directory of LAMMPS is: " $lammpsDir

echo "Copying packages to LAMMPS.."
cp -rf $sedifoamDir/interfaceToLammps/MAKE $lammpsSRC/
cp -rf $sedifoamDir/interfaceToLammps/USER-CFDDEM $lammpsSRC/
cp -rf $sedifoamDir/interfaceToLammps/Makefile $lammpsSRC/
cp -rf $sedifoamDir/interfaceToLammps/lib/* $lammpsDir/lib/
cp -rf $nufebDir/src/* $lammpsSRC/
cp -rf $nufebDir/lib/* $lammpsDir/lib/

# Make STUBS 
cd $lammpsSRC/STUBS
make
cd $lammpsSRC

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
    make -j4 shanghailinux mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_shanghailinux.so .
    cd $sedifoamDir/lammpsFoam
    touch Make/options
    echo "LAMMPS_DIR ="$lammpsSRC > Make/options
    #change to 'cat Make/options-linux-openmpi >> Make/options' is you are using other linux version
    cat Make/options-ubuntu-openmpi >> Make/options
elif [ $version == "Darwin" ]
then
    echo "The version you choose is mac version"
    make -j4 shanghaimac mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_shanghaimac.so .
    cd $sedifoamDir/lammpsFoam
    touch Make/options
    echo "LAMMPS_DIR ="$lammpsSRC > Make/options
    cat Make/options-mac-openmpi >> Make/options
else
    echo "Sorry, we haven't got the required version."
    echo "Please contact the developer (sunrui@vt.edu) for help."
fi

wmake libso dragModels
wmake libso chPressureGrad
wmake libso lammpsFoamTurbulenceModels
wmake 

