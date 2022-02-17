#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing lammpsFoam (for mac/linux).."
currentDir=$PWD
echo "Enter the directory of your LAMMPS and press [ENTER] "
echo -n "(default directory ./lammps-5Nov16): "
read lammpsDir

# Determine if the directory of LAMMPS exists or not.
# If not, look for LAMMPS in the default directory.
if [ ! -d "$lammpsDir" ]
then
    echo "Directory NOT found! Use default directory instead."
    lammpsDir="$PWD/lammps-5Nov16"
fi

cd $lammpsDir
lammpsDir=$PWD

echo "Directory of LAMMPS is: " $lammpsDir

# Copy/link all the extra implementations
cd $lammpsDir/src
lammpsSRC=$PWD

echo "Copying packages to LAMMPS.."
cp -rf $currentDir/interfaceToLammps/MAKE $lammpsSRC/
cp -rf $currentDir/interfaceToLammps/USER-CFDDEM $lammpsSRC/
cp -rf $currentDir/interfaceToLammps/Makefile $lammpsSRC/
cp -rf $currentDir/interfaceToLammps/lib/* $lammpsDir/lib/

# Make STUBS 
cd $lammpsDir/src/STUBS
make
cd $lammpsDir/src

# Make packages
make yes-GRANULAR
make yes-KSPACE
make yes-MANYBODY
make yes-MOLECULE
make yes-FLD # lubrication
make yes-RIGID # freeze
make yes-MISC # deposit
make yes-VORONOI # ??

version=`uname`
# Use different options according to different versions
if [ $version == "Linux" ]
then
    echo "The version you choose is openmpi version"
    make -j4 shanghailinux mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_shanghailinux.so .
    cd $currentDir/lammpsFoam
    touch Make/options
    echo "LAMMPS_DIR ="$lammpsSRC > Make/options
    cat Make/options-linux-openmpi >> Make/options
elif [ $version == "Darwin" ]
then
    echo "The version you choose is mac version"
    make -j4 shanghaimac mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_shanghaimac.so .
    cd $currentDir/lammpsFoam
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

