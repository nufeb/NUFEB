#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing SediFOAM (for mac/linux).."

echo "*******************************************"
echo "select the system you are running, then press enter"
echo "  1) Ubuntu14.x - Ubuntu16.x"
echo "  2) Ubuntu17.x - Ubuntu20.x" 
echo "  3) Centos"
echo "  4) Mac" 
echo "  5) Other (install openmpi locally, this will take a while)" 
echo "*******************************************"
read n

case $n in
  1) echo "You chose 1) Ubuntu14.x - Ubuntu16.x";;
  2) echo "You chose 2) Ubuntu17.x - Ubuntu20.x";;
  3) echo "You chose 3) Centos";;
  4) echo "You chose 4) Mac";;
  5) echo "You chose 5) Other";;
  *) echo "Unknown option"; exit;;
esac

sedifoamDir=$PWD/sedifoam
cd ..
ompiDir=$PWD/thirdparty/openmpi-3.0.6/ompi-build
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
make yes-USER-VTK

# Use different options according to different versions
if [ $n == 4 ]
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
    make -j4 shanghailinux mode=shlib
    cd $FOAM_USER_LIBBIN
    ln -sf $lammpsDir/src/liblammps_shanghailinux.so .
    cd $sedifoamDir/lammpsFoam
    touch Make/options
    echo "LAMMPS_DIR ="$lammpsSRC > Make/options

    if [ $n == 1 ] 
    then 
	cat Make/options-ubuntu16-openmpi >> Make/options
    elif [ $n == 2 ] 
    then 
	cat Make/options-ubuntu18-openmpi >> Make/options
    elif [ $n == 3 ] 
    then 
	cat Make/options-linux-openmpi >> Make/options
    elif [ $n == 5 ] 
    then 
        sh $nufebDir/thirdparty/./install-openmpi.sh
        echo "OMPI_DIR ="$ompiDir >> Make/options
	cat Make/options-local-openmpi >> Make/options
    fi
fi 

wmake libso dragModels 
wmake libso chPressureGrad 
wmake libso lammpsFoamTurbulenceModels
wmake 

