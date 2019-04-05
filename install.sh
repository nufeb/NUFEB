#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

echo "Installing NUFEB.."
currentDir=$PWD

#### Copy package and lib files to LAMMPS directory #####
echo "Copying packages to LAMMPS.."
cp -rf $currentDir/src/* $currentDir/lammps/src/
cp -rf $currentDir/lib/* $currentDir/lammps/lib/

echo "Configuring vtk-6.3 library.."
cd $currentDir/lammps/lib/nufeb
cp Makefile.lammps_vtk6.3 Makefile.lammps
cd $currentDir/lammps/lib/vtk
cp Makefile.lammps_vtk6.3 Makefile.lammps


#### Build LAMMPS with NUFEB and VTK packages#####
cd $currentDir/lammps/src
make yes-user-nufeb
make yes-user-vtk
make yes-granular

make -j4 mpi
