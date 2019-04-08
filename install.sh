#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

echo "Installing NUFEB.."
currentDir=$PWD

#### Copy package and lib files to LAMMPS directory #####
echo "Copying packages to LAMMPS.."
cp -rf $currentDir/src/* $currentDir/lammps/src/
cp -rf $currentDir/lib/* $currentDir/lammps/lib/

echo "Configuring Makefile.lammps.."

cd $currentDir/lammps/lib/nufeb

for var in "$@"
do 
    if [ $var == "--enable-vtk" ] ; then
       cp Makefile.lammps_vtk6.3 Makefile.lammps
       cd ../vtk
       cp Makefile.lammps_vtk6.3 Makefile.lammps
    elif [ $var == "--enable-hdf5" ]; then
       cp Makefile.lammps_hdf5 Makefile.lammps
    elif [ $var == "--enable-essential" ]; then
       cp Makefile.lammps_essential Makefile.lammps
    elif [ $var == "--enable-vtk-hdf5" ]; then
       cp Makefile.lammps_hdf5_vtk6.3 Makefile.lammps
       cd ../vtk
       cp Makefile.lammps_vtk6.3 Makefile.lammps
    else
       echo "Unknown parameter"
    fi
done


#### Build LAMMPS with NUFEB and VTK packages#####
echo "Installing required packages.."

cd $currentDir/lammps/src
make yes-user-nufeb
make yes-granular

for var in "$@"
do 
    if [ $var == "--enable-vtk" ] || [ $var == "--enable-vtk-hdf5" ]; then
	echo "export LD_LIBRARY_PATH=$currentDir/thirdparty/vtk-build/vtk-6.3/" >> ~/.bashrc
	make yes-user-vtk
    fi
    if [ $var == "--enable-hdf5" ] || [ $var == "--enable-vtk-hdf5" ]; then
	echo "export LD_LIBRARY_PATH=$currentDir/thirdparty/hdf5/hdf5/lib/" >> ~/.bashrc
    fi
done

make -j4 mpi

