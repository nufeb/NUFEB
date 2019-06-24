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
       cp Makefile.lammps_vtk8.0 Makefile.lammps
       cd ../vtk
       cp Makefile.lammps_vtk8.0 Makefile.lammps
    elif [ $var == "--enable-hdf5" ]; then
       cp Makefile.lammps_hdf5 Makefile.lammps
    elif [ $var == "--enable-vtk-hdf5" ]; then
       cp Makefile.lammps_hdf5_vtk8.0 Makefile.lammps
       cd ../vtk
       cp Makefile.lammps_vtk8.0 Makefile.lammps
    elif [ $var == "--serial" ]; then continue
    else
       echo "Unknown parameter"
       exit 1
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
	make yes-user-vtk
    fi
done

echo "Building NUFEB.."
for var in "$@"
do 
    if [ $var == "--serial" ]; then
	cd STUBS
        make
        cd ..
        make -j4 serial
        exit 1
   fi
done

make -j4 mpi
exit 1

#echo "Writing path to .bashrc"
#echo "export PATH=\$PATH:$currentDir/lammps/src/" >> ~/.bashrc
