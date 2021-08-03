#!/bin/bash -l
export LAMMPS=~/NUFEB/lammps/src/lmp_png

ldd $LAMMPS

base=$PWD

#run NUFEB simulations
for dir in runs/*/
do
cd "$dir"
mpirun -np 8 $LAMMPS -in *.lammps > nufeb.log
cd "$base"
done

#check if the previous run went ok, exit if not
if [ $? -ne 0 ]
then
    echo "Something went wrong while running simulations, exiting"
    exit
fi

date

#do the post-processing tasks here

#create tarballs for the VTK files
for dir in runs/*/
do
cd "$dir"
tar -zcf VTK.tar.gz *.vtr *.vtu *.vti
rm *.vtr *.vtu *.vti
cd "$base"
done

#check if the previous run went ok, exit if not
if [ $? -ne 0 ]
then
    echo "Something went wrong while creating tarballs, exiting"
    exit
fi

