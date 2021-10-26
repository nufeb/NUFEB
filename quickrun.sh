#!/bin/bash -l
export LAMMPS=~/NUFEB/lammps/src/lmp_png

base=$PWD

#run NUFEB simulations
for dir in runs/*/
do
cd "$dir"
mpirun -np 6 $LAMMPS -in *.lammps > nufeb.log
cd "$base"
done
