#!/bin/bash

mpirun -np 1 ../../../lammps/src/lmp_mpi -in Inputscript.lammps > Inputscript_50_50_100.log
cd Run_50_50_100_1
tar -zcf VTK.tar.gz *.vtr *.vtu *.vti
cd ..
