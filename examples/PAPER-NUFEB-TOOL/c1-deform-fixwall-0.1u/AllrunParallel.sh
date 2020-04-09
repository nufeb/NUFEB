#!/bin/bash
cd ${0%/*} || exit 1 # Run from this directory

./Allclean.sh

blockMesh > log.blockMesh
decomposePar > log.decomposePar
mpirun -np 4 lammpsFoam -parallel
