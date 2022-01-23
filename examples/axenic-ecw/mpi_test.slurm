#!/bin/bash -l
#SBATCH -J ecw_mpi_test
#SBATCH -p testing
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --exclusive
#SBATCH --mem=50g
#SBATCH -A cnms
#SBATCH -t 00:30:00

module purge
module load env/cades-cnms
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load openmpi/qnfab5m
spack load vtk%gcc@8.2.0

#export OMP_NUM_THREADS=1

export LAMMPS=/lustre/or-hydra/cades-cnms/12x/NUFEB/lammps/src/lmp_mpi

ldd $LAMMPS

mpirun -np 16 $LAMMPS -in Inputscript.lammps
