#!/bin/bash -l
#SBATCH -J NUFEB_cyano1
#SBATCH -p batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --exclusive
#SBATCH --mem=100g
#SBATCH -A cnms
#SBATCH -t 10:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sakkosjo@msu.edu

module purge
module load env/cades-cnms
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load openmpi/qnfab5m
spack load vtk%gcc@8.2.0

#export OMP_NUM_THREADS=1

export LAMMPS=/lustre/or-hydra/cades-cnms/12x/NUFEB/lammps/src/lmp_mpi

ldd $LAMMPS

mpirun -np 32 $LAMMPS -in Inputscript_1.lammps
tar -zcf Run_1.tar.gz Run_1_* run_1.pkl atom_1_1.in Inputscript_1.lammps
