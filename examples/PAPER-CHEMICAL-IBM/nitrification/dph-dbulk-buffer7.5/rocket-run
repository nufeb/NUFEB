#!/bin/bash
#SBATCH -A tcnufeb
#SBATCH --ntasks=60
#SBATCH -c 1
#SBATCH -N 1-2
module load intel
module load HDF5/1.10.1-intel-2017.03-GCC-6.3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PATH=$PATH:/mnt/nfs/home/nbl21/nufeb/code/lammps5Nov16/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/nfs/home/nbl21/nufeb/code/lammps5Nov16/lib:/mnt/nfs/home/nbl21/local/lib

srun lmp_rocket_intel -in Inputscript.lammps 
