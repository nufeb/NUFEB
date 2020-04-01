#!/bin/bash
#SBATCH -A tcnufeb
#SBATCH --ntasks=100
#SBATCH -c 1
#SBATCH -N 1-3
module load intel

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PATH=$PATH:/mnt/nfs/home/nbl21/nufeb/code/lammps5Nov16/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/nfs/home/nbl21/nufeb/code/lammps5Nov16/lib:/mnt/nfs/home/nbl21/local/lib

srun lmp_rocket_intel -in Inputscript.lammps 
#srun lmp_g++_openmpi -in Inputscript.lammps
