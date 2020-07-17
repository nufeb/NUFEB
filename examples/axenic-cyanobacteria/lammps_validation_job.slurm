#!/bin/bash -l
#SBATCH -J sample_job
#SBATCH -p high_mem
#SBATCH --nodes 2
#SBATCH --ntasks-per-node 32
#SBATCH --cpus-per-task 1
#SBATCH --exclusive
#SBATCH --mem=100g
#SBATCH -A cnms
#SBATCH -t 04:00:00

#launch from a fresh shell with the following loaded

#module load env/cades-cnms
#. $SOFTWARECNMS/spack/share/spack/setup-env.sh
#spack load openmpi/qnfab5m
#spack load vtk%gcc@8.2.0

# The job shell will inherit this, somehow modules are messed up
# in the slurm script shell so you don't get what you
# want if you load here in the script.

export OMP_NUM_THREADS=1


export LAMMPS=/lustre/or-hydra/cades-cnms/12x/NUFEB/lammps/src/lmp_mpi
module list
printenv
ldd $LAMMPS

echo "srun -n 8 $LAMMPS -in Inputscript-validation.lammps"

srun -n 8 $LAMMPS -in Inputscript-validation.lammps
