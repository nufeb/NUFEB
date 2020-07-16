# Instructions for NUFEB + phototroph on CADES CONDO

# To build

module load env/cades-cnms
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load openmpi/qnfab5m
spack load vtk%gcc@8.2.0

./install_cades.sh

# To test
## interactively
salloc -A cnms -p high_mem --nodes=1 --mem=80G --exclusive -t 00:30:00
srun --ntasks-per-node 32 -n 32 ../../lammps/src/lmp_mpi -in Inputscript.lammps

