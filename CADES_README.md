# Instructions for NUFEB + phototroph on CADES CONDO


# clone repo
module load env/cades-cnms
cd $SCRATCH
git clone https://github.com/Jsakkos/NUFEB --recursive
git checkout add_phototroph
git pull

# To build
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load openmpi/qnfab5m
spack load vtk%gcc@8.2.0

./install_cades.sh --enable-vtk

# To test
## interactively
salloc -A cnms -p high_mem --nodes=1 --mem=80G --exclusive -t 00:30:00
srun --ntasks-per-node 32 -n 32 ../../lammps/src/lmp_mpi -in Inputscript.lammps

