# Instructions for NUFEB + phototroph on CADES CONDO


# clone repo
module load env/cades-cnms
```shell
cd $SCRATCH
git clone https://github.com/Jsakkos/NUFEB --recursive
git checkout add_phototroph
git pull
```
# To build
```shell
. $SOFTWARECNMS/spack/share/spack/setup-env.sh
spack load openmpi/qnfab5m
spack load vtk%gcc@8.2.0

./install_cades.sh --enable-vtk
```
# To test
## interactively
```shell
salloc -A cnms -p high_mem --nodes=1 --mem=80G --exclusive -t 00:30:00
srun --ntasks-per-node 32 -n 32 ../../lammps/src/lmp_mpi -in Inputscript.lammps
```
## batches
```shell
module load env/cades-cnms
module load anaconda3
cd $SCRATCH
cd NUFEB/examples/cyanobacteria-sucrose
python NUFEBatom --n 3 --r 3
./slurmRun.sh
```
# To login remotely
https://support.cades.ornl.gov/user-documentation/_book/external-access-ucams.html
```shell
ssh ucams@login1.ornl.gov
ssh username@mod-condo-login.ornl.gov
```
