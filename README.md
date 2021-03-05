Fork of [NUFEB](https://github.com/nufeb/NUFEB) to simulation sucrose-secreting cyanobacteria with heterotrophic partners.

# Instructions for compiling and testing NUFEB phototroph on CADES

## To login to CADES remotely
[See CADES Documentation](https://docs.cades.ornl.gov/#external-access-ucams/)
```shell
ssh userID@login1.ornl.gov
ssh userID@or-slurm-login.ornl.gov
```

## clone repo into your home directory
```shell
git clone https://github.com/Jsakkos/NUFEB --recursive
git checkout add_phototroph
git pull
```
## To build
```shell
module purge
module load PE-gnu/3.0
cd ~/NUFEB/thirdparty/
./install-hdf5.sh
cd ~/NUFEB/./install.sh --enable-hdf5
./install_cades.sh --enable-vtk
```
## To test

### interactively
```shell
module purge
module load env/cades-cnms
module load PE-gnu/3.0
module load anaconda3
cd $SCRATCH
python NUFEBatom.py
salloc -A cnms -p high_mem --nodes=1 --mem=80G --exclusive -t 00:30:00
srun --ntasks-per-node 32 -n 32 ~/NUFEB/lammps/src/lmp_png -in ~/NUFEB/examples/cyanobacteria-sucrose/Inputscript_1.lammps
```
### batches
```shell
module load env/cades-cnms
module load anaconda3
cd $SCRATCH
cd ~/NUFEB/examples/cyanobacteria-sucrose
python NUFEBatom.py --n 3 --r 3
./slurmRun.sh
```

