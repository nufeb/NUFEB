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
git checkout cyano
git pull
```
## To build

### With Ansible
Install Ansible
#### On CADES
```
module load python/3.6.3
pip3 install ansible --user
```
#### On Ubuntu
```
sudo apt update
sudo apt install software-properties-common
sudo apt-add-repository --yes --update ppa:ansible/ansible
sudo apt install ansible
```
Run the playbook
```
ansible-playbook playbook.yml
```
### Manually
```shell
module purge
module load PE-gnu/3.0
module load cmake
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
python ./tools/GenerateAtom.py --u your_user_name
salloc -A cnms -p high_mem --nodes=1 --mem=80G --exclusive -t 00:30:00
srun --ntasks-per-node 32 -n 32 ~/NUFEB/lammps/src/lmp_png -in ~/NUFEB/runs/Inputscript*.lammps
```
### batches
```shell
module load env/cades-cnms
module load anaconda3
cd $SCRATCH
python ./tools/GenerateAtom.py --u your_user_name
./tools/slurmRun.sh
```
