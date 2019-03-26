# NUFEB
NUFEB is an open source tool for Individual-based Model (IbM) simulation. 
The tool is based on the Individual-based Modelling (IbM) approach, 
where microbes are represented as discrete units and their
behaviour changes over time due to a variety of processes. 

NUFEB is built on top of the classical molecular dynamics simulator
LAMMPS, extended with IbM features. A wide range of biological, physical and
chemical processes are implemented to explicitly model microbial systems. 
NUFEB is fully parallelised and allows for the simulation of large numbers of microbes 
(10^7 microbes and beyond).

NUFEB is a freely-available open-source code, distributed under the terms
of the GNU Public License.

NUFEB development has been funded by the UK’s EPSRC EP/K039083/1 
Newcastle University Frontiers in Engineering Biology (NUFEB) project.

### Building

NUFEB requires GCC/G++ and OpenMPI libraries for a successful build. 

To compile this code, go to the lammps5Nov16/src directory:

$ cd NUFEB/lammps5Nov16/src/

Then, install the NUFEB and granular packages in /src directory with the following instruction:

$ make yes-USER-NUFEB

$ make yes-GRANULAR

Finally, execute the following command to compile the NUFEB executable:

$ make mpi

### Running

You can run NUFEB by going to one of the sub-directories in /examples/ and run:

$  mpirun -np 4 PATH_TO_SRC/./lmp_mpi -in Inputscript.lammps
