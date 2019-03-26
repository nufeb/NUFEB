FROM ubuntu:16.04

#add new user
RUN useradd --create-home --shell /bin/bash admin

#install required libs
RUN apt-get update
RUN apt-get -y install build-essential sudo make git openmpi-bin openmpi-common libopenmpi-dev libvtk6-dev vim 

#switch to user admin, download NUFEB
USER admin
WORKDIR /home/admin
RUN git clone https://nufebdocker:nufeb123@github.com/shelllbw/nufeb.git
WORKDIR /home/admin/nufeb

#setup vtk makefile
WORKDIR /home/admin/nufeb/code/lammps5Nov16/lib/vtk
RUN cp Makefile.lammps.ubuntu16.04_vtk6.2 Makefile.lammps
WORKDIR /home/admin/nufeb/code/lammps5Nov16/lib/nufeb
RUN cp Makefile.lammps.ubuntu16.04_vtk6.2 Makefile.lammps

#compile NUFEB code
WORKDIR /home/admin/nufeb/code/lammps5Nov16/src
RUN make yes-user-nufeb
RUN make yes-granular
RUN make yes-user-vtk
RUN make -j4 mpi

#add NUFEB executable to bashrc
RUN echo "export PATH=$PATH:/home/admin/nufeb/code/lammps5Nov16/src/" >> /home/admin/.bashrc

