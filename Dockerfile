FROM ubuntu:20.04
#MAINTAINER Jonathan Sakkos <sakkosjo@msu.edu>
SHELL ["/bin/bash", "-c"]
#  && apt-get upgrade -y 
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \ 
    libpng-dev \
    libvtk6-dev \
    python3-pip
#WORKDIR /NUFEB
RUN pip install nufeb-tools -U
RUN ls
RUN pwd
ADD . /nufeb
RUN echo "Installing NUFEB.."
RUN echo "Copying packages to LAMMPS.."
ADD ./src /nufeb/lammps/src
ADD ./lib /nufeb/lammps/lib
RUN echo "Configuring Makefile.lammps.."
ADD ./lib/nufeb/Makefile.lammps /nufeb/lammps/lib/nufeb/Makefile.lammps
RUN echo "Installing required packages.."
WORKDIR /nufeb/lammps/src
RUN make yes-user-nufeb
RUN make yes-granular
#RUN make yes-user-vtk
RUN make -j4 mpi
#RUN chmod 755 -R /nufeb
#WORKDIR /nufeb/thirdparty
#RUN chmod +x ./install-hdf5.sh
#RUN ./install-hdf5.sh
#RUN chmod 755 install-hdf5.sh
#RUN ./install-hdf5.sh 
#RUN chmod 755 install-vtk.sh
#RUN ./install-vtk.sh
#ENTRYPOINT ./install.sh
#RUN chmod +x ./install.sh


