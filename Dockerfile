FROM ubuntu:20.04
MAINTAINER Jonathan Sakkos <sakkosjo@msu.edu>
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    cmake \
    git-core \
    g++ \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \ 
    libpng-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    mesa-common-dev \
ADD /nufeb /nufeb
WORKDIR nufeb/thirdparty
RUN ./install-hdf5.sh && ./install-vtk.sh
WORKDIR ../
RUN ./install.sh --enable-vtk --enable-hdf5
