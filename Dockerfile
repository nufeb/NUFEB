FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    libpng-dev \
    libvtk6-dev \
    python3-pip
ADD . /nufeb
RUN useradd --create-home --shell /bin/bash admin
WORKDIR /nufeb/thirdparty
RUN chmod +x ./install-hdf5.sh
RUN ./install-hdf5.sh
RUN chmod +x ./install-vtk.sh
RUN ./install-vtk.sh
WORKDIR /nufeb
RUN chmod +x ./install.sh
RUN ./install.sh --enable-hdf5 --enable-vtk
USER admin
RUN pip install nufeb-tools -U