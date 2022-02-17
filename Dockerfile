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
ADD /nufeb /nufeb
#RUN chmod 755 -R /nufeb
#WORKDIR /nufeb/thirdparty
#RUN chmod 755 install-hdf5.sh
#RUN ./install-hdf5.sh 
#RUN chmod 755 install-vtk.sh
#RUN ./install-vtk.sh
RUN sed -i 's/\r$//' nufeb/install.sh  && \  
    chmod +x nufeb/install.sh
WORKDIR /nufeb
RUN ./install.sh

#ENTRYPOINT ./install.sh
#RUN chmod +x ./install.sh


