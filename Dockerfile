FROM ubuntu:20.04
RUN useradd --create-home --shell /bin/bash admin

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    cmake \
    git \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    libpng-dev \
    python3-pip \
    bash \
    nano

USER admin
WORKDIR /home/admin
RUN git clone https://github.com/Jsakkos/NUFEB nufeb --recursive

WORKDIR /home/admin/nufeb
RUN git checkout cyano

WORKDIR /home/admin/nufeb/thirdparty
RUN chmod +x ./install-hdf5.sh
RUN ./install-hdf5.sh

WORKDIR /home/admin/nufeb
RUN chmod +x ./install.sh
RUN ./install.sh --enable-hdf5
RUN pip install nufeb-tools -U --user
ENV PATH $HOME/.local/bin:$PATH"

