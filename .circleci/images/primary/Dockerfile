FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y git ssh tar gzip ca-certificates gcc g++ gfortran make cmake zlib1g-dev libgsl-dev openmpi-bin openmpi-common openssh-client openssh-server libopenmpi-dev

RUN wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
RUN tar -xvf hdf5-*.tar.gz
RUN cd hdf5-*; CC=mpicc ./configure --enable-parallel --prefix=/opt/hdf5/
RUN cd hdf5-*; make
RUN cd hdf5-*; make install
