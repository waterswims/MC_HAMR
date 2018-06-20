FROM ubuntu:latest

RUN apt-get update
RUN apt-get install -y git ssh tar gzip ca-certificates gcc g++ gfortran make cmake zlib1g-dev libgsl-dev

RUN wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz
RUN tar -xvf openmpi-3.0.0.tar.gz

RUN cd openmpi-3.0.0; ./configure --prefix=/opt/openmpi/ CC=gcc CXX=g++ F77=gfortran FC=gfortran
RUN cd openmpi-3.0.0; make all install

RUN wget ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-*.tar.gz
RUN tar -xvf hdf5-*.tar.gz
RUN cd hdf5-*; CC=/opt/openmpi/bin/mpicc ./configure --enable-parallel --prefix=/opt/hdf5/
RUN cd hdf5-*; make
RUN cd hdf5-*; make install
