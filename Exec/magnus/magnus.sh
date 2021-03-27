#!/bin/bash

# load modules
#module swap PrgEnv-cray PrgEnv-intel
module swap PrgEnv-cray PrgEnv-gnu
#module swap gcc gcc/4.9.3
module load boost
module load cray-hdf5-parallel

export CRAY_CPU_TARGET=haswell

# now build
#make clean
make -j8
