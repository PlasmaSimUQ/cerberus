#!/bin/bash

# load modules
module purge
module load gnu
module load openmpi_ib


# now build
#make clean
make -j8
