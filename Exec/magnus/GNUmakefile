AMREX_HOME := ../../../amrex
TOP := ../../
COMP := GNU

DIM = 2

USE_EB = FALSE

NO_MPI_CHECKING = TRUE

TINY_PROFILE = FALSE
PROFILE = FALSE
DEBUG = FALSE
VERBOSE = TRUE

USE_HDF5 = TRUE
HDF5_HOME   = /opt/cray/pe/hdf5-parallel/default/gnu/5.1
INCLUDE_LOCATIONS += ${HDF5_HOME}/include
VPATH_LOCATIONS   += ${HDF5_HOME}/include

USE_EB = FALSE
AMREX_PARTICLES = FALSE

include ../Make.MFP

CXX := CC
CC  := cc
FC  := ftn
F90 := ftn

