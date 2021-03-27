AMREX_HOME := ../../amrex
TOP := ../../
COMP := GNU

DIM = 2

TINY_PROFILE = FALSE
PROFILE = FALSE
DEBUG = FALSE
VERBOSE = TRUE

DEFINES += -DLUA_USE_LINUX
LIBRARIES += -ldl -lreadline

CVODE_LIB_DIR=../../../sundials/instdir
USE_CVODE_LIBS=FALSE

HDF5_HOME   = /usr/lib/x86_64-linux-gnu/hdf5/openmpi
USE_HDF5 = FALSE

USE_EB = FALSE
AMREX_PARTICLES = FALSE

PYTHON_PLOT = FALSE
PYTHON_INCLUDE=/usr/include/python3.8
PYTHON_LIB=/usr/lib/python3.8/config-3.8-x86_64-linux-gnu

include ../Make.MFP

FFLAGS   += -fcheck=array-temps -fbacktrace -fbounds-check -cpp
F90FLAGS += -fcheck=array-temps -fbacktrace -fbounds-check -cpp

CXXFLAGS += -fno-omit-frame-pointer -no-pie
CFLAGS += -fno-omit-frame-pointer -no-pie

