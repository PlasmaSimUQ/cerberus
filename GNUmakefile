# Initialize the amrex submodule (if not already done)
$(info "Initializing/updating the amrex submodule")
$(shell git submodule update --init --recursive)

TOP := $(realpath .)

AMREX_HOME := $(TOP)/amrex
COMP := GNU

DIM = 2

TINY_PROFILE = FALSE
PROFILE = FALSE
DEBUG = FALSE
VERBOSE = TRUE

DEFINES += -DLUA_USE_LINUX
LIBRARIES += -ldl -lreadline

USE_EB = FALSE
AMREX_PARTICLES = FALSE

PYTHON_PLOT = FALSE
PYTHON_INCLUDE = /usr/include/python3.10/
PYTHON_LIB = /usr/lib/python3.10/config-3.10-x86_64-linux-gnu

PPROF = FALSE
PPROF_INCLUDE = /usr/include/gperftools
PPROF_LIB = /usr/lib/x86_64-linux-gnu

EILMER_GAS = TRUE
DLANG_LIB := /usr/lib

EILMER_HOME = $(TOP)/gdtkinst

# this url will be sent to the "eilmer" make rule in cerberus.rules
EILMER_URL = https://github.com/gdtk-uq/gdtk

include $(TOP)/cerberus.rules

FFLAGS   += -fcheck=array-temps -fbacktrace -fbounds-check -cpp
F90FLAGS += -fcheck=array-temps -fbacktrace -fbounds-check -cpp

CXXFLAGS += -fno-omit-frame-pointer -no-pie
CFLAGS += -fno-omit-frame-pointer -no-pie
