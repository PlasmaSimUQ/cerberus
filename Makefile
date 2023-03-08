# Initialize the amrex submodule (if not already done)
$(info "Initializing/updating submodules")
$(shell git submodule update --init --recursive)

TOP := $(shell git rev-parse --show-toplevel)

# include the local makefile (ignore if it doesn't exist)
-include $(TOP)/Make.local


AMREX_HOME ?= $(TOP)/amrex
COMP ?= GNU

# use ccache if available
USE_CCACHE ?= TRUE

DIM ?= 2

TINY_PROFILE ?= FALSE
PROFILE ?= FALSE
DEBUG ?= FALSE
VERBOSE ?= FALSE
LAZY ?= TRUE
USE_MPI ?= TRUE
USE_OMP ?= FALSE
WARN_ALL ?= FALSE

# disable fortran
BL_NO_FORT=TRUE

DEFINES += -DLUA_USE_LINUX
LIBRARIES += -ldl -lreadline

USE_EB ?= TRUE
AMREX_PARTICLES ?= TRUE

PYTHON_PLOT ?= FALSE
PYTHON_INCLUDE ?= /usr/include/python3.10/
PYTHON_LIB ?= /usr/lib/python3.10/config-3.10-x86_64-linux-gnu

PPROF ?= FALSE
PPROF_INCLUDE ?= /usr/include/gperftools
PPROF_LIB ?= /usr/lib/x86_64-linux-gnu

EILMER_GAS ?= FALSE
DLANG_LIB ?= /usr/lib

EILMER_HOME ?= $(TOP)/gdtkinst
EILMER_LIB ?= $(EILMER_HOME)/lib

SYMPLECTIC ?= FALSE

CXXFLAGS += -fno-omit-frame-pointer -no-pie
CFLAGS += -fno-omit-frame-pointer -no-pie

include $(TOP)/cerberus.rules

