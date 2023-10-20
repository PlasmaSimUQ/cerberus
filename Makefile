# Colors for pretty makefile output
GREEN := $(shell echo -e "\033[0;32m")
YELLOW := $(shell echo -e "\033[0;33m")
RED := $(shell echo -e "\033[0;31m")
END := $(shell echo -e "\033[0m")

# set the root path
TOP := $(shell git rev-parse --show-toplevel)
INSTALLDIR ?= ./MFP
BUILDDIR ?= ./tmp_build_dir


# include the local makefile (ignore if it doesn't exist)
-include $(TOP)/Make.local

DEBUG ?= FALSE

CERBERUS_GIT_NAME := $(shell git branch --show-current)
CERBERUS_GIT_VERSION := $(shell git describe --abbrev --dirty --always --tags)

CHECK_UPDATES ?= FALSE

AMREX_HOME := $(TOP)/amrex
UPDATE_AMREX ?= TRUE

ifeq ($(UPDATE_AMREX), TRUE) # and
ifneq ($(BUILD_ACTUAL), TRUE)
  $(info ${GREEN}Initialising/updating the AMReX submodule${END})
  $(shell git submodule update --init --recursive -- $(AMREX_HOME))
endif
endif


BinDir := ./
BASE := MFP
EBASE := $(BinDir)/$(BASE).

$(info $(f_DIM_CHECK))

ifndef DIM
  fc_DIM_NO_USER = TRUE
  DIM ?= 2
endif

ifneq ($(DIM), 1)
ifneq ($(DIM), 2)
ifneq ($(DIM), 3)
  $(error "${RED}DIM must be one of: [1, 2, 3]!${END}")
endif
endif
endif

COMP ?= GNU

VERBOSE ?= FALSE

LAZY ?= TRUE

USE_MPI ?= TRUE

# force this to FALSE as currently unsupported
ifdef USE_OMP
  fc_OMP_ENABLED = TRUE
endif

override USE_OMP = FALSE

WARN_ALL ?= FALSE

#=== EB ===

ifndef USE_EB
  ifeq ($(DIM), 1)
    fc_EB_NO_USER_WRONG_DIM=TRUE
    USE_EB = FALSE
  else
    fc_EB_NO_USER=TRUE
    USE_EB = TRUE
endif
endif

USE_EB ?= TRUE
ifeq ($(USE_EB), TRUE)
  ifeq ($(DIM), 1)
    fc_EB_WRONG_DIM=TRUE
    override USE_EB = FALSE
  else
    USERSuffix := $(USERSuffix).EB
  endif
endif


#=== AMREX PARTICLES ===

ifndef AMREX_PARTICLES
  ifeq ($(DIM), 1)
    fc_PARTICLES_NO_USER_WRONG_DIM=TRUE
    AMREX_PARTICLES = FALSE
  else
    fc_PARTICLES_NO_USER=TRUE
    AMREX_PARTICLES = TRUE
  endif
endif

AMREX_PARTICLES ?= TRUE
ifeq ($(AMREX_PARTICLES), TRUE)
  ifeq ($(DIM), 1)
    fc_PARTICLES_WRONG_DIM=TRUE
    override AMREX_PARTICLES = FALSE
  else
    DEFINES += -DAMREX_PARTICLES
    USERSuffix := $(USERSuffix).PARTICLES
  endif
endif


#=== EILMER GAS ===

ifndef EILMER_GAS
  fc_EILMER_NO_USER=TRUE
  EILMER_GAS ?= FALSE
endif

EILMER_SRC ?= $(TOP)/gdtk_git
EILMER_URL ?= https://github.com/gdtk-uq/gdtk
EILMER_HOME ?= $(TOP)/gdtkinst

ifeq ($(DEBUG), TRUE)
  EILMER_COMPILE_FLAVOUR = debug
else
  EILMER_COMPILE_FLAVOUR = fast
endif


DLANG_LIB ?= /usr/lib

# ORIGIN is the directory where the MFP executable is located
ifeq ($(EILMER_GAS), TRUE)
  LIBRARY_LOCATIONS += $(DLANG_LIB)
  LIBRARY_LOCATIONS += $(EILMER_HOME)/lib
  LIBRARIES += -llua
  LIBRARIES += -lgas
  LIBRARIES += -Wl,-rpath='$$ORIGIN/../lib'
  LIBRARIES += -Wl,-rpath='$$ORIGIN/lib'
  LIBRARIES += -Wl,-rpath='$$ORIGIN'
  LIBRARIES += -Wl,-rpath=$(EILMER_HOME)/lib
  USERSuffix := $(USERSuffix).EILMER
  DEFINES += -DEILMER_GAS
else
  # if eilmer doesn't exist, then use our lua lib
  include $(TOP)/Source/extern/lua-5.4.3/Make.package
endif

#=== PPROF ===

TINY_PROFILE ?= FALSE

PPROF ?= FALSE
PPROF_INCLUDE ?= /usr/include/gperftools
PPROF_LIB ?= /usr/lib/x86_64-linux-gnu

ifeq ($(PPROF), TRUE)
  INCLUDE_LOCATIONS += $(PPROF_INCLUDE)
  LIBRARY_LOCATIONS += $(PPROF_LIB)
  LIBRARIES         += -lprofiler
  USERSuffix := $(USERSuffix).PPROF
  DEFINES += -DPPROF
endif

#=== PYTHON ===

PYTHON_PLOT ?= FALSE
PYTHON_INCLUDE ?= /usr/include/python3.10 #/usr/include/python3.8
PYTHON_LIB ?= /usr/lib/python3.10/config-3.10-x86_64-linux-gnu

#/usr/lib/python3.8/config-3.8-x86_64-linux-gnu

ifeq ($(PYTHON_PLOT), TRUE)
  INCLUDE_LOCATIONS += $(PYTHON_INCLUDE)
  LIBRARY_LOCATIONS += $(PYTHON_LIB)
  LIBRARIES         += -l$(notdir $(PYTHON_INCLUDE))
  USERSuffix := $(USERSuffix).PYTHON
  DEFINES += -DPYTHON
endif


#=== SYMPLECTIC INTERPOLANT ===
SYMPLECTIC ?= NONE

ifeq ($(SYMPLECTIC),P8R2)
  DEFINES += -DINTERPOLATION_P8R2 -DWRANGE=2 -DSYMPLECTIC
  USERSuffix := $(USERSuffix).SYMPLECTIC_P8R2
else ifeq ($(SYMPLECTIC),PWL)
  DEFINES += -DINTERPOLATION_PWL -DWRANGE=1 -DSYMPLECTIC
  USERSuffix := $(USERSuffix).SYMPLECTIC_PWL
else ifeq ($(SYMPLECTIC), NONE)
 # do nothing
else
  $(error ${RED}(experimental) SYMPLECTIC=$(SYMPLECTIC) is not valid! Must be one of: [PR82, PWL, NONE]${END})
endif

# Disable Fortran
BL_NO_FORT = TRUE

DEFINES += -DLUA_USE_LINUX

LIBRARIES += -ldl -lreadline

CERBERUS_GIT_NAME := $(shell git branch --show-current)
CERBERUS_GIT_VERSION := $(shell git describe --abbrev --dirty --always --tags)

DEFINES += -DCERBERUS_GIT_VERSION=\"$(CERBERUS_GIT_VERSION)\"
DEFINES += -DCERBERUS_GIT_NAME=\"$(CERBERUS_GIT_NAME)\"


# ccache
USE_CCACHE ?= TRUE
ifeq ($(USE_CCACHE), TRUE)
  ifeq (, $(shell command -v ccache))
    CCACHE_EXEC_NOT_FOUND=TRUE
    override USE_CCACHE = FALSE
  endif
endif

CXXFLAGS += -fno-omit-frame-pointer -no-pie
CFLAGS += -fno-omit-frame-pointer -no-pie

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(TOP)/Source/Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Amr/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package

ifeq ($(USE_EB), TRUE)
  include $(AMREX_HOME)/Src/EB/Make.package
endif

ifeq ($(AMREX_PARTICLES), TRUE)
  include $(AMREX_HOME)/Src/Particle/Make.package
endif

# possibly override the default amrex "tmp_build_dir" directory with our own
TmpBuildDir := $(BUILDDIR)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
include $(TOP)/cerberus.rules
