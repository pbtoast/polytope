#-*-makefile-*-
# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug          = not-set
MPI            = not-set
CC             = not-set
CXX            = not-set
prefix         = not-set
boost_root     = not-set
hdf5_root      = not-set
use_silo       = 1
use_python     = 0
use_C          = 1
python_exe     = not-set
python_version = not-set
build_tests    = 1
c_real_type    = double

# This proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=1 -DUNIX=1

# Process configuration options.

# Did the user specify compilers?
ifneq ($(CC), not-set)
  CC = '$(CC)'
else
  ifeq ($(MPI), 1)
    CC = mpicc
  else
    CC = cc
  endif
endif

ifneq ($(CXX), not-set)
  CXX = $(CXX)
else
  ifeq ($(MPI), 1)
    CXX = mpicxx
  else
    CXX = c++
  endif
endif

# MPI
ifeq ($(MPI), 1)
  BUILDDIR := ${BUILDDIR}-MPI
  CONFIG_FLAGS += -DUSE_MPI=1
else
  CONFIG_FLAGS += -DUSE_MPI=0
endif

CONFIG_FLAGS += -DCC='${CC}' -DCXX='${CXX}'

# Debugging symbols
ifneq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Debug
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
else
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
endif

# prefix path for installation
ifneq ($(prefix), not-set)
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX=$(prefix)
else
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX=/usr/local
endif

# explicit boost path
ifneq ($(boost_root), not-set)
  CONFIG_FLAGS += -DBOOST_ROOT=$(boost_root)
endif

# Choose to build silo or not if available
CONFIG_FLAGS += -DUSE_SILO=$(use_silo)

# Explicit HDF5 path
ifneq ($(hdf5_root), not-set)
  CONFIG_FLAGS += -DHDF5_ROOT=$(hdf5_root)
endif

# Choose to build the test set or not
CONFIG_FLAGS += -DTESTING=$(build_tests)

# Explicit path for PyBindGen
ifeq ($(use_python), 1)
  ifneq ($(python_exe), not-set)
    ifneq ($(python_version), not-set)
       CONFIG_FLAGS += -DPYTHON_EXE=$(python_exe)
       CONFIG_FLAGS += -DPYTHON_VERSION=$(python_version)
    else
       use_python = 0
    endif
  else
    use_python = 0
  endif
endif

# Choose to build python bindings with pybindgen
CONFIG_FLAGS += -DUSE_PYTHON=$(use_python)

# Choose if we're building the C interface
CONFIG_FLAGS += -DBUILD_C_INTERFACE=$(use_C)

# real number type for C library
ifeq ($(use_C), 1)
  CONFIG_FLAGS += -DC_REAL_TYPE=$(c_real_type)
endif

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else 
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

define run-config
mkdir -p $(BUILDDIR)
cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all test clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		make -C $(BUILDDIR) $@ $(MAKEFLAGS); \
	fi

config: distclean
	$(run-config)

distclean:
		rm -rf $(BUILDDIR)

#dist:
#	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall 
