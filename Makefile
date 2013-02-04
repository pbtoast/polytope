#-*-makefile-*-
# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug          = not-set
MPI            = not-set
CC             = not-set
CXX            = not-set
prefix         = not-set
boost_root     = not-set
use_silo       = 1
use_python     = 0
python_exe     = not-set
pybindgen_path = not-set

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

# Explicit path for PyBindGen
# IDEA: cmake can find the python interpreter and libraries for you. If you set use_python=1,
#       but a valid pybindgen root is not given, do nothing. 
ifeq ($(use_python), 1)
  ifneq ($(pybindgen_path), not-set)
    CONFIG_FLAGS += -DPYBINDGEN_PATH=$(pybindgen_path)
  else
    use_python = 0
  endif
  ifneq ($(python_exe), not-set)
    CONFIG_FLAGS += -DPYTHON_EXE=$(python_exe)
  endif
endif

# Choose to build python bindings with pybindgen
CONFIG_FLAGS += -DUSE_PYTHON=$(use_python)

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
