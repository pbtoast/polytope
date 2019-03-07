[![License: BSD](https://img.shields.io/badge/License-BSD%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

# Polytope

Copyright (c) 2013, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory
Written by Mike Owen, David Starinshak, and Jeffrey Johnson.
LLNL-CODE-647432
All rights reserved.

Polytope is a C/C++ library for generating polygonal and polyhedral meshes.
It makes use of various 2D and 3D tessellation techniques, but provides a
single representation for these tessellations, and a single interface for
generating them.

Polytope includes bindings for Python, so you can write Python scripts to
easily generate your tessellations, and incorporate them into your own
mesh generation tools.

## Installation

Polytope works on most Linux and Mac systems.

### Software Requirements

+ A C++11 compiler.
+ MPI for parallelism
+ The Bourne Again SHell (bash), for bootstrapping
+ CMake 3.10+, for configuring and generating build files
+ GNU Make or Ninja, for performing the actual build
+ A recent version of Perl, for HDF5's installation process

If you want to build Python bindings, you also need the following:
+ A [Python 2](https://www.python.org/downloads) interpreter.
+ The [pybind11](https://github.com/pybind/pybind11) Python/C++11
  interoperability layer
+ The [PYB11Generator](https://github.com/jmikeowen/PYB11Generator) tool,
  a Python code generator that generates input for pybind11.

### Building

To build polytope on a UNIX-like system, change to your `polytope` source
directory and type the following commands:

```
./bootstrap build_dir
```

where `build_dir` is the directory in which you want to build. Then just
follow the onscreen directions: you change to that build directory, edit
`config.sh` to define your build, and then start the build using your
generator's build process. For the default generator (UNIX makefiles), this
is just `make`. For Ninja (recommended if you have it), it's `ninja`.

### Installing

To install polymec, use the install command for the generator you've selected.
For example, if you're using a generator that writes UNIX makefiles, run

```
make install [-j #threads]
```

from your build directory.

### Other Targets

These targets all work with Make and Ninja.

+ `test` - Runs all unit tests for the library. Use `ctest -j #threads` instead, though, to run the tests in parallel.
+ `clean` - Removes all build assets but retains configuration options.
+ `distclean` - Performs clean and completely removes the build directory.

## Other Considerations

Polytope provides interfaces for a number of geometry-related tools:

+ Voro++ by Chris Rycroft at LBL
+ Triangle by Jonathan Shewchuck at Berkeley
+ Tetgen by Hang Si
+ Boost.Polygon.Voronoi, part of the Boost C++ Library

Support for Triangle and Tetgen is available with minimal effort, but you
must respect the licenses of their authors and make arrangements with them
if you want to use them in a commercial application. As a result, we don't
distribute the source for either of these tools. If you want to build in
support for Triangle and/or Tetgen, you will need to do the following:

+ For Triangle, copy `triangle.h` and `triangle.c` to `src/`.
+ For Tetgen, copy `tetgen.h` and `tetgen.cxx` to `src/`.

