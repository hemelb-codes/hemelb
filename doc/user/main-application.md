# HemeLB Users' Guide

A simple guide to compile and run HemeLB

## Dependencies

The main HemeLB application requires:

- Unix (-ish) OS: modern Linux or MacOS
- C++17 compiler
- CMake version >= 3.13
- MPI implementation supporting version 3.0 of the standard or later
- Boost* header-only libraries version >= 1.54
- CTemplate*
- Catch2* (if building tests)
- ParMETIS* graph partitioning library
- TinyXML* XML parsing library
- zlib* compression library (likely included in your OS)
- MPWide* (if building for multiscale)

If using resolved blood cells, you also require:
- HDF5*
- VTK*>=9 (for RBC mode)

The dependencies marked with * can be built automatically by the build system.

## Download

Clone the repository:
```bash
git clone https://github.com/hemelb-codes/hemelb.git
```

## Compile

TL;DR:
```
mkdir build
cmake -S path/to/source -B build -DCMAKE_INSTALL_PREFIX=/where/to/install -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build
```

HemeLB is configured with CMake. For a basic guide see:
https://cmake.org/runningcmake/

You can choose to use the top level "superbuild" which will
automatically download and install any missing dependencies (of those
marked with an asterisk above). If you have any of these installed on
your machine (e.g. via a module system) then they can be used here
also.

The system has three options for each package:
- `System`: never build and use the system one instead
- `Build`: always build and use that
- `Auto`: look for the system package and use is available, otherwise
build.
The default is `Auto`. This is controlled through variables named
`DEPS_${PACKAGE_NAME}`, e.g. `DEPS_TINYXML`.

Alternatively, you can use the code-only build, where you must
ensure that all the dependencies are available. To use this, simply
pass the path to the `Code` subdirectory to CMake as the source.

For development, we suggest that you first install any missing
dependencies using the CMakeLists.txt in the dependencies directory
(supply that path to cmake) and then do a code-only build.

In all cases, the options you can pass to CMake are given in
<CMakeOptions.md>.

We are adding machine-specific instructions to the folder
<machine-specific-build-notes/>.
