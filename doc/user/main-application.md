# HemeLB Users' Guide

A simple guide to compile and run HemeLB

## Dependencies

The main HemeLB application requires:

- Unix (-ish) OS: modern Linux or MacOS
- C++17 compiler
- CMake version >= 3.10
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
- VTK* version 9 (for RBC mode)

The dependencies marked with * can be built automatically by the build
system.

## Download

Clone the repository:
```bash
git clone https://github.com/hemelb-codes/hemelb.git
```

## Compile

You can choose to use the top level "superbuild" which will
automatically download and install any missing dependencies (of those
marked with an asterisk above). If you have any of these installed on
your machine (e.g. via a module system) then they can be used here
also. Alternatively, you can use the code-only build, where you must
ensure that all the dependencies are availble.

For superbuild pass the root directory of the repository to CMake and
for a code-only build pass the `Code` subdirectory instead. In either
case, we strongly recommend doing a standard CMake out-of-source
build.

```bash
mkdir build
cmake -B build -S path/to/hemelb/Code <other cmake options>
```

The options you can pass to CMake are given in <CMakeOptions.md>.
