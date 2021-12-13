# HemeLB's CMake options

## Superbuild vs code-only
You can choose to use the top level "superbuild" which will
automatically download and install any missing dependencies (of those
marked with an asterisk above). If you have any of these installed on
your machine (e.g. via a module system) then they can be used here
also.

Alternatively, you can use the code-only build, where you must
ensure that all the dependencies are availble.

## Installation locations

As well as the usual `CMAKE_INSTALL_PREFIX` variable, you can also set
`HEMELB_DEPENDENCIES_INSTALL_PREFIX`, which will tell the superbuild
where to install any dependencies that it compiles and give an extra
search location to the code-only build.

## Dependency variables

You can tell either build where to look for dependencies if CMake
cannot find them automatically.

- Boost: `BOOST_ROOT` for the prefix or `BOOST_INCLUDEDIR` / `BOOST_LIBRARYDIR`

- Catch2:

- CTemplate: `CTEMPLATE_INCLUDE_DIR` / `CTEMPLATE_LIBRARIES`

- HDF5

- METIS (required by ParMETIS and sometimes installed independently):
  `METIS_ROOT`/ `METIS_DIR` for prefix or `METIS_INCLUDE_DIR`/
  `METIS_LIBRARY`

- ParMETIS: as above but substitute `PARMETIS` for `METIS`

- TinyXML: `TINYXML_INCLUDE_DIR` / `TINYXML_LIBRARIES`

- VTK: `VTK_DIR`

You can see these by using the CMake interactive CLI `ccmake`.

## Resolved red blood cells

To include modelling of fully resolved red blood cells via the
immersed boundary method, set `HEMELB_BUILD_RBC=ON`.


## Lattice Boltzmann options
The HemeLB-specific options and variables are all given and briefly
documented in Code/cmake/options.cmake - please see that file for details.

HemeLB supports multiple lattice Boltzmann velocity sets, collisions,
boundary condition etc, but which is active is chosen at compile time
by setting these options.

- `HEMELB_EXECUTABLE` sets the name of the produced application. By
  default this is `hemelb` but you might wish to add a suffix if you are
  experimenting with multiple LB models etc.

- The lattice or velocity set is chosen with `HEMELB_LATTICE` and can
  be one of D3Q15 (default), D3Q19, D3Q27, D3Q15i

- The collision kernel is chosen with `HEMELB_KERNEL` from LBGK
  (default), EntropicAnsumali, EntropicChik, MRT, TRT, NNCY, NNCYMOUSE,
  NNC, NNTPL

- The no-slip solid wall boundary is selected with
  `HEMELB_WALL_BOUNDARY` from BFL, GZS, SIMPLEBOUNCEBACK (default),
  JUNKYANG

- The inlet and outlet boundary conditions are chosen by
  `HEMELB_INLET_BOUNDARY` and `HEMELB_OUTLET_BOUNDARY` respectively from
  NASHZEROTHORDERPRESSUREIOLET (default), LADDIOLET. It is *very
  important* that you also select `HEMELB_WALL_INLET_BOUNDARY` and
  `HEMELB_WALL_INLET_BOUNDARY` to match the combination of your
  selected wall and in/outlet boundaries. (Options are:
  NASHZEROTHORDERPRESSURESBB, NASHZEROTHORDERPRESSUREBFL, LADDIOLETSBB,
  LADDIOLETBFL)

- `HEMELB_BUILD_MULTISCALE`: enable HemeLB's multiscale coupling mode.
   Requires MPWIde.

## Performance options

- `HEMELB_SUBPROJECT_MAKE_JOBS`: enable parallel builds for the
  superbuild. Set to approximately twice the number of cores available
  for your use. (On some HPC login nodes you may hit limits on number
  of allowed processes/filehandles/etc quite quickly so if your build
  dies mysteriously, especially with internal compiler errors, set
  this to a small number, e.g. 4)

- `HEMELB_USE_SSE3`: this is on by default and enables use of SSE3
  intrinsics. This may not work on your architecture (e.g. ARM)


## Developer

- HemeLB has an option to make small-scale parallel debugging with GDB
  and LLDB easier (typically used on a developers workstation). Turn
  this on with `HEMELB_BUILD_DEBUGGER` and supply the `-d 1` command
  line option.

- `HEMELB_VALIDATE_GEOMETRY`: the code can validate that a geometry file
  is self-consistent on loading.
