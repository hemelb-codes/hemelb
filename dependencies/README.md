# HemeLB dependency builder

This CMake project automates installing most of the main application's dependencies.

For each dependency, there is a folder with the same name (case sensitive as the CMake package).

In the folder there are two required files: `find.cmake` and `build.cmake`, which do what you might expect.

The for a package called `Foo`, the find file should call `find_package(Foo ${DEP_FIND_MODE_Foo})` where `DEP_FIND_MODE_Foo` is filled in to zero or more of `QUIET` and `REQUIRED` by the calling code. You can add version and component specifiers in the usual way, see the [CMake docs for find_package](https://cmake.org/cmake/help/latest/command/find_package.html).

For the same package, the build file should call `ExternalProject_add` with a target name of `dep_Foo` and whatever options needed. The source tarball should searched for in `dependendencies/distributions/` first before downloading to allow this to work on machines that don't allow outgoing connections.
