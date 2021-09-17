# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required(VERSION 3.10)

project(HemeLbGmyTool)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

# Debug print a variable nicely
macro(dvar name)
  message("${name} = ${${name}}")
endmacro()

function(dvar_re regex)
  get_cmake_property(all_var_names VARIABLES)
  list (SORT all_var_names)
  foreach (varname ${all_var_names})
    string(REGEX MATCH "${regex}" matched "${varname}")
    if (matched)
      dvar(${varname})
    endif()
  endforeach()
endfunction()

include(GNUInstallDirs)

find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED COMPONENTS headers)
find_package(VTK REQUIRED COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMath
  vtkFiltersCore
  vtkFiltersGeneral
  vtkFiltersSources
  vtkIOXML
  vtkWrappingPythonCore
)
find_package(CGAL REQUIRED)

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific
  # directory here.
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" -c
      "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _pybind11_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmake_dir}")
endif()

find_package(pybind11 CONFIG REQUIRED)

set(HEMELB_CODE_DIR ${PROJECT_SOURCE_DIR}/../Code)
set(model_dir HlbGmyTool/Model)
set(generation_dir HlbGmyTool/Model/Generation)

# First, the bits of HemeLB main we re-use (micro-hemelb geometry).
set(hemelb_cpps
  util/Vector3D.cc
  geometry/SiteDataBare.cc
  lb/lattices/D3Q27.cc
  io/formats/geometry.cc
  io/writers/xdr/XdrFileWriter.cc
  io/writers/xdr/XdrMemWriter.cc
  io/writers/xdr/XdrWriter.cc
  io/writers/Writer.cc
)
list(TRANSFORM hemelb_cpps PREPEND "${HEMELB_CODE_DIR}/")

# Now the C++ implementation of the extension
set(gen_cpps
  Neighbours.cpp
  BuildCGALPolygon.cpp
  Block.cpp
  BlockWriter.cpp
  BufferPool.cpp
  GeometryGenerator.cpp
  GeometryWriter.cpp
  Index.cpp
  Domain.cpp
  Site.cpp
  InconsistentFluidnessError.cpp
  CylinderGenerator.cpp
  PolyDataGenerator.cpp
  SquareDuctGenerator.cpp
  Debug.cpp
)
list(TRANSFORM gen_cpps PREPEND "${generation_dir}/")

pybind11_add_module(Generation ${hemelb_cpps} ${gen_cpps} "${generation_dir}/Wrap.cpp")
target_include_directories(Generation PRIVATE ${generation_dir} ${HEMELB_CODE_DIR})
target_link_libraries(Generation PRIVATE
  Boost::headers
  ZLIB::ZLIB
  Boost::headers
  CGAL::CGAL
  ${VTK_LIBRARIES}
)

install(TARGETS Generation DESTINATION ${model_dir})
