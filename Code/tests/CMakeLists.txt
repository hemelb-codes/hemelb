# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

hemelb_dependency(catch2 use)
add_subdirectory(helpers)

function(create_test_executable target)
  add_executable(${target})
  target_sources(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

  target_compile_definitions(${target} PRIVATE HEMELB_DOING_UNITTESTS)
  # ReporterTests directly use ctemplate
  # FlowExtensionTests directly uses tinyxml
  target_link_libraries(${target} PRIVATE
    Catch2::Catch2
    CTemplate::CTemplate
    TinyXML::TinyXML
    ${heme_libraries}
    ${MPI_LIBRARIES}
    test-helpers
    )
  if (HEMELB_BUILD_RBC)
    # redblood tests need VTK
    target_link_libraries(${target} PRIVATE VTK::CommonCore)
  endif()

  install(TARGETS ${target} RUNTIME DESTINATION bin)
endfunction()

create_test_executable(hemelb-tests)

add_to_resources(resources/four_cube.gmy resources/four_cube.xml resources/four_cube_multiscale.xml
  resources/config.xml resources/config0_2_0.xml
  resources/config_file_inlet.xml resources/iolet.txt
  resources/config-velocity-iolet.xml resources/config_new_velocity_inlets.xml
  resources/velocity_inlet.txt.weights.txt
  resources/xmltest.xml resources/config_file_velocity_inlet.xml resources/velocity_inlet.txt
  )
if (HEMELB_BUILD_MULTISCALE)
  add_to_resources(multiscale/mpwide/MPWSettings.cfg)
endif()

# TODO: when cmake min version increases to 3.13, can drop ${CMAKE_CURRENT_SOURCE_DIR}
target_sources(hemelb-tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/SimulationMasterTests.cc)
add_subdirectory(configuration)
add_subdirectory(extraction)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(lb)
add_subdirectory(multiscale)
add_subdirectory(net)
add_subdirectory(reporting)
add_subdirectory(util)
add_subdirectory(vis)

if (HEMELB_BUILD_RBC)
  create_test_executable(mpi_redblood_tests)
  add_subdirectory(redblood)
endif()
