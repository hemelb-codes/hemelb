# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
if (HEMELB_BUILD_RBC)
hemelb_dependency(cppunit use)
# hemelb_dependency(vtk use)

function(create_test_executable target includename)
  set(HEMELB_UNITTEST_INCLUDE ${includename})
  configure_file(main.in.cc "${CMAKE_CURRENT_BINARY_DIR}/${target}.cc")
  add_executable(${target} "${CMAKE_CURRENT_BINARY_DIR}/${target}.cc")
  target_link_libraries(${target}
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    TinyXML::TinyXML
    CPPUnit::CPPUnit
    ${Boost_LIBRARIES}
    ${CTEMPLATE_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${MPWide_LIBRARIES}
    #Because on some systems CPPUNIT needs to be linked to libdl
    ${CMAKE_DL_LIBS}
    ${VTK_LIBRARIES}
  )  
  set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
endfunction()

function(add_test_executable target includename)
  create_test_executable(${target} ${includename})
  add_test(
    NAME ${target}
    COMMAND
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
        $<TARGET_FILE:${target}> -o ${HEMELB_TESTDIR}/${target}.xml
  )
endfunction()

function(add_mpi_test_executable target includename)
  create_test_executable(mpi_${target} ${includename})
  if(MPIEXEC_MAX_NUMPROCS GREATER 3 AND "${HEMELB_STEERING_LIB}" STREQUAL "none")
    add_test(NAME mpi_${target}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
        $<TARGET_FILE:mpi_${target}> -o ${HEMELB_TESTDIR}/mpi_${target}.xml
    )
    set_target_properties(mpi_${target} PROPERTIES LABEL mpi)
    set_tests_properties (mpi_${target} PROPERTIES TIMEOUT 60)
  endif()
  if(NOT MPIEXEC_MAX_NUMPROCS GREATER 3)
    message(WARNING "MPIEXEC_MAX_NUMPROCS = ${MPIEXEC_MAX_NUMPROCS}"
      " yet redblood cells mpi tests require > 4: "
      " disabling test mpi_${target}")
  endif()
  if(NOT "${HEMELB_STEERING_LIB}" STREQUAL "none")
    message(WARNING "Cannot run mpi tests unless steering is disabled."
      " Please set HEMELB_STEERING_LIB to none."
      " Disabling test mpi_${target} until then.")
  endif()
endfunction()

if(NOT HEMELB_TESTDIR)
  set(HEMELB_TESTDIR "${PROJECT_BINARY_DIR}/unittests")
  file(MAKE_DIRECTORY "${HEMELB_TESTDIR}")
endif()

add_test_executable(unittests_redblood redblood.h)
# add_test_executable(unittests_redblood_integrations redblood_integrations.h)
add_mpi_test_executable(unittests_redblood redblood_mpi.h)
add_mpi_test_executable(unittests_redblood_integration redbloodMPIIntegrations.h)

add_to_resources(
  resources/red_blood_cell.txt resources/red_blood_cube.txt
  resources/large_cylinder.gmy resources/large_cylinder.xml
  resources/large_cylinder_rbc.xml resources/iccb_capillary_network.xml
  resources/P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_388889_1000_3.gmy
  resources/P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_875_1000_3.gmy
  resources/P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_1_75_1000_3.gmy
  resources/rbc_ico_1280.msh resources/rbc_ico_720.msh
  resources/fedosov1c.xml resources/fedosov1c.gmy
  resources/sad.msh
  resources/cyl_l100_r5.xml
  resources/cyl_l100_r5.gmy
  resources/rbc_ico_2880.msh
  resources/rbc_ico_720_correct.msh
  resources/rbc_ico_720.msh
  resources/rbc_ico_720.vtp
  resources/992Particles_rank3_26_t992.msh
  resources/992Particles_rank3_26_t992.vtp
  resources/empty_for_relative_paths.xml
)
endif()
