find_package(CPPUnit REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIR})

add_executable(unittests_hemelb main.cc)
# A function to specify different tests using the same main.cc
function(add_include_to_test target filename)
  get_target_property(alldefs ${target} COMPILE_DEFINITIONS)
  set(property HEMELB_UNITTEST_INCLUDE=\"unittests/${filename}\")
  if(NOT "${alldefs}" MATCHES "${property}")
    if(${alldefs})
      set(alldefs "${alldefs} ${property}")
    else()
      set(alldefs "${property}")
    endif()
    set_target_properties(${target}
       PROPERTIES
       COMPILE_DEFINITIONS "${alldefs}"
    )
  endif()
endfunction()

macro(add_test_executable target includename)
  target_link_libraries(${target}
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${PARMETIS_LIBRARIES}
    ${TINYXML_LIBRARIES}
    ${CPPUNIT_LIBRARY}
    ${Boost_LIBRARIES}
    ${CTEMPLATE_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${MPWide_LIBRARIES}
    #Because on some systems CPPUNIT needs to be linked to libdl
    ${CMAKE_DL_LIBS}
  )
  add_include_to_test(${target} ${includename})
  set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
  INSTALL(TARGETS ${target} RUNTIME DESTINATION bin)
  set(testdir "${PROJECT_BINARY_DIR}/test_results/")
  file(MAKE_DIRECTORY "${testdir}")
  add_test(NAME ${target} COMMAND ${target} -o ${testdir}/${target}.xml)
endmacro()

add_test_executable(unittests_hemelb hemelb_unittests.h)


add_to_resources(resources/four_cube.gmy resources/four_cube.xml
  resources/four_cube_multiscale.xml resources/config.xml
  resources/config0_2_0.xml resources/config_file_inlet.xml
  resources/iolet.txt resources/config-velocity-iolet.xml
  resources/config_new_velocity_inlets.xml resources/xmltest.xml
  resources/config_file_velocity_inlet.xml resources/velocity_inlet.txt
  resources/red_blood_cell.txt resources/red_blood_cube.txt
  resources/large_cylinder.gmy resources/large_cylinder.xml
  resources/large_cylinder_rbc.xml resources/iccb_capillary_network.xml
  resources/P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_388889_1000_3.gmy
  resources/P6_190514_x63x0.6_0_RED_BW_corrected_tubed_smoothed_0_777778_1000_3.gmy
)
