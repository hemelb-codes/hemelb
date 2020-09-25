hemelb_dependency(MPWide find)

macro(hemelb_add_target_dependency_mpwide tgt)
  target_include_directories(${tgt} PRIVATE ${MPWide_INCLUDE_DIR})
  target_compile_definitions(${tgt} PRIVATE -DHEMELB_BUILD_MULTISCALE)
  target_link_libraries(${tgt} PRIVATE ${MPWide_LIBRARIES})
endmacro()
