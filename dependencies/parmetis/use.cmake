hemelb_dependency(parmetis find)

macro(hemelb_add_target_dependency_parmetis tgt)
  target_include_directories(${tgt} PRIVATE ${PARMETIS_INCLUDE_DIRS})
  target_link_libraries(${tgt} PRIVATE ${ZLIB_LIBRARIES})
endmacro()
