hemelb_dependency(zlib find)

macro(hemelb_add_target_dependency_zlib tgt)
  target_include_directories(${tgt} PRIVATE ${ZLIB_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${PARMETIS_LIBRARIES})
endmacro()
