hemelb_dependency(ctemplate find)

macro(hemelb_add_target_dependency_ctemplate tgt)
  target_include_directories(${tgt} PRIVATE ${CTEMPLATE_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${CTEMPLATE_LIBRARIES})
endmacro()
