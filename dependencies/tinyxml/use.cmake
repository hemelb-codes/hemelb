hemelb_dependency(tinyxml find)

macro(hemelb_add_target_dependency_tinyxml tgt)
  # STL use is required by HemeLB and the build options 
  target_compile_definitions(${tgt}  PRIVATE -DTIXML_USE_STL)
  target_include_directories(${tgt} PRIVATE ${TINYXML_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${TINYXML_LIBRARIES})
endmacro()

