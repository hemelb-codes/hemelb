hemelb_dependency(cppunit find)

macro(hemelb_add_target_dependency_cppunit tgt)
  target_include_directories(${tgt} PRIVATE ${CPPUNIT_INCLUDE_DIR})
endmacro()
