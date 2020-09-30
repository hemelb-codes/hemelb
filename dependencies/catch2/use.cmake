find_package(Catch2 REQUIRED)

macro(hemelb_add_target_dependency_catch2 tgt)
  target_link_libraries(${tgt} PRIVATE Catch2::Catch2)
endmacro()
