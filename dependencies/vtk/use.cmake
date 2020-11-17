find_package(VTK 8.2 REQUIRED NO_MODULE)

macro(hemelb_add_target_dependency_vtk tgt)
  target_link_libraries(${tgt} PRIVATE VTK::CommonCore VTK::CommonDataModel VTK::IOXML VTK::FiltersCore)
endmacro()
