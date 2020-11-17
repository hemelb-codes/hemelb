find_package(HDF5 1.8 REQUIRED COMPONENTS C HL)

macro(hemelb_add_target_dependency_hdf5 tgt)
  target_link_libraries(${tgt} PRIVATE hdf5::hdf5_hl)
endmacro()
