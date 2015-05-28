add_library(hemelb_redblood Mesh.h Mesh.cc Cell.cc Interpolation.cc stencil.cc
  CellCell.cc FlowExtension.cc FaderCell.cc RBCInserter.cpp
  VertexBag.cc
  buffer/Buffer.cc buffer/Columns.cc)
