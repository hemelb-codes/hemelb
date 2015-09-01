add_library(hemelb_redblood Mesh.h Mesh.cc Cell.cc Interpolation.cc
  CellCell.cc FlowExtension.cc FaderCell.cc io.cc RBCInserter.cc
  VertexBag.cc
  buffer/Buffer.cc buffer/Columns.cc)
