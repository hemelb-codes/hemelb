add_library(hemelb_redblood Mesh.h Mesh.cc CellBase.cc Cell.cc Interpolation.cc
  CellCell.cc FlowExtension.cc FaderCell.cc io.cc RBCInserter.cc
  VertexBag.cc Borders.cc
  WallCellPairIterator.cc
  buffer/Buffer.cc buffer/Columns.cc
  parallel/ParticleShuffler.cc)
