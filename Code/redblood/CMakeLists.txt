add_library(hemelb_redblood Mesh.h Mesh.cc Cell.cc Interpolation.cc stencil.cc
  CellCell.cc FlowExtension.cc FaderCell.cc parallelization/ParticleShuffler.cc)
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  target_compile_features( hemelb_redblood
    PUBLIC
      cxx_auto_type cxx_lambdas cxx_range_for cxx_override
      cxx_delegating_constructors cxx_inheriting_constructors
  )
endif()
