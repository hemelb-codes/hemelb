# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

file(GLOB_RECURSE HEMELB_CU_FILES ${CMAKE_SOURCE_DIR}/*.cc)
set_source_files_properties( ${HEMELB_CU_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )

cuda_add_library(hemelb_colloidsBCs
  BoundaryConditions.cc
  )
target_link_libraries(hemelb_colloidsBCs PUBLIC
  hemelb_io
  )
cuda_add_library(hemelb_colloidsForces
  BodyForces.cc
  )
target_link_libraries(hemelb_colloidsForces PUBLIC
  hemelb_io
  )
cuda_add_library(hemelb_colloids
  ParticleMpiDatatypes.cc
  ParticleSet.cc
  Particle.cc
  PersistedParticle.cc
  ColloidController.cc
  )
target_link_libraries(hemelb_colloids PUBLIC
  hemelb_io
  hemelb_lb
  hemelb_colloidsForces
  hemelb_colloidsBCs
  )
