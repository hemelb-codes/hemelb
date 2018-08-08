# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_colloidsBCs
  BoundaryConditions.cc
  )
target_link_libraries(hemelb_colloidsBCs PRIVATE
  hemelb_io
  )
hemelb_add_target_dependency_parmetis(hemelb_colloidsBCs)
add_library(hemelb_colloidsForces
  BodyForces.cc
  )
target_link_libraries(hemelb_colloidsForces PRIVATE
  hemelb_io
  )
hemelb_add_target_dependency_parmetis(hemelb_colloidsForces)
add_library(hemelb_colloids
  ParticleMpiDatatypes.cc
  ParticleSet.cc
  Particle.cc
  PersistedParticle.cc
  ColloidController.cc
  )
target_link_libraries(hemelb_colloids PRIVATE
  hemelb_io
  hemelb_lb
  hemelb_colloidsForces
  hemelb_colloidsBCs
  )
hemelb_add_target_dependency_parmetis(hemelb_colloids)
