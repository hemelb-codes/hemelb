//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <vector>
#include "redblood/MeshAndParticle.h"
// Separated some of the implementation so it can be unit-tested and still
// remain unexported within an anonymous namespace
#include "redblood/MeshAndParticle.impl.cc"

namespace hemelb { namespace redblood {

//! Spreads the forces from the particle to the lattice
Dimensionless forces_on_grid(
    Particle const &_particle,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  std::vector<LatticeForceVector> forces(_particle.GetNumberOfNodes(), 0);
  Dimensionless const energy = _particle(forces);
  spread_forces_to_grid(_particle, forces, _latticeData, _stencil);
  return energy;
}

}} // namespace hemelb::redblood
