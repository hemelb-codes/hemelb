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

namespace hemelb { namespace redblood {

//! Spreads the forces from the particle to the lattice
Dimensionless forcesOnGrid(
    Particle const &_particle,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  std::vector<LatticeForceVector> forces(_particle.GetNumberOfNodes(), 0);
  Dimensionless const energy = _particle(forces);
  details::spreadForce2Grid(
      _particle, details::SpreadForces(forces, _latticeData), _stencil);
  return energy;
}

}} // namespace hemelb::redblood
