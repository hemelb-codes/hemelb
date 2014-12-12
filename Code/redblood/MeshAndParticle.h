//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_MESH_AND_PARTICLE_H
#define HEMELB_REDBLOOD_MESH_AND_PARTICLE_H

#include <vector>
#include "units.h"
#include "redblood/Particle.h"
#include "redblood/stencil.h"
#include "redblood/VelocityInterpolation.h"
#include "geometry/LatticeData.h"


namespace hemelb { namespace redblood {
// Implementation details
#include "redblood/MeshAndParticle.impl.h"

//! Displacement of the particles nodes interpolated from lattice velocities
template<class T_KERNEL> void velocitiesOnMesh(
    Particle const &_particle,
    geometry::LatticeData const &_latDat,
    stencil::types _stencil,
    std::vector<LatticePosition> &_displacements
) {
  _displacements.resize(_particle.GetNumberOfNodes());
  details::VelocityNodeLoop<T_KERNEL>(_stencil, _particle, _latDat)
    .loop(details::transform_iterator(_displacements.begin()));
}

//! Computes and Spreads the forces from the particle to the lattice
//! Returns the energy
Dimensionless forcesOnGrid(
    Particle const &_particle,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
);

//! Computes and Spreads the forces from the particle to the lattice
//! Adds in the node-wall interaction. It is easier to add here since we
//! already have a loop over neighboring grid nodes. Assumption is that the
//! interaction distance is smaller or equal to stencil.
//! Returns the energy (excluding node-wall interaction)
template<class LATTICE> Dimensionless forcesOnGridWithWallInteraction(
    Particle const &_particle,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  std::vector<LatticeForceVector> forces(_particle.GetNumberOfNodes(), 0);
  Dimensionless const energy = _particle(forces);

  details::spreadForce2Grid(
      _particle,
      details::SpreadForcesAndWallForces<LATTICE>(
        _particle, forces, _latticeData
      ),
      _stencil
  );
  return energy;
}

}} // hemelb::redblood

#endif
