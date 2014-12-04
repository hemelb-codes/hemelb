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

namespace hemelb { namespace redblood { namespace {
//! Iterates over vertices of a mesh and the nearby nodes of a grid
//! The functor argument is called with the current vertex index, the
//! global site index triplet, and the associated interpolation weight.
template<class T_FUNCTOR> void spread_forces_to_grid(
    Particle const &_particle,
    T_FUNCTOR _functor,
    stencil::types _stencil
) {
  typedef MeshData::t_Vertices::const_iterator const_iterator;
  // Spread them onto lattice
  const_iterator i_vertex = _particle.GetVertices().begin();
  const_iterator const i_end = _particle.GetVertices().end();
  for(size_t i(0); i_vertex != i_end; ++i_vertex, ++i) {
    InterpolationIterator spreader
      = interpolationIterator(*i_vertex, _stencil);
    for(; spreader; ++spreader)
      _functor(i, *spreader, spreader.weight());
  }
}

class SpreadForces {
  public:
    SpreadForces(
        std::vector<LatticePosition> const &_forces,
        geometry::LatticeData &_latticeData
    ) : latticeData_(_latticeData), forces_(_forces) {}

    void operator()(
        size_t _vertex, LatticeVector const &_site,
        Dimensionless _weight) {
      if(latticeData_.IsValidLatticeSite(_site)) {
        geometry::Site<geometry::LatticeData>
          site(latticeData_.GetSite(_site));
        site.AddToForce(forces_[_vertex] * _weight);
      }
    }

  protected:
    geometry::LatticeData &latticeData_;
    std::vector<LatticeForceVector> const & forces_;
};

void spread_forces_to_grid(
    Particle const &_particle,
    std::vector<LatticeForceVector> const &_forces,
    geometry::LatticeData &_latticeData,
    stencil::types _stencil
) {
  assert(_forces.size() == _particle.GetNumberOfNodes());
  spread_forces_to_grid(
      _particle, SpreadForces(_forces, _latticeData), _stencil);
}
}}} // namespace hemelb::redblood::anonymous
