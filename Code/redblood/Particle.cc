//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "Particle.h"
// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
#include "ParticleImpl.cc"

namespace hemelb { namespace redblood {

PhysicalEnergy Particle::operator()() const {
  return facetBending_()
    + volumeEnergy(*mesh_, *template_, moduli.volume)
    + surfaceEnergy(*mesh_, *template_, moduli.surface)
    + strainEnergy(*mesh_, *template_, moduli.dilation, moduli.strain);
}
PhysicalEnergy Particle::operator()(
    std::vector<LatticeForceVector> &_forces) const {
  assert(_forces.size() == mesh_->vertices.size());
  return facetBending_(_forces)
    + volumeEnergy(*mesh_, *template_, moduli.volume, _forces)
    + surfaceEnergy(*mesh_, *template_, moduli.surface, _forces)
    + strainEnergy(*mesh_, *template_,
        moduli.dilation, moduli.strain, _forces);
}

PhysicalEnergy Particle::facetBending_() const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::t_FacetNeighbors::const_iterator t_FacetIterator;
  t_FacetIterator i_facet = topology_->facetNeighbors.begin();
  t_FacetIterator const i_facetEnd = topology_->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += facetBending(*mesh_, *template_, current, (*i_facet)[i],
            moduli.bending);
  }
  return result;
}

PhysicalEnergy Particle::facetBending_(
    std::vector<LatticeForceVector> &_forces) const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::t_FacetNeighbors::const_iterator t_FacetIterator;
  t_FacetIterator i_facet = topology_->facetNeighbors.begin();
  t_FacetIterator const i_facetEnd = topology_->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += facetBending(*mesh_, *template_, current, (*i_facet)[i],
            moduli.bending, _forces);
  }
  return result;
}

}} // hemelb::rbc
