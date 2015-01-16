//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include "redblood/Cell.h"
// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
#include "redblood/Cell.impl.cc"

namespace hemelb { namespace redblood {

PhysicalEnergy Cell::operator()() const {
  return facetBending_() // facet bending unaffected by template scale
    + volumeEnergy(vertices_, *template_.GetData(), moduli.volume, scale_)
    + surfaceEnergy(vertices_, *template_.GetData(), moduli.surface, scale_)
    + strainEnergy(
        vertices_, *template_.GetData(),
        moduli.dilation, moduli.strain,
        scale_
      );
}
PhysicalEnergy Cell::operator()(
    std::vector<LatticeForceVector> &_forces) const {
  assert(_forces.size() == vertices_.size());
  return facetBending_(_forces)
    + volumeEnergy(
        vertices_, *template_.GetData(), moduli.volume, _forces, scale_)
    + surfaceEnergy(
        vertices_, *template_.GetData(), moduli.surface, _forces, scale_)
    + strainEnergy(
        vertices_, *template_.GetData(),
        moduli.dilation, moduli.strain,
        _forces,
        scale_
      );
}

PhysicalEnergy Cell::facetBending_() const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::t_FacetNeighbors::const_iterator t_FacetIterator;
  t_FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
  t_FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += facetBending(
            vertices_, *template_.GetData(),
            current, (*i_facet)[i],
            moduli.bending
        );
  }
  return result;
}

PhysicalEnergy Cell::facetBending_(
    std::vector<LatticeForceVector> &_forces) const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::t_FacetNeighbors::const_iterator t_FacetIterator;
  t_FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
  t_FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += facetBending(
            vertices_, *template_.GetData(),
            current, (*i_facet)[i],
            moduli.bending, _forces
        );
  }
  return result;
}

void CellBase::operator*=(Dimensionless const &_scale) {
  auto const barycenter = GetBarycenter();
  for(auto &vertex: vertices_)
    vertex = (vertex - barycenter) * _scale + barycenter;
}
void CellBase::operator*=(util::Matrix3D const &_scale) {
  auto const barycenter = GetBarycenter();
  for(auto &vertex: vertices_) {
    _scale.timesVector(vertex - barycenter, vertex);
    vertex += barycenter;
  }
}
void CellBase::operator+=(LatticePosition const &_offset) {
  for(auto &vertex: vertices_)
    vertex += _offset;
}
void CellBase::operator+=(std::vector<LatticePosition> const &_displacements) {
  assert(_displacements.size() == vertices_.size());
  auto i_disp = _displacements.begin();
  for(auto &vertex: vertices_)
      vertex += *(i_disp++);
}

}} // hemelb::rbc
