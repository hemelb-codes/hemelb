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
  return facetBending() // facet bending unaffected by template scale
    + volumeEnergy(vertices, *templateMesh.GetData(), moduli.volume, scale)
    + surfaceEnergy(vertices, *templateMesh.GetData(), moduli.surface, scale)
    + strainEnergy(
        vertices, *templateMesh.GetData(),
        moduli.dilation, moduli.strain,
        scale
      );
}
PhysicalEnergy Cell::operator()(
    std::vector<LatticeForceVector> &forces) const {
  assert(forces.size() == vertices.size());
  return facetBending(forces)
    + volumeEnergy(
        vertices, *templateMesh.GetData(), moduli.volume, forces, scale)
    + surfaceEnergy(
        vertices, *templateMesh.GetData(), moduli.surface, forces, scale)
    + strainEnergy(
        vertices, *templateMesh.GetData(),
        moduli.dilation, moduli.strain,
        forces,
        scale
      );
}

PhysicalEnergy Cell::facetBending() const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::FacetNeighbors::const_iterator FacetIterator;
  FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
  FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += hemelb::redblood::facetBending(
            vertices, *templateMesh.GetData(),
            current, (*i_facet)[i],
            moduli.bending
        );
  }
  return result;
}

PhysicalEnergy Cell::facetBending(
    std::vector<LatticeForceVector> &forces) const {
  if(std::abs(moduli.bending) < 1e-8) return 0e0;

  PhysicalEnergy result(0);
  typedef MeshTopology::FacetNeighbors::const_iterator FacetIterator;
  FacetIterator i_facet = GetTopology()->facetNeighbors.begin();
  FacetIterator const i_facetEnd = GetTopology()->facetNeighbors.end();
  for(size_t current(0); i_facet != i_facetEnd; ++i_facet, ++current) {
    for(size_t i(0); i < 3; ++i)
      if((*i_facet)[i] > current)
        result += hemelb::redblood::facetBending(
            vertices, *templateMesh.GetData(),
            current, (*i_facet)[i],
            moduli.bending, forces
        );
  }
  return result;
}

void CellBase::operator*=(Dimensionless const &scaleIn) {
  auto const barycenter = GetBarycenter();
  for(auto &vertex: vertices)
    vertex = (vertex - barycenter) * scaleIn + barycenter;
}
void CellBase::operator*=(util::Matrix3D const &rotation) {
  auto const barycenter = GetBarycenter();
  for(auto &vertex: vertices) {
    rotation.timesVector(vertex - barycenter, vertex);
    vertex += barycenter;
  }
}
void CellBase::operator+=(LatticePosition const &offset) {
  for(auto &vertex: vertices)
    vertex += offset;
}
void CellBase::operator+=(std::vector<LatticePosition> const &displacements) {
  assert(displacements.size() == vertices.size());
  auto i_disp = displacements.begin();
  for(auto &vertex: vertices)
      vertex += *(i_disp++);
}

}} // hemelb::rbc
