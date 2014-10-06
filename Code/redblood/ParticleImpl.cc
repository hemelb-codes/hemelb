//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <cmath>
#include <iomanip>
#include "Mesh.h"
#include "constants.h"
#include "facet.h"

// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
namespace hemelb { namespace redblood { namespace {

// Facet bending energy between two facets
PhysicalEnergy facetBending(MeshData const& _mesh, MeshData const& _orig,
            size_t _facet_index, size_t _neigh_index,
            PhysicalForce _intensity) {
  Angle const theta0 = orientedAngle(_orig, _facet_index, _neigh_index);
  Angle const theta = orientedAngle(_mesh, _facet_index, _neigh_index);
  return _intensity * (theta - theta0) * (theta - theta0);
}

// Facet bending energy and force between neighboring facets
PhysicalEnergy facetBending(MeshData const& _mesh, MeshData const& _orig,
  size_t _facet_index, size_t _neighbor_index, PhysicalForce _intensity,
  std::vector<LatticeForceVector> &_forces) {

  ForceFacet facetA(_mesh, _facet_index, _forces);
  ForceFacet facetB(_mesh, _neighbor_index, _forces);

  t_IndexPair const commons = commonNodes(facetA, facetB);
  t_IndexPair const singles = singleNodes(facetA, facetB);

  LatticePosition const normalA = facetA.normal();
  LatticePosition const normalB = facetB.normal();

  PhysicalDistance const inverseAreaA = 1e0 / normalA.GetMagnitude();
  PhysicalDistance const inverseAreaB = 1e0 / normalB.GetMagnitude();

  // Orthogonalize normal vectors and normalize to inverse area of other facet
  LatticePosition const unitA(normalA.GetNormalised());
  LatticePosition const unitB(normalB.GetNormalised());
  LatticePosition const vecA(
      (unitB - unitA * unitA.Dot(unitB)).GetNormalised() * inverseAreaA
  );
  LatticePosition const vecB(
      (unitA - unitB * unitA.Dot(unitB)).GetNormalised() * inverseAreaB
  );

  // NOTE: the two lines below could make use of stuff computed previously
  Angle const theta0 = orientedAngle(_orig, _facet_index, _neighbor_index);
  Angle const theta = orientedAngle(_mesh, _facet_index, _neighbor_index);

  const PhysicalForce strength = -2.0 * _intensity * (theta - theta0);
  // forces on nodes that are in common
  facetA.forces(commons.first) += (
      facetA(singles.first, commons.second).Cross(vecA)
      + (facetA(commons.second) - facetB(singles.second)).Cross(vecB)
  ) * strength;
  facetA.forces(commons.second) += (
      (facetB(singles.second) - facetA(commons.first)).Cross(vecB)
      + (facetA(commons.first, singles.first)).Cross(vecA)
  ) * strength;
  // forces on nodes that are *not* in common
  facetA.forces(singles.first) += (
      facetA(commons.second, commons.first).Cross(vecA)
  ) * strength;
  facetB.forces(singles.second) += (
      facetA(commons.first, commons.second).Cross(vecB)
  ) * strength;

  return _intensity * (theta - theta0) * (theta - theta0);
}

PhysicalEnergy volumeEnergy(MeshData const &_mesh, MeshData const &_orig,
    PhysicalForce _intensity) {
  PhysicalVolume const vol0 = volume(_orig);
  PhysicalVolume const deltaV = volume(_mesh) - vol0;
  return _intensity * 0.5 * deltaV * deltaV / vol0;
}

PhysicalEnergy volumeEnergy(MeshData const &_mesh, MeshData const &_orig,
    PhysicalForce _intensity, std::vector<LatticeForceVector> &_forces) {
  MeshData::t_Facets::const_iterator i_facet = _mesh.facets.begin();
  MeshData::t_Facets::const_iterator const i_facet_end = _mesh.facets.end();
  assert(_mesh.facets.size() == _orig.facets.size());

  PhysicalVolume const vol0 = volume(_orig);
  PhysicalVolume const deltaV = volume(_mesh) - vol0;
  double const strength(_intensity / 6.0 * deltaV / vol0);
  for(; i_facet != i_facet_end; ++i_facet) {
    // Come aliases to make it easier to refer to vertices
    size_t const i0((*i_facet)[0]), i1((*i_facet)[1]), i2((*i_facet)[2]);
    LatticePosition const & a(_mesh.vertices[i0]);
    LatticePosition const & b(_mesh.vertices[i1]);
    LatticePosition const & c(_mesh.vertices[i2]);

    _forces[i0] += b.Cross(c) * strength;
    _forces[i1] += c.Cross(a) * strength;
    _forces[i2] += a.Cross(b) * strength;
  }
  return 0.5 * _intensity * deltaV * deltaV / vol0;
}

PhysicalEnergy surfaceEnergy(MeshData const &_mesh, MeshData const &_orig,
    PhysicalForce _intensity) {
  PhysicalSurface const surf0 = surface(_orig);
  PhysicalSurface const deltaS = surface(_mesh) - surf0;
  return _intensity * 0.5 * deltaS * deltaS / surf0;
}

PhysicalEnergy surfaceEnergy(MeshData const &_mesh, MeshData const &_orig,
    PhysicalForce _intensity, std::vector<LatticeForceVector> &_forces) {

  PhysicalSurface const surf0 = surface(_orig);
  PhysicalSurface const deltaS = surface(_mesh) - surf0;
  double const strength = _intensity * 0.5 * deltaS / surf0;
  for(size_t facetIndex(0); facetIndex < _mesh.facets.size(); ++facetIndex) {
    ForceFacet facet(_mesh, facetIndex, _forces);
    LatticePosition const n0 = facet.unitNormal();

    facet.forces(0) += n0.Cross(facet(2, 1)) * strength;
    facet.forces(1) += n0.Cross(facet(0, 2)) * strength;
    facet.forces(2) += n0.Cross(facet(1, 0)) * strength;
  }
  return _intensity * 0.5 * deltaS * deltaS / surf0;
}

PhysicalEnergy strainEnergyDensity(
    std::pair<Dimensionless, Dimensionless> const &_strainParams,
    PhysicalForce _shearModulus, PhysicalForce _dilationModulus) {
  Dimensionless const I1 = _strainParams.first, I2 = _strainParams.second;
  return _shearModulus / 12. * (I1 * I1 + 2. * I1 - 2. * I2)
        + _dilationModulus / 12. * I2 * I2;
}
PhysicalEnergy strainEnergy(
    Facet const &_deformed, Facet const &_undeformed,
    PhysicalForce _shearModulus, PhysicalForce _dilationModulus
) {
  return strainEnergyDensity(
      strainInvariants(_deformed, _undeformed),
      _shearModulus, _dilationModulus
  ) * _undeformed.area();
}

PhysicalEnergy strainEnergy(
    ForceFacet const &_deformed, Facet const &_undeformed,
    PhysicalForce _shearModulus, PhysicalForce _dilationModulus
) {
  // Shape function parameters
  Dimensionless const
    b0 = _undeformed.length(0) * 0.5,
    b1 = (
      _undeformed.length(1) * _undeformed.cosine() - _undeformed.length(0)
    ) * 0.5,
    a1 = -0.5 * _undeformed.length(1) * _undeformed.sine();

  LatticePosition const disps(displacements(_deformed, _undeformed));
  LatticePosition const squaredDisps(squaredDisplacements(disps));
  std::pair<Dimensionless, Dimensionless> const
    strainInvs(strainInvariants(squaredDisps));
  Dimensionless const I1 = strainInvs.first, I2 = strainInvs.second;
  Dimensionless const w = strainEnergyDensity(
      strainInvs, _shearModulus, _dilationModulus);

  // Skalak Parameters
  PhysicalForce const
    dw_dI1 = _shearModulus / 6 * (I1 + 1),
    dw_dI2 = -_shearModulus / 6. + _dilationModulus / 6. * I2;

  size_t const xx = 0, yy = 1, xy = 2;

  // Derivatives of strain invariants
  Dimensionless const
    dI1_dGxx = 1.,
    dI1_dGyy = 1.,
    dI2_dGxx = squaredDisps[yy],
    dI2_dGyy = squaredDisps[xx],
    dI2_dGxy = -2. * squaredDisps[xy];

  // Derivatives of squared deformation tensor
  Dimensionless const
    dGxx_du1x = 2. * a1 * disps[xx],
    dGxy_du0x = b0 * disps[xx],
    dGxy_du1x = a1 * disps[xy] + b1 * disps[xx],
    dGxy_du1y = a1 * disps[yy],
    dGyy_du0x = 2. * b0 * disps[xy],
    dGyy_du0y = 2. * b0 * disps[yy],
    dGyy_du1x = 2. * b1 * disps[xy],
    dGyy_du1y = 2. * b1 * disps[yy];

  PhysicalForce const
    force0x = (
        dw_dI1 * dI1_dGyy * dGyy_du0x
      + dw_dI2 * (dI2_dGyy * dGyy_du0x + dI2_dGxy * dGxy_du0x)
    ),
    force0y = (
        dw_dI1 * dI1_dGyy * dGyy_du0y
      + dw_dI2 * dI2_dGyy * dGyy_du0y
    ),
    force1x = (
        dw_dI1 * (dI1_dGxx * dGxx_du1x + dI1_dGyy * dGyy_du1x)
      + dw_dI2 * (
          dI2_dGxx * dGxx_du1x + dI2_dGyy * dGyy_du1x + dI2_dGxy * dGxy_du1x
      )
    ),
    force1y = (
        dw_dI1 * dI1_dGyy * dGyy_du1y
      + dw_dI2 * (dI2_dGyy * dGyy_du1y + dI2_dGxy * dGxy_du1y)
    );

    /// Coordinate system
    LatticePosition const ex = _deformed.edge(0).GetNormalised(),
      ez = _deformed.edge(0).Cross(_deformed.edge(1)).GetNormalised(),
      ey = ez.Cross(ex);

    LatticeForceVector const
      force0 = ex * force0x + ey * force0y,
      force1 = ex * force1x + ey * force1y;
    _deformed.forces(0) -= force0;
    _deformed.forces(1) -= force1;
    _deformed.forces(2) += force0 + force1;

    return w * _undeformed.area();
}

PhysicalEnergy strainEnergy(MeshData const &_mesh, MeshData const &_origin,
      PhysicalForce _shearModulus, PhysicalForce _dilationModulus
    ) {
  PhysicalEnergy result(0);
  for(size_t i(0); i < _mesh.facets.size(); ++i)
    result += strainEnergy(Facet(_mesh, i), Facet(_origin, i),
        _shearModulus, _dilationModulus);
  return result;
}
PhysicalEnergy strainEnergy(MeshData const &_mesh, MeshData const &_origin,
      PhysicalForce _shearModulus, PhysicalForce _dilationModulus,
      std::vector<LatticeForceVector> &_forces
    ) {
  PhysicalEnergy result(0);
  for(size_t i(0); i < _mesh.facets.size(); ++i)
    result += strainEnergy(ForceFacet(_mesh, i, _forces), Facet(_origin, i),
        _shearModulus, _dilationModulus);
  return result;
}

}}}

