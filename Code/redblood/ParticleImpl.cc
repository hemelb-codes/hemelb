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

// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
namespace hemelb { namespace redblood { namespace {

bool contains(MeshData::t_Facet const &_a,
        MeshData::t_Facet::value_type _v) {
  return _a[0] == _v or _a[1] == _v or _a[2] == _v;
}

inline bool is_negative_or_null(double _in) {
  return _in < 1e-12;
}

// Helper class to avoid explicit indexing over vertices
struct Facet {
  // References nodes of a facet
  LatticePosition const * nodes[3];
  // Indices of nodes in original array
  MeshData::t_Facet const & indices;
  Facet   (MeshData const &_mesh, size_t _index)
        : indices(_mesh.facets[_index]) {
    nodes[0] = &_mesh.vertices[indices[0]];
    nodes[1] = &_mesh.vertices[indices[1]];
    nodes[2] = &_mesh.vertices[indices[2]];
  }

  // returns an edge nodes[i] - nodes[j]
  LatticePosition operator()(size_t i, size_t j) const {
    return (*this)(i) - (*this)(j);
  }
  // returns node i
  LatticePosition const &operator()(size_t i) const {
    return *(nodes[i]);
  }
};

// Facet that also includes forces
struct ForceFacet : public Facet {
  // References forces on a node
  LatticeForceVector * forces_[3];
  ForceFacet   (MeshData const &_mesh, size_t _index,
          std::vector<LatticeForceVector> &_forces)
        : Facet(_mesh, _index) {
    forces_[0] = &_forces[indices[0]];
    forces_[1] = &_forces[indices[1]];
    forces_[2] = &_forces[indices[2]];
  }
  LatticeForceVector & forces(size_t i) { return *(forces_[i]); }
  LatticeForceVector const & forces(size_t i) const { return *(forces_[i]); }
};


typedef std::pair<size_t, size_t> t_IndexPair;
// Computes common nodes for neighboring facets
// This routine will report nonsense if facets are not neighbors
t_IndexPair commonNodes(Facet const &_a, Facet const &_b) {
  // First node differs, hence other two nodes in common
  if(not contains(_b.indices, _a.indices[0]))
    return t_IndexPair(_a.indices[1], _a.indices[2]);
  // First node in common, second node differs
  if(not contains(_b.indices, _a.indices[1]))
    return t_IndexPair(_a.indices[0], _a.indices[2]);
  // First node and second in common
  return t_IndexPair(_a.indices[0], _a.indices[1]);
}
// Returns common edge with specific direction
LatticePosition commonEdge(Facet const &_a, Facet const &_b) {
  t_IndexPair common(commonNodes(_a, _b));
  return _a(common.first, common.second);
}

// Figures out nodes that are not in common
// Returns non-sense if the nodes are not neighbors.
t_IndexPair singleNodes(Facet const &_a, Facet const &_b) {
  t_IndexPair result;
  size_t mappingB[3] = {4, 4, 4};
  for(size_t i(0); i < 3; ++i)
    if(_a.indices[i] == _b.indices[0]) mappingB[0] = i;
    else if(_a.indices[i] == _b.indices[1]) mappingB[1] = i;
    else if(_a.indices[i] == _b.indices[2]) mappingB[2] = i;
    else result.first = i;
  for(size_t i(0); i < 3; ++i)
    if(mappingB[i] == 4) {result.second = i; break;}
  return result;
}

// Computes vector normal to facet
// Order of edges determines direction of normal
// This is taken straight from T. Krueger's code
LatticePosition normal(Facet const &_facet) {
  LatticePosition const edgeA(_facet(0, 1));
  LatticePosition const edgeB(_facet(2, 1));
  return edgeA.Cross(edgeB);
}

LatticePosition unit_normal(Facet const &_facet) {
  return normal(_facet).Normalise();
}

// Computes angle between two facets
Angle angle(LatticePosition const &_a, LatticePosition const &_b) {
  Angle const cosine(_a.Dot(_b));
  if(cosine >= (1e0 - 1e-6) )       return 0e0;
  else if(cosine <= -(1e0 - 1e-6) ) return PI;
  return std::acos(cosine);
}
Angle angle(Facet const &_a, Facet const &_b) {
  return angle(unit_normal(_a), unit_normal(_b));
}
// Angle angle(MeshData const &_mesh, size_t _facet, size_t _neighbor) {
//   return angle(Facet(_mesh, _facet), Facet(_mesh, _neighbor));
// }

// Angle with orientation
// Computes angle between two vectors, including orientation.
// The orientation should be a vector with an out-of-plane component (eg
// parallel to the cross product of the two normals)
Angle orientedAngle(LatticePosition const &_a, LatticePosition const &_b,
        LatticePosition const &_orient) {
  Angle const result(angle(_a, _b));
  return _a.Cross(_b).Dot(_orient) <= 0e0 ? result: -result;
}
Angle orientedAngle(Facet const &_a, Facet const &_b,
        LatticePosition const &_orient) {
  return orientedAngle(unit_normal(_a), unit_normal(_b), _orient);
}
Angle orientedAngle(Facet const &_a, Facet const &_b) {
  return orientedAngle(_a, _b, commonEdge(_a, _b));
}
Angle orientedAngle(MeshData const &_mesh, size_t _facet, size_t _neighbor) {
  return orientedAngle(Facet(_mesh, _facet), Facet(_mesh, _neighbor));
}



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

  LatticePosition const normalA = normal(facetA);
  LatticePosition const normalB = normal(facetB);

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

  const PhysicalForce strength = _intensity * (theta - theta0);
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

}}}

