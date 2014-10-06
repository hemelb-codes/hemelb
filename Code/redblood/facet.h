#ifndef HEMELB_REDBLOOD_FACET_H
#define HEMELB_REDBLOOD_FACET_H

#include <cmath>
#include <iomanip>
#include "constants.h"

// Helper functions in anonymous namespace.
// These are located in separate file so we can easily unit-test them.
// They are in an anonymous namespace since they likely not needed as a public
// API.
namespace hemelb { namespace redblood { namespace {

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
  LatticePosition operator()(size_t i, size_t j) const;
  // returns node i
  LatticePosition const &operator()(size_t i) const;
  // Returns edges
  LatticePosition edge(size_t _i) const;
  // Returns edge length
  PhysicalDistance length(size_t _i) const;
  // Returns angle cosine
  Dimensionless cosine() const;
  // Returns angle sine
  Dimensionless sine() const;
  // Computes vector normal to facet
  // Order of edges determines direction of normal
  // This is taken straight from T. Krueger's code
  LatticePosition normal() const;
  // Unit vector normal to facet
  LatticePosition unitNormal() const;
  // Area of the facet
  PhysicalSurface area() const;
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
  LatticeForceVector & forces(size_t i) const { return *(forces_[i]); }
};

LatticePosition Facet::operator()(size_t i, size_t j) const {
  return (*this)(i) - (*this)(j);
}
LatticePosition const &Facet::operator()(size_t i) const {
  return *(nodes[i]);
}
LatticePosition Facet::edge(size_t _i) const {
   switch(_i) {
     case 0: return operator()(2, 1);
     case 1: return operator()(0, 1);
     case 2: return operator()(2, 0);
     default: return LatticePosition(0, 0, 0);
   };
}
PhysicalDistance Facet::length(size_t _i) const {
  return edge(_i).GetMagnitude();
}
Dimensionless Facet::cosine() const {
  return edge(0).Dot(edge(1)) / (length(0) * length(1));
}
Dimensionless Facet::sine() const {
  // 0 < angle < 180 in a triangle, so sine always positive in case of
  // interest. I think.
  return edge(0).Cross(edge(1)).GetMagnitude() / (length(0) * length(1));
}
LatticePosition Facet::normal() const {
  LatticePosition const edgeA(operator()(0, 1));
  LatticePosition const edgeB(operator()(2, 1));
  return edgeA.Cross(edgeB);
}
LatticePosition Facet::unitNormal() const { return normal().Normalise(); }
PhysicalSurface Facet::area() const { return normal().GetMagnitude() * 0.5; }


// Just a helper function to check whether _v is in _a
bool contains(MeshData::t_Facet const &_a,
        MeshData::t_Facet::value_type _v) {
  return _a[0] == _v or _a[1] == _v or _a[2] == _v;
}

// Computes common nodes for neighboring facets
// This routine will report nonsense if facets are not neighbors
typedef std::pair<size_t, size_t> t_IndexPair;
t_IndexPair commonNodes(Facet const &_a, Facet const &_b) {
  // First node differs, hence other two nodes in common
  if(not contains(_b.indices, _a.indices[0]))
    return t_IndexPair(1, 2);
  // First node in common, second node differs
  if(not contains(_b.indices, _a.indices[1]))
    return t_IndexPair(0, 2);
  // First node and second in common
  return t_IndexPair(0, 1);
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


// Computes angle between two facets
Angle angle(LatticePosition const &_a, LatticePosition const &_b) {
  Angle const cosine(_a.Dot(_b));
  if(cosine >= (1e0 - 1e-6) )       return 0e0;
  else if(cosine <= -(1e0 - 1e-6) ) return PI;
  return std::acos(cosine);
}
Angle angle(Facet const &_a, Facet const &_b) {
  return angle(_a.unitNormal(), _b.unitNormal());
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
  return orientedAngle(_a.unitNormal(), _b.unitNormal(), _orient);
}
Angle orientedAngle(Facet const &_a, Facet const &_b) {
  return orientedAngle(_a, _b, commonEdge(_a, _b));
}

// Returns Dxx, Dyy, Dxy packed in vector
LatticePosition displacements(Facet const &_deformed, Facet const &_ref) {
  PhysicalDistance const
    dlength0 = _deformed.length(0),
    rlength0 = _ref.length(0),
    dlength1 = _deformed.length(1),
    rlength1 = _ref.length(1);
  Dimensionless const dsine = _deformed.sine(), rsine = _ref.sine();
  return LatticePosition(
      // Dxx
      dlength0 / rlength0,
      // Dyy
      (dlength1 * dsine) / (rlength1 * rsine),
      // Dxy
      (
         dlength1 / rlength1 * _deformed.cosine()
         - dlength0 / rlength0 * _ref.cosine()
      ) / rsine
  );
}

// Returns Gxx, Gyy, Gxy packed in vector
LatticePosition squaredDisplacements(LatticePosition const &_disp) {
  return LatticePosition(
      // Gxx
      _disp[0] * _disp[0],
      // Gyy
      _disp[2] * _disp[2] + _disp[1] * _disp[1],
      // Gxy
      _disp[0] * _disp[2]
  );
}
LatticePosition squaredDisplacements(Facet const &_deformed, Facet const &_ref)
  { return squaredDisplacements(displacements(_deformed, _ref)); }

// Strain invariants I1 and I2
std::pair<Dimensionless, Dimensionless>
  strainInvariants(LatticePosition const &_squaredDisp) {
    return std::pair<Dimensionless, Dimensionless>(
        _squaredDisp[0] + _squaredDisp[1] - 2.0,
        _squaredDisp[0] * _squaredDisp[1]
        - _squaredDisp[2] * _squaredDisp[2] - 1e0
    );
  }

std::pair<Dimensionless, Dimensionless>
  strainInvariants(Facet const &_deformed, Facet const &_ref) {
    return strainInvariants(squaredDisplacements(_deformed, _ref));
}

}}}

#endif
