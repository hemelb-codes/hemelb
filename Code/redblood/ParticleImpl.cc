//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <cmath>
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

//! Helper class to avoid explicit indexing over vertices
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
};

// Computes common nodes for neighboring facets
// This routine will report nonsense if facets are not neighbors
std::pair<size_t, size_t> common_nodes(Facet const &_a, Facet const &_b) {
    // First node differs, hence other two nodes in common
    if(not contains(_b.indices, _a.indices[0]))
        return std::pair<size_t, size_t>(_a.indices[1], _a.indices[2]);
    // First node in common, second node differs
    if(not contains(_b.indices, _a.indices[1]))
        return std::pair<size_t, size_t>(_a.indices[0], _a.indices[2]);
    // First node and second in common
    return std::pair<size_t, size_t>(_a.indices[0], _a.indices[1]);
}
// Returns common edge with specific direction
LatticePosition common_edge(Facet const &_a, Facet const &_b) {
    std::pair<size_t, size_t> common(common_nodes(_a, _b));
    return *_a.nodes[common.first] - *_a.nodes[common.second];
}

// Computes vector normal to facet
// Order of edges determines direction of normal
// This is taken straight from T. Krueger's code
LatticePosition normal(Facet const &_facet) {
    LatticePosition const edgeA(*_facet.nodes[0] - *_facet.nodes[1]);
    LatticePosition const edgeB(*_facet.nodes[2] - *_facet.nodes[1]);
    return edgeA.Cross(edgeB);
}

LatticePosition unit_normal(Facet const &_facet) {
    return normal(_facet).Normalise();
}

// Computes angle between two facets
Angle angle(LatticePosition const &_a, LatticePosition const &_b) {
    Angle const cosine(_a.Dot(_b));
    if(cosine >= (1e0 - 1e-6) )
        return 0e0;
    else if(cosine <= -(1e0 - 1e-6) )
        return PI;
    return std::acos(cosine);
}
Angle angle(Facet const &_a, Facet const &_b) {
    return angle(unit_normal(_a), unit_normal(_b));
}
Angle angle(MeshData const &_mesh, size_t _facet, size_t _neighbor) {
    return angle(Facet(_mesh, _facet), Facet(_mesh, _neighbor));
}

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
    return orientedAngle(_a, _b, common_edge(_a, _b));
}
Angle orientedAngle(MeshData const &_mesh, size_t _facet, size_t _neighbor) {
    return orientedAngle(Facet(_mesh, _facet), Facet(_mesh, _neighbor));
}



// Facet bending energy between two facets
PhysicalEnergy facetBending(MeshData const& _mesh, MeshData const& _orig,
            size_t _facet_index, size_t _neighbor_index) {
    Angle const theta0 = orientedAngle(_orig, _facet_index, _neighbor_index);
    Angle const theta = orientedAngle(_mesh, _facet_index, _neighbor_index);
    return (theta - theta0) * (theta - theta0);
}


}}}

