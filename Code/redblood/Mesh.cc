//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <fstream>
#include <cassert>
#include <numeric>

#include "redblood/Mesh.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "Exception.h"

namespace hemelb { namespace redblood {

boost::shared_ptr<redblood::MeshData> read_mesh(std::string const &_filename) {
  log::Logger::Log<log::Debug, log::Singleton>(
          "Reading red blood cell from %s", _filename.c_str());

  // Open file if it exists
  std::string line;
  std::ifstream file;
  if (!util::file_exists(_filename.c_str()))
    throw Exception() << "Red-blood-cell mesh file '"
      << _filename.c_str() << "' does not exist";
  file.open(_filename.c_str());

  // Drop header
  for(int i(0); i < 4; ++i) std::getline(file, line);

  // Number of vertices
  unsigned int num_vertices;
  file >> num_vertices;


  // Create Mesh data
  boost::shared_ptr<MeshData> result(new MeshData);
  result->vertices.resize(num_vertices);

  // Then read in first and subsequent lines
  MeshData::t_Facet::value_type offset;
  file >> offset >> result->vertices[0].x
    >> result->vertices[0].y
    >> result->vertices[0].z;
  log::Logger::Log<log::Trace, log::Singleton>(
    "Vertex 0 at %d, %d, %d", result->vertices[0].x,
    result->vertices[0].y, result->vertices[0].z
  );
  for(unsigned int i(1), index(0); i < num_vertices; ++i) {
    file >> index >> result->vertices[i].x
      >> result->vertices[i].y
      >> result->vertices[i].z;
    log::Logger::Log<log::Trace, log::Singleton>(
      "Vertex %i at %d, %d, %d", i, result->vertices[i].x,
      result->vertices[i].y, result->vertices[i].z
    );
  }

  // Drop mid-file headers
  for(int i(0); i < 3; ++i) std::getline(file, line);

  // Read facet indices
  unsigned int num_facets;
  MeshData::t_Facet indices;
  file >> num_facets;
  result->facets.resize(num_facets);
  for(unsigned int i(0), dummy(0); i < num_facets; ++i) {
    file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
      >> indices[0] >> indices[1] >> indices[2];
    result->facets[i][0] = indices[0] - offset;
    result->facets[i][1] = indices[1] - offset;
    result->facets[i][2] = indices[2] - offset;
    log::Logger::Log<log::Trace, log::Singleton>(
        "Facet %i with %i, %i, %i", i, indices[0] - offset,
        indices[1] - offset, indices[2] - offset
    );
  }

  log::Logger::Log<log::Debug, log::Singleton>(
    "Read %i vertices and %i triangular facets", num_vertices, num_facets);

  return result;
}

LatticePosition barycenter(MeshData const &_mesh) {
  typedef MeshData::t_Vertices::value_type t_Vertex;
  return std::accumulate(
    _mesh.vertices.begin(), _mesh.vertices.end(), t_Vertex(0, 0, 0)
  ) / t_Vertex::value_type(_mesh.vertices.size());
}
PhysicalVolume volume(MeshData const &_mesh) {
    MeshData::t_Facets::const_iterator i_facet = _mesh.facets.begin();
    MeshData::t_Facets::const_iterator const i_facet_end(_mesh.facets.end());
    PhysicalVolume result(0);
    for(; i_facet != i_facet_end; ++i_facet) {
        LatticePosition const &v0(_mesh.vertices[(*i_facet)[0]]);
        LatticePosition const &v1(_mesh.vertices[(*i_facet)[1]]);
        LatticePosition const &v2(_mesh.vertices[(*i_facet)[2]]);
        result += v0.Cross(v1).Dot(v2);
    }
    // Minus sign comes from outward facing facet orientation
    return -result / PhysicalVolume(6);
}
PhysicalSurface surface(MeshData const &_mesh) {
    MeshData::t_Facets::const_iterator i_facet = _mesh.facets.begin();
    MeshData::t_Facets::const_iterator const i_facet_end(_mesh.facets.end());
    PhysicalVolume result(0);
    for(; i_facet != i_facet_end; ++i_facet) {
        LatticePosition const &v0(_mesh.vertices[(*i_facet)[0]]);
        LatticePosition const &v1(_mesh.vertices[(*i_facet)[1]]);
        LatticePosition const &v2(_mesh.vertices[(*i_facet)[2]]);
        result += (v0 - v1).Cross(v2 - v1).GetMagnitude();
    }
    return result * 0.5;
}


namespace {
  bool contains(MeshData::t_Facet const &_a,
          MeshData::t_Facet::value_type _v) {
    return _a[0] == _v or _a[1] == _v or _a[2] == _v;
  }
  bool edge_sharing(MeshData::t_Facet const &_a, 
      MeshData::t_Facet const &_b ) {
    return (contains(_b, _a[0]) ? 1: 0)
        + (contains(_b, _a[1]) ? 1: 0)
        + (contains(_b, _a[2]) ? 1: 0) == 2;
  }
  // Adds value as first non-negative number, if value not in array yet
  void insert(MeshData::t_Facet &_container,
    MeshData::t_Facet::value_type _value,
    MeshData::t_Facet::value_type _max) {
    for(size_t i(0); i < _container.size(); ++i)
        if(_container[i] == _max) { _container[i] = _value; return; }
        else if(_container[i] == _value) return;
  }
}

MeshTopology::MeshTopology(MeshData const &_mesh) {
  vertexToFacets.resize(_mesh.vertices.size());
  facetNeighbors.resize(_mesh.facets.size());

  // Loop over facets to create map from vertices to facets
  MeshData::t_Facets::const_iterator i_facet = _mesh.facets.begin();
  MeshData::t_Facets::const_iterator const i_facet_end = _mesh.facets.end();
  for(unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i) {
    vertexToFacets.at((*i_facet)[0]).insert(i);
    vertexToFacets.at((*i_facet)[1]).insert(i);
    vertexToFacets.at((*i_facet)[2]).insert(i);
  }

  // Now creates map of neighboring facets
  unsigned int const Nmax = _mesh.facets.size();
  boost::array<unsigned int, 3> const neg = {{ Nmax, Nmax, Nmax }};
  for(unsigned int i(0); i < facetNeighbors.size(); ++i)
    facetNeighbors[i] = neg;
  i_facet = _mesh.facets.begin();
  for(unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i) {
    for(size_t node(0); node != i_facet->size(); ++node) {
      // check facets that this node is attached to
      MeshTopology::t_VertexToFacets::const_reference
        facets = vertexToFacets.at((*i_facet)[node]);
      MeshTopology::t_VertexToFacets::value_type::const_iterator
        i_neigh = facets.begin();
      for(; i_neigh != facets.end(); ++i_neigh) {
        if(edge_sharing(*i_facet, _mesh.facets.at(*i_neigh)))
          insert(facetNeighbors.at(i), *i_neigh, Nmax);
      }
    }
  }
# ifndef NDEBUG
  // Checks there are no uninitialized values
  for(unsigned int i(0); i < facetNeighbors.size(); ++i)
    for(unsigned int j(0); j < 3; ++j)
      assert(facetNeighbors[i][j] < Nmax);
# endif
}


}} // hemelb::rbc
