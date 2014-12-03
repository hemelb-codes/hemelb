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
#include <map>

#include "redblood/Mesh.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "Exception.h"
#include "Constants.h"

namespace hemelb { namespace redblood {

boost::shared_ptr<MeshData> read_mesh(std::string const &_filename) {
  log::Logger::Log<log::Debug, log::Singleton>(
          "Reading red blood cell from %s", _filename.c_str());

  // Open file if it exists
  std::ifstream file;
  if (!util::file_exists(_filename.c_str()))
    throw Exception() << "Red-blood-cell mesh file '"
      << _filename.c_str() << "' does not exist";
  file.open(_filename.c_str());
  return read_mesh(file);
}

boost::shared_ptr<MeshData> read_mesh(std::istream &_stream) {
  log::Logger::Log<log::Debug, log::Singleton>(
          "Reading red blood cell from stream");

  std::string line;
  // Drop header
  for(int i(0); i < 4; ++i)
    std::getline(_stream, line);

  // Number of vertices
  unsigned int num_vertices;
  _stream >> num_vertices;


  // Create Mesh data
  boost::shared_ptr<MeshData> result(new MeshData);
  result->vertices.resize(num_vertices);

  // Then read in first and subsequent lines
  MeshData::t_Facet::value_type offset;
  _stream >> offset >> result->vertices[0].x
    >> result->vertices[0].y
    >> result->vertices[0].z;
  log::Logger::Log<log::Trace, log::Singleton>(
    "Vertex 0 at %d, %d, %d", result->vertices[0].x,
    result->vertices[0].y, result->vertices[0].z
  );
  for(unsigned int i(1), index(0); i < num_vertices; ++i) {
    _stream >> index >> result->vertices[i].x
      >> result->vertices[i].y
      >> result->vertices[i].z;
    log::Logger::Log<log::Trace, log::Singleton>(
      "Vertex %i at %d, %d, %d", i, result->vertices[i].x,
      result->vertices[i].y, result->vertices[i].z
    );
  }

  // Drop mid-file headers
  for(int i(0); i < 3; ++i) std::getline(_stream, line);

  // Read facet indices
  unsigned int num_facets;
  MeshData::t_Facet indices;
  _stream >> num_facets;
  result->facets.resize(num_facets);
  for(unsigned int i(0), dummy(0); i < num_facets; ++i) {
    _stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
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

void write_mesh(
    std::string const &_filename, MeshData const & _data) {
  log::Logger::Log<log::Debug, log::Singleton>(
          "Writing red blood cell from %s", _filename.c_str());
  std::ofstream file(_filename.c_str());
  write_mesh(file, _data);
}

void write_vtkmesh(
    std::string const &_filename, MeshData const & _data) {
  log::Logger::Log<log::Debug, log::Singleton>(
          "Writing red blood cell from %s", _filename.c_str());
  std::ofstream file(_filename.c_str());
  write_vtkmesh(file, _data);
}

void write_mesh(std::ostream &_stream, MeshData const &_data) {
  // Write Header
  _stream << "$MeshFormat\n2 0 8\n$EndMeshFormat\n"
    << "$Nodes\n"
    << _data.vertices.size() << "\n";

  typedef MeshData::t_Vertices::const_iterator VertexIterator;
  VertexIterator i_vertex = _data.vertices.begin();
  VertexIterator const i_vertex_end = _data.vertices.end();
  for(unsigned i(1); i_vertex != i_vertex_end; ++i_vertex, ++i)
    _stream << i << " "
      << (*i_vertex)[0] << " "
      << (*i_vertex)[1] << " "
      << (*i_vertex)[2] << "\n";
  _stream << "$EndNode\n"
       << "$Elements\n"
       << _data.facets.size() << "\n";

  typedef MeshData::t_Facets::const_iterator FacetIterator;
  FacetIterator i_facet = _data.facets.begin();
  FacetIterator const i_facet_end = _data.facets.end();
  for(unsigned i(1); i_facet != i_facet_end; ++i_facet, ++i)
    _stream << i << " 1 2 3 4 5 "
      << (*i_facet)[0] + 1 << " "
      << (*i_facet)[1] + 1 << " "
      << (*i_facet)[2] + 1 << "\n";
  _stream << "$EndElement\n";
}

void write_vtkmesh(std::ostream &_stream, MeshData const &_data) {
  // Write Header
  _stream << "<?xml version=\"1.0\"?>\n"
    << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
    << "  <PolyData>\n"
    << "    <Piece NumberOfPoints=\""
    << _data.vertices.size() << "\" NumberOfPolys=\""
    << _data.facets.size() << "\">\n"
    << "      <Points>\n"
    << "        <DataArray NumberOfComponents=\"3\" type=\"Float32\">\n";

  typedef MeshData::t_Vertices::const_iterator VertexIterator;
  VertexIterator i_vertex = _data.vertices.begin();
  VertexIterator const i_vertex_end = _data.vertices.end();
  for(unsigned i(1); i_vertex != i_vertex_end; ++i_vertex, ++i)
    _stream << (*i_vertex)[0] << " "
      << (*i_vertex)[1] << " "
      << (*i_vertex)[2] << " ";
  _stream << "\n        </DataArray>\n" << "      </Points>\n"
    << "      <Polys>\n"
    << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";

  typedef MeshData::t_Facets::const_iterator FacetIterator;
  FacetIterator i_facet = _data.facets.begin();
  FacetIterator const i_facet_end = _data.facets.end();
  for(; i_facet != i_facet_end; ++i_facet)
    _stream << (*i_facet)[0] << " " << (*i_facet)[1] << " "
      << (*i_facet)[2] << " ";
  _stream << "\n        </DataArray>\n"
    << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";
  for(unsigned i(0); i < _data.facets.size(); ++i)
    _stream << (i+1) * 3 << " ";
  _stream << "\n        </DataArray>\n" << "      </Polys>\n";
  _stream << "    </Piece>\n  </PolyData>\n</VTKFile>\n";
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
        + (contains(_b, _a[2]) ? 1: 0) >= 2;
  }
  // Adds value as first non-negative number, if value not in array yet
  void insert(MeshData::t_Facet &_container,
    MeshData::t_Facet::value_type _value,
    MeshData::t_Facet::value_type _max) {
    for(size_t i(0); i < _container.size(); ++i)
        if(_container[i] >= _max) { _container[i] = _value; return; }
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
  boost::array<unsigned int, 3> const neg = {{Nmax, Nmax,  Nmax}};
  std::fill(facetNeighbors.begin(), facetNeighbors.end(), neg);
  i_facet = _mesh.facets.begin();
  for(unsigned int i(0); i_facet != i_facet_end; ++i_facet, ++i) {
    for(size_t node(0); node != i_facet->size(); ++node) {
      // check facets that this node is attached to
      MeshTopology::t_VertexToFacets::const_reference
        neighboringFacets = vertexToFacets.at((*i_facet)[node]);
      MeshTopology::t_VertexToFacets::value_type::const_iterator
        i_neigh = neighboringFacets.begin();
      for(; i_neigh != neighboringFacets.end(); ++i_neigh) {
        if(i == *i_neigh) continue;
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

void Mesh::operator*=(Dimensionless const &_scale) {
  LatticePosition const barycenter = GetBarycenter();
  MeshData::t_Vertices::iterator i_first = mesh_->vertices.begin();
  MeshData::t_Vertices::iterator const i_end = mesh_->vertices.end();
  for(; i_first != i_end; ++i_first)
    (*i_first) = (*i_first - barycenter) * _scale + barycenter;
}

void Mesh::operator+=(LatticePosition const &_offset) {
  MeshData::t_Vertices::iterator i_first = mesh_->vertices.begin();
  MeshData::t_Vertices::iterator const i_end = mesh_->vertices.end();
  for(; i_first != i_end; ++i_first)
    (*i_first) += _offset;
}

namespace {
  boost::shared_ptr<MeshData> initial_tetrahedron() {
    boost::shared_ptr<MeshData> data(new MeshData);

    // facets at something degrees from one another
    data->vertices.push_back(LatticePosition(0, 0, 0));
    data->vertices.push_back(LatticePosition(1, 0, 1));
    data->vertices.push_back(LatticePosition(1, 1, 0));
    data->vertices.push_back(LatticePosition(0, 1, 1));

    redblood::MeshData::t_Facet indices;
    indices[0] = 0; indices[1] = 1; indices[2] = 2;
    data->facets.push_back(indices);
    indices[0] = 0; indices[1] = 2; indices[2] = 3;
    data->facets.push_back(indices);
    indices[0] = 0; indices[1] = 3; indices[2] = 1;
    data->facets.push_back(indices);
    indices[0] = 1; indices[1] = 3; indices[2] = 2;
    data->facets.push_back(indices);
    return data;
  }

  size_t vertex(boost::shared_ptr<MeshData> &_data,
      std::map<std::pair<size_t, size_t>, size_t> &_vertices,
      size_t const &_i0, size_t const &_i1) {
    std::pair<size_t, size_t> const indices(_i0, _i1);
    std::map<std::pair<size_t, size_t>, size_t>:: const_iterator i_found
      = _vertices.find(indices);
    if(i_found == _vertices.end())
      i_found = _vertices.find(std::pair<size_t, size_t>(_i1, _i0));
    if(i_found == _vertices.end()) {
      _data->vertices.push_back(
          (_data->vertices[_i0] + _data->vertices[_i1]) * 0.5);
      _vertices[indices] = _data->vertices.size() - 1;
      return _data->vertices.size() - 1;
    }
    return i_found->second;
  }

  void refine(boost::shared_ptr<MeshData> &_data) {
    MeshData::t_Facets const facets(_data->facets);
    _data->facets.clear();
    _data->facets.resize(facets.size() * 4);
    _data->vertices.reserve(_data->facets.size() * 3 + _data->vertices.size());

    // Container with midpoint indices, so midpoints are only added once
    typedef std::pair<size_t, size_t> t_Pair;
    std::map<t_Pair, size_t> new_vertices;

    MeshData::t_Facets::const_iterator i_orig_facet(facets.begin());
    MeshData::t_Facets::const_iterator const i_orig_facet_end(facets.end());
    MeshData::t_Facets::iterator i_facet = _data->facets.begin();
    for(; i_orig_facet != i_orig_facet_end; ++i_orig_facet) {
      MeshData::t_Facet::value_type const i0 = (*i_orig_facet)[0];
      MeshData::t_Facet::value_type const i1 = (*i_orig_facet)[1];
      MeshData::t_Facet::value_type const i2 = (*i_orig_facet)[2];

      // Adds new vertices halfway through edges
      MeshData::t_Facet::value_type
        mid0(vertex(_data, new_vertices, i0, i1)),
        mid1(vertex(_data, new_vertices, i1, i2)),
        mid2(vertex(_data, new_vertices, i2, i0));

      // Adds all four new faces
      (*i_facet)[0] = i0;
      (*i_facet)[1] = mid0;
      (*i_facet)[2] = mid2;
      ++i_facet;

      (*i_facet)[0] = mid0;
      (*i_facet)[1] = i1;
      (*i_facet)[2] = mid1;
      ++i_facet;

      (*i_facet)[0] = mid1;
      (*i_facet)[1] = i2;
      (*i_facet)[2] = mid2;
      ++i_facet;

      (*i_facet)[0] = mid0;
      (*i_facet)[1] = mid1;
      (*i_facet)[2] = mid2;
      ++i_facet;
    }
  }
}

Mesh refine(Mesh _data, unsigned int _depth) {
  if(_depth == 0) return _data.clone();
  boost::shared_ptr<MeshData> data(new MeshData(*_data.GetData()));
  for(unsigned int i(0); i < _depth; ++i)
    refine(data);
  return Mesh(data);
}


Mesh tetrahedron(unsigned int depth) {
  boost::shared_ptr<MeshData> result(initial_tetrahedron());
  for(unsigned int i(0); i < depth; ++i)
    refine(result);
  return Mesh(result);
}

Mesh pancakeSamosa(unsigned int depth) {
  boost::shared_ptr<redblood::MeshData> mesh(new redblood::MeshData);

  // facets at something degrees from one another
  mesh->vertices.push_back(LatticePosition(0, 0, 0));
  mesh->vertices.push_back(LatticePosition(1, 0, 1));
  mesh->vertices.push_back(LatticePosition(1, 1, 0));

  redblood::MeshData::t_Facet indices;
  indices[0] = 0; indices[1] = 1; indices[2] = 2;
  mesh->facets.push_back(indices);
  indices[0] = 2; indices[1] = 1; indices[2] = 0;
  mesh->facets.push_back(indices);

  // Create topology by hand cos we generally don't allow for this kind of
  // ambiguous self-referencing shape.
  boost::shared_ptr<redblood::MeshTopology> topo(new redblood::MeshTopology);
  MeshTopology::t_VertexToFacets::value_type v2f;
  v2f.insert(0); v2f.insert(1);
  topo->vertexToFacets.resize(3, v2f);

  MeshTopology::t_FacetNeighbors::value_type neighbors[2]
    = {{{0, 0, 0}}, {{1, 1, 1}}};
  topo->facetNeighbors.push_back(neighbors[1]);
  topo->facetNeighbors.push_back(neighbors[0]);

  return refine(Mesh(mesh, topo), depth);
}

}} // hemelb::rbc
