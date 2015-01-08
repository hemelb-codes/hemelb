// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_REDBLOOD_MESH_H
#define HEMELB_REDBLOOD_MESH_H

#include <memory>
#include <array>
#include <vector>
#include <set>
#include <string>
#include "util/Vector3D.h"
#include "util/Matrix3D.h"
#include "units.h"

namespace hemelb { namespace redblood {


//! Holds raw mesh data
//! Data is separated into vertices and triangular facets
struct MeshData {
  //! Type of containers over indices
  typedef std::array<size_t, 3> t_Facet;
  //! Facet container type
  typedef std::vector<t_Facet> t_Facets;
  //! Vertex container type
  typedef std::vector<LatticePosition> t_Vertices;
  //! Vertex container
  t_Vertices vertices;
  //! Facet container
  t_Facets facets;
};

LatticePosition barycenter(MeshData const &_mesh);
PhysicalVolume volume(MeshData const &_mesh);
PhysicalSurface surface(MeshData const &_mesh);

//! Holds raw topology data
struct MeshTopology {
  //! Type for map from vertices to facets
  typedef std::vector< std::set<size_t> > t_VertexToFacets;
  //! Type for map from facets to its neighbors
  typedef std::vector< std::array<size_t, 3> > t_FacetNeighbors;
  //! For each vertex, lists the facet indices
  t_VertexToFacets vertexToFacets;
  //! For each facet, lists the neighboring facets
  t_FacetNeighbors facetNeighbors;

  // Creates mesh topology from mesh data
  MeshTopology(MeshData const &_mesh);
  // Empty topology
  MeshTopology() {}
};

//! Triangular mesh
class Mesh {

  //! Differentiates between two shallow and deep copy constructors
  struct deepcopy_tag {};

  //! Deep copy constructor
  Mesh(Mesh const &_c, deepcopy_tag const &)
    : mesh_(new MeshData(*_c.mesh_)),
      topology_(new MeshTopology(*_c.topology_)) {}

public:
  //! Initializes mesh from mesh data
  Mesh(std::shared_ptr<MeshData> const & _mesh)
       : mesh_(_mesh), topology_(new MeshTopology(*_mesh)) {}
  //! Initializes mesh from mesh data and topology
  Mesh(std::shared_ptr<MeshData> const & _mesh,
      std::shared_ptr<MeshTopology> const &_topo)
       : mesh_(_mesh), topology_(_topo) {}
  //! Initialize mesh by copying data
  Mesh(MeshData const &_data)
       : mesh_(new MeshData(_data)), topology_(new MeshTopology(_data)) {}
  //! Shallow copy constructor
  Mesh(Mesh const &_in) : mesh_(_in.mesh_), topology_(_in.topology_) {}

  //! Determines barycenter of mesh
  LatticePosition GetBarycenter() const { return barycenter(*mesh_); }
  //! Computes volume of the mesh
  PhysicalVolume GetVolume() const { return volume(*mesh_); }
  //! Computes surface of the mesh
  PhysicalSurface GetSurface() const { return surface(*mesh_); }
  //! Connectivity data
  std::shared_ptr<const MeshTopology> GetTopology() const
    { return topology_; }
  //! Returns mesh data
  std::shared_ptr<MeshData> GetData() { return mesh_; }
  //! Returns mesh data
  std::shared_ptr<MeshData const> GetData() const { return mesh_; }
  //! Returns mesh data
  void SetData(MeshData const &_data) { mesh_.reset(new MeshData(_data)); }

  //! Makes a copy of the mesh
  Mesh clone() const { return Mesh(*this, deepcopy_tag()); }

  //! Scale mesh around barycenter
  void operator*=(Dimensionless const &_scale);
  //! Scale by matrix around barycenter
  void operator*=(util::Matrix3D const &_scale);
  //! Translate mesh
  void operator+=(LatticePosition const &_offset);
  //! Transform mesh
  void operator+=(std::vector<LatticePosition> const &_displacements);

  //! Number of nodes
  size_t GetNumberOfNodes() const { return mesh_->vertices.size(); }
  //! Const access to vertices
  MeshData::t_Vertices & GetVertices() { return mesh_->vertices; }
  //! Const access to vertices
  MeshData::t_Vertices const & GetVertices() const { return mesh_->vertices; }
  //! Const access to vertices
  MeshData::t_Vertices::const_reference GetVertex(size_t _site) const {
    return mesh_->vertices[_site];
  }
  //! Const access to facets
  MeshData::t_Facets const & GetFacets() const { return mesh_->facets; }

protected:
  //! Holds actual data about the mesh
  std::shared_ptr<MeshData> mesh_;
  //! Holds topology information;
  std::shared_ptr<MeshTopology> topology_;
};

//! Read mesh from file
//! Format is from T. Krueger's thesis
std::shared_ptr<MeshData> read_mesh(std::string const &_filename);
//! Read mesh from file
//! Format is from T. Krueger's thesis
std::shared_ptr<MeshData> read_mesh(std::istream &_stream);
//! Write mesh from file
//! Format is from T. Krueger's thesis
void write_mesh(std::ostream &_stream, MeshData const &_data);
//! Write mesh from file
//! Format is from T. Krueger's thesis
void write_mesh(std::string const &_filename, MeshData const &_data);
//! Write mesh from file in VTK XML format
void write_vtkmesh(std::ostream &_stream, MeshData const &_data);
//! Write mesh from file in VTK XML format
void write_vtkmesh(std::string const &_filename, MeshData const &_data);

//! Tetrahedron of a depth
//! Depth refers to the number of triangular subdivision in each facet
Mesh tetrahedron(unsigned int depth=0);
//! Flat triangular shape
//! Depth refers to the number of triangular subdivision in each facet
//! There are two initial facets, referring to the same three vertices.
Mesh pancakeSamosa(unsigned int depth=0);
//! Refine a mesh by decomposing each facet into four triangles
Mesh refine(Mesh _mesh, unsigned int depth=0);

}}
#endif
