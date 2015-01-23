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
  typedef std::array<size_t, 3> Facet;
  //! Facet container type
  typedef std::vector<Facet> Facets;
  //! Vertex container type
  typedef std::vector<LatticePosition> Vertices;
  //! Vertex container
  Vertices vertices;
  //! Facet container
  Facets facets;
};

LatticePosition barycenter(MeshData const &mesh);
LatticePosition barycenter(MeshData::Vertices const &vertices);
PhysicalVolume volume(MeshData const &mesh);
PhysicalVolume volume(
    MeshData::Vertices const &vertices,
    MeshData::Facets const &facets
);
PhysicalSurface surface(MeshData const &mesh);
PhysicalVolume surface(
    MeshData::Vertices const &vertices,
    MeshData::Facets const &facets
);

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
  MeshTopology(MeshData const &mesh);
  // Empty topology
  MeshTopology() {}
};

//! Triangular mesh
class Mesh {

  //! Differentiates between two shallow and deep copy constructors
  struct deepcopy_tag {};

  //! Deep copy constructor
  Mesh(Mesh const &c, deepcopy_tag const &)
    : mesh_(new MeshData(*c.mesh_)),
      topology_(new MeshTopology(*c.topology_)) {}

public:
  //! Initializes mesh from mesh data
  Mesh(std::shared_ptr<MeshData> const & mesh)
       : mesh_(mesh), topology_(new MeshTopology(*mesh)) {}
  //! Initializes mesh from mesh data and topology
  Mesh(std::shared_ptr<MeshData> const & mesh,
      std::shared_ptr<MeshTopology> const &topo)
       : mesh_(mesh), topology_(topo) {}
  //! Initialize mesh by copying data
  Mesh(MeshData const &data)
       : mesh_(new MeshData(data)), topology_(new MeshTopology(data)) {}
  //! Shallow copy constructor
  Mesh(Mesh const &in) : mesh_(in.mesh_), topology_(in.topology_) {}

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
  void SetData(MeshData const &data) { mesh_.reset(new MeshData(data)); }

  //! Makes a copy of the mesh
  Mesh clone() const { return Mesh(*this, deepcopy_tag()); }

  //! Scale mesh around barycenter
  void operator*=(Dimensionless const &scale);
  //! Scale by matrix around barycenter
  void operator*=(util::Matrix3D const &scale);
  //! Translate mesh
  void operator+=(LatticePosition const &offset);
  //! Transform mesh
  void operator+=(std::vector<LatticePosition> const &displacements);

  //! Number of nodes
  size_t GetNumberOfNodes() const { return mesh_->vertices.size(); }
  //! Const access to vertices
  MeshData::Vertices & GetVertices() { return mesh_->vertices; }
  //! Const access to vertices
  MeshData::Vertices const & GetVertices() const { return mesh_->vertices; }
  //! Const access to vertices
  MeshData::Vertices::const_reference GetVertex(size_t site) const {
    return mesh_->vertices[site];
  }
  //! Const access to facets
  MeshData::Facets const & GetFacets() const { return mesh_->facets; }

protected:
  //! Holds actual data about the mesh
  std::shared_ptr<MeshData> mesh_;
  //! Holds topology information;
  std::shared_ptr<MeshTopology> topology_;
};

//! Read mesh from file
//! Format is from T. Krueger's thesis
std::shared_ptr<MeshData> read_mesh(std::string const &filename);
//! Read mesh from file
//! Format is from T. Krueger's thesis
std::shared_ptr<MeshData> read_mesh(std::istream &stream);
//! Write mesh from file
//! Format is from T. Krueger's thesis
void write_mesh(std::ostream &stream, MeshData const &data);
//! Write mesh from file
//! Format is from T. Krueger's thesis
void write_mesh(std::string const &filename, MeshData const &data);
//! Write mesh from file in VTK XML format
void write_vtkmesh(std::ostream &stream, MeshData const &data);
//! Write mesh from file in VTK XML format
void write_vtkmesh(std::string const &filename, MeshData const &data);

//! Tetrahedron of a depth
//! Depth refers to the number of triangular subdivision in each facet
Mesh tetrahedron(unsigned int depth=0);
//! Flat triangular shape
//! Depth refers to the number of triangular subdivision in each facet
//! There are two initial facets, referring to the same three vertices.
Mesh pancakeSamosa(unsigned int depth=0);
//! Refine a mesh by decomposing each facet into four triangles
Mesh refine(Mesh mesh, unsigned int depth=0);

}}
#endif
