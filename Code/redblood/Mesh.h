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

#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <vector>
#include <set>
#include <string>
#include <array>
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb { namespace redblood {


//! Holds raw mesh data
//! Data is separated into vertices and triangular facets
struct MeshData {
  //! Type of containers over indices
  typedef boost::array<size_t, 3> t_Facet;
  //! Facet container type
  typedef std::vector<t_Facet> t_Facets;
  //! Vertex container type
  typedef std::vector<LatticePosition> t_Vertices;
  //! Vertex container
  t_Vertices vertices;
  //! Facet container
  t_Facets facets;
};

//! Holds data about facet-facet connections, e.g. neighbors
struct MeshTopology {
  //! Type for map from vertices to facets
  typedef std::vector< std::set<size_t> > t_VertexToFacets;
  //! Type for map from facets to its neighbors
  typedef std::vector< boost::array<size_t, 3> > t_FacetNeighbors;
  //! For each vertex, lists the facet indices
  t_VertexToFacets vertexToFacets;
  //! For each facet, lists the neighboring facets
  t_FacetNeighbors facetNeighbors;

  // Creates mesh topology from mesh data
  MeshTopology(MeshData const &_mesh);
};

//! Performs 
class Mesh {

public:
  //! Initializes mesh from mesh data
  Mesh   (boost::shared_ptr<MeshData> const & _mesh)
       : mesh_(_mesh), connectivity_(new MeshTopology(*_mesh)) {}

  //! Determines barycenter of mesh
  LatticePosition GetBarycenter() const;
  //! Connectivity data
  boost::shared_ptr<const MeshTopology> GetConnectivity() const
    { return connectivity_; }

protected:
  //! Holds actual data about the mesh
  boost::shared_ptr<MeshData> mesh_;
  //! Holds connectivity information;
  boost::shared_ptr<MeshTopology> connectivity_;
};

//! Read mesh from file
//! Format is from T. Krueger's thesis
boost::shared_ptr<MeshData> read_mesh(std::string const &_filename);


}}
#endif
