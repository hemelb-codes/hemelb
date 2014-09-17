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

namespace hemelb { namespace redblood {

//! Type of 3D vectors over real space
typedef util::Vector3D<double> t_Real3D;
//! Type of containers over indices
typedef std::set<unsigned int> t_Indices;

//! Holds raw mesh data
//! Data is separated into vertices and triangular facets
struct MeshData {
    //! Facet container type
    typedef std::vector<t_Indices> t_Facets;
    //! Vertex container type
    typedef std::vector<t_Real3D> t_Vertices;
    //! Vertex container
    t_Vertices vertices;
    //! Facet container
    t_Facets facets;
};

//! Holds data about facet-facet connections, e.g. neighbors
struct MeshTopology {
    //! For each vertex, lists the facet indices
    std::vector< std::set<unsigned int> > vertex_to_facets;
    //! For each facet, lists the neighboring facets
    std::vector< boost::array<unsigned int, 3> > facet_neighbors;

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
    t_Real3D barycenter() const;
    //! Connectivity data
    boost::shared_ptr<const MeshTopology> connectivity() const
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
