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
#include <vector>
#include <string>
#include "util/Vector3D.h"

namespace hemelb { namespace redblood {

//! Holds raw mesh data
//! Data is separated into vertices and triangular facets
struct MeshData {
    std::vector< util::Vector3D<double> > vertices;
    std::vector< util::Vector3D<unsigned int> > facets;
};

class Mesh {

    public:
    Mesh(boost::shared_ptr<MeshData> const & _data) : data_(_data) {};

    protected:
    //! Holds actual data
    boost::shared_ptr<MeshData> data_;
};

//! Read mesh from file
//! Format is from T. Krueger's thesis
boost::shared_ptr<MeshData> read_mesh(std::string const &_filename);


}}
#endif
