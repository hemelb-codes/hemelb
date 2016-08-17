//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELBSETUPTOOL_MKCGALMESH_H
#define HEMELBSETUPTOOL_MKCGALMESH_H

#include <vector>
#include "Vector.h"
#include "Cgal.h"

// Create a CGAL 3D Polyhedral surface based on the input
// Labels each face with its index in the input - use that to index the
// surface type vector
CgalMeshPtr MkCgalMesh(const std::vector<Vector>& ptsIn,
		const std::vector<Index>& polysIn);


#endif //HEMELBSETUPTOOL_MKCGALMESH_H
