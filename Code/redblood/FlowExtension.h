//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_FLOWEXTENSION_H
#define HEMELB_REDBLOOD_FLOWEXTENSION_H

#include <redblood/Mesh.h>
#include "util/Vector3D.h"
#include "units.h"

namespace hemelb
{
  namespace redblood
  {

    struct FlowExtension
    {
      util::Vector3D<Dimensionless> normal;
      LatticePosition position;
      LatticeDistance length, radius;
    };

    //! Checks whether a cell is inside a flow extension
    bool contains(const FlowExtension &, const MeshData::Vertices::value_type &);

  } // namespace hemelb::redblood
} // namespace hemelb

#endif
