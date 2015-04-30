//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#include <iterator>
#include "redblood/buffer/Buffer.h"
#include "Exception.h"

namespace hemelb
{
  namespace redblood
  {
    namespace buffer
    {
      CellContainer::value_type Buffer::drop()
      {
        if(virtuals.size() == 0)
          throw Exception() << "No cell left to drop";

        std::vector<LatticeDistance> dists(virtuals.size(), 0);
        auto const &origin = this->geometry->origin;
        auto const &normal = this->geometry->normal;
        // normal should point from nearest to drop point to furthest from drop point
        // The drop point is *always* (0, 0, 0) in the coordinate system of the virtual buffer.
        auto getdist = [&origin, &normal](CellContainer::value_type c)
        {
          return c->GetBarycenter().Dot(normal);
        };
        std::transform(virtuals.begin(), virtuals.end(), dists.begin(), getdist);
        auto const n = std::min_element(dists.begin(), dists.end()) - dists.begin();
        auto min_iter = virtuals.begin(); std::advance(min_iter, n);
        justDropped = *min_iter;
        virtuals.erase(min_iter);

        *justDropped += offset - origin;
        return justDropped;
      }
    } // buffer
  }
} // hemelb::redblood
