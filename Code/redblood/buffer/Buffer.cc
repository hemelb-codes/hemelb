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
      namespace
      {
        LatticeDistance maxCellRadius(CellContainer::value_type cell)
        {
          auto const barycenter = cell->GetBarycenter();
          auto const &vertices = cell->GetVertices();
          auto const first = vertices.begin();
          auto dist = [&barycenter](LatticePosition const &a, LatticePosition const &b)
          {
            return (a-barycenter).GetMagnitudeSquared() < (b-barycenter).GetMagnitudeSquared();
          };
          LatticePosition const &max = *std::max_element(first, vertices.end(), dist);
          return 4e0 * (max - barycenter).GetMagnitudeSquared();
        }
      }

      CellContainer::value_type Buffer::nextCell() const
      {
        if (virtuals.size() == 0)
          throw Exception() << "No cell left to drop";

        std::vector<LatticeDistance> dists(virtuals.size(), 0);
        auto const &origin = geometry->origin;
        auto const &normal = geometry->normal;
        // normal should point from nearest to drop point to furthest from drop point
        // The drop point is *always* (0, 0, 0) in the coordinate system of the virtual buffer.
        auto getdist = [&origin, &normal](CellContainer::value_type c)
        {
          return c->GetBarycenter().Dot(normal);
        };
        std::transform(virtuals.begin(), virtuals.end(), dists.begin(), getdist);
        auto const n = std::min_element(dists.begin(), dists.end()) - dists.begin();
        auto min_iter = virtuals.begin();
        std::advance(min_iter, n);
        return *min_iter;
      }

      CellContainer::value_type Buffer::drop()
      {
        auto next = nextCell();
        justDropped = next;
        virtuals.erase(next);

        lastZ = justDropped->GetBarycenter().Dot(geometry->normal);
        *justDropped += geometry->origin + geometry->normal * offset;
        return justDropped;
      }

      void Buffer::updateOffset()
      {
        // if no cell, then reset offset to zero
        if (virtuals.size() == 0)
        {
          offset = 0;
          return;
        }

        // Makes sure there is an interaction distance, otherwise compute it
        if (interactionRadius <= 0e0)
          interactionRadius = 1.25 * maxCellRadius(*virtuals.begin());

        auto const normal = geometry->normal;
        auto const zCell = normal.Dot(nextCell()->GetBarycenter());
        if (not justDropped)
        {
          offset = -zCell;
          return;
        }
        auto const zDropped = normal.Dot(justDropped->GetBarycenter() - geometry->origin) - offset;
        // dropped cell moved forward but distance with next cell is small
        // Then movement must the smallest of:
        // - distance moved by dropped cell over last LB iteration
        // - subject to the condition that the next cell does not move past the drop point.
        // The former will be negative. We've already checked that it is.
        if (zCell - interactionRadius < zDropped)
        {
          auto const delta = lastZ - zDropped;
          lastZ = zDropped;
          // Only move if dropped cell moved forward
          if (delta > 0)
          {
            offset -= std::min(delta, zCell);
          }
          return;
        }

        // Now we can move by as much as zCell, as long as we don't reach the interaction radius
        offset -= std::min(zCell, zCell - zDropped - interactionRadius);
      }
    } // buffer
  }
} // hemelb::redblood
