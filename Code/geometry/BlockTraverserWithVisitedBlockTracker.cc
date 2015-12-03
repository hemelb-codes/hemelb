
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/BlockTraverserWithVisitedBlockTracker.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace geometry
  {
    BlockTraverserWithVisitedBlockTracker::BlockTraverserWithVisitedBlockTracker(const geometry::LatticeData& iLatDat) :
      BlockTraverser(iLatDat),
      //Initially no blocks have been visited
          mBlockVisited(iLatDat.GetBlockCount(), false)
    {
    }

    BlockTraverserWithVisitedBlockTracker::~BlockTraverserWithVisitedBlockTracker()
    {
    }

    bool BlockTraverserWithVisitedBlockTracker::GoToNextUnvisitedBlock()
    {
      do
      {
        bool validBlock = GoToNextBlock();
        if (!validBlock)
        {
          return false;
        }
      }
      while (IsCurrentBlockVisited());

      return true;
    }

    bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(size_t iN) const
    {
      return mBlockVisited[iN];
    }

    bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(util::Vector3D<site_t> iLocation) const
    {
      return IsBlockVisited(GetIndexFromLocation(iLocation));
    }

    bool BlockTraverserWithVisitedBlockTracker::IsCurrentBlockVisited() const
    {
      return IsBlockVisited(CurrentBlockNumber());
    }

    void BlockTraverserWithVisitedBlockTracker::MarkCurrentBlockVisited()
    {
      MarkBlockVisited(CurrentBlockNumber());
    }

    void BlockTraverserWithVisitedBlockTracker::MarkBlockVisited(size_t iBlockId)
    {
      mBlockVisited[iBlockId] = true;
    }

    void BlockTraverserWithVisitedBlockTracker::MarkBlockVisited(util::Vector3D<site_t> iLocation)
    {
      MarkBlockVisited(GetIndexFromLocation(iLocation));
    }
  }
}
