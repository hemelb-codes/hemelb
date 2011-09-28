#include "geometry/BlockTraverserWithVisitedBlockTracker.h"
#include "geometry/LatticeData.h"
#include "util/Vector3D.h"
#include "vis/rayTracer/RayTracer.h"

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

    bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(size_t iN)
    {
      return mBlockVisited[iN];
    }

    bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(util::Vector3D<site_t> iLocation)
    {

      return mBlockVisited[GetIndexFromLocation(iLocation)];
    }

    bool BlockTraverserWithVisitedBlockTracker::IsCurrentBlockVisited()
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
      site_t lNumber = GetIndexFromLocation(iLocation);
      MarkBlockVisited(lNumber);
    }
  }
}
