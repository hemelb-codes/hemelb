//#define NDEBUG;
#include <assert.h>

#include "geometry/LatticeData.h"
#include "vis/Vector3D.h"
#include "vis/rayTracer/BlockTraverserWithVisitedBlockTracker.h"
#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      BlockTraverserWithVisitedBlockTracker::BlockTraverserWithVisitedBlockTracker
      (const geometry::LatticeData& iLatDat)
	: BlockTraverser(iLatDat)
      {
	//Initially no blocks have been visited
	mBlockVisited = new bool[mLatticeData.GetBlockCount()];
	for (site_t n = 0; 
	     n < mLatticeData.GetBlockCount(); 
	     n++)
	{
	  mBlockVisited[n] = false;
	}		
      }

      BlockTraverserWithVisitedBlockTracker::~BlockTraverserWithVisitedBlockTracker()
      {
	delete[] mBlockVisited;
      }
      
           		
      bool BlockTraverserWithVisitedBlockTracker::GoToNextUnvisitedBlock()
      {
	assert(IsCurrentBlockVisited());
	do 
	{
	  bool validBlock = GoToNextBlock();
	  if(!validBlock)
	  {
	    return false;
	  }
	}
	while (IsCurrentBlockVisited());

	return true;
      }

      bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(site_t iN)
      {
	return mBlockVisited[iN];
      }

      bool BlockTraverserWithVisitedBlockTracker::IsBlockVisited(Vector3D<site_t>iLocation)
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
	    
      void BlockTraverserWithVisitedBlockTracker::MarkBlockVisited(site_t iBlockId)
      {
	mBlockVisited[iBlockId] = true;
      }

      void BlockTraverserWithVisitedBlockTracker::MarkBlockVisited(Vector3D<site_t> iLocation)
      {
	site_t lNumber = GetIndexFromLocation(iLocation);
	assert(lNumber < mLatticeData.GetBlockCount());
	MarkBlockVisited(lNumber);
      }	   
    }
  }
}
