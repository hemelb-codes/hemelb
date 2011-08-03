#include <assert.h>

#include "geometry/LatticeData.h"
#include "vis/rayTracer/Location.h"
#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      RayTracer::ClusterBuilder::BlockIterator::BlockIterator
      (const geometry::LatticeData* iLatDat)
      {
	mLatticeData = iLatDat;
		
	//Initially no blocks have been visited
	mBlockVisited = new bool[mLatticeData->GetBlockCount()];
	for (site_t n = 0; 
	     n < mLatticeData->GetBlockCount(); 
	     n++)
	{
	  mBlockVisited[n] = false;
	}		
      }

      RayTracer::ClusterBuilder::BlockIterator::~BlockIterator()
      {
	delete[] mBlockVisited;
      }
      
      site_t RayTracer::ClusterBuilder::BlockIterator::CurrentBlockNumber()
      {
	return CurrentNumber();
      }

      Location<site_t> RayTracer::ClusterBuilder::BlockIterator::GetLowestSiteCoordinates()
      {
	return GetLocation()*mLatticeData->GetBlockSize();
      }
      		
      bool RayTracer::ClusterBuilder::BlockIterator::GoToNextUnvisitedBlock()
      {
	while (CurrentBlockVisited())
	{
	  bool validBlock = GoToNextBlock();
	  if(!validBlock)
	  {
	    return false;
	  }
	}

	return true;
      }
	    
      geometry::LatticeData::BlockData *
      RayTracer::ClusterBuilder::BlockIterator::GetBlockData()
      {
	return mLatticeData->GetBlock(mCurrentNumber);
      }

      geometry::LatticeData::BlockData *
      RayTracer::ClusterBuilder::BlockIterator::GetBlockData
      (const Location<site_t>& iLocation)
      {
	return mLatticeData->GetBlock(GetNumberFromLocation(iLocation));
      }      

      site_t RayTracer::ClusterBuilder::BlockIterator::GetBlockSize()
      {
	return mLatticeData->GetBlockSize();
      }

      RayTracer::ClusterBuilder::SiteIterator 
      RayTracer::ClusterBuilder::BlockIterator::GetSiteIterator()
      {
	return SiteIterator(mLatticeData, CurrentBlockNumber());
      }	

      RayTracer::ClusterBuilder::SiteIterator 
      RayTracer::ClusterBuilder::BlockIterator::GetSiteIterator(const Location<site_t>& iLocation)
      {
	return SiteIterator(mLatticeData, GetNumberFromLocation(iLocation));
      }	

      

      bool RayTracer::ClusterBuilder::BlockIterator::BlockValid(Location<site_t> iBlock)
      {
	return mLatticeData->IsValidBlockSite
	  (iBlock.i,
	   iBlock.j,
	   iBlock.k);
      }

      bool RayTracer::ClusterBuilder::BlockIterator::BlockVisited(site_t iN)
      {
	return mBlockVisited[iN];
      }

      bool RayTracer::ClusterBuilder::BlockIterator::BlockVisited(Location <site_t>iLocation)
      {

	return mBlockVisited[GetNumberFromLocation(iLocation)];
      }
	    

      bool RayTracer::ClusterBuilder::BlockIterator::CurrentBlockVisited()
      {
	return BlockVisited(CurrentBlockNumber());
      }

      void RayTracer::ClusterBuilder::BlockIterator::MarkBlockVisited()
      {
	MarkBlockVisited(CurrentBlockNumber());
      }
	    
      void RayTracer::ClusterBuilder::BlockIterator::MarkBlockVisited(site_t iBlockId)
      {
	mBlockVisited[iBlockId] = true;
      }

      void RayTracer::ClusterBuilder::BlockIterator::MarkBlockVisited(Location<site_t> iLocation)
      {
	site_t lNumber = GetNumberFromLocation(iLocation);
	assert(lNumber < mLatticeData->GetBlockCount());
	MarkBlockVisited(lNumber);
      }


      bool RayTracer::ClusterBuilder::BlockIterator::GoToNextBlock()
      {
	return Iterate();
      }
	  
      site_t RayTracer::ClusterBuilder::BlockIterator::GetXCount() 
      {
	return mLatticeData->GetXBlockCount();
      }

      site_t RayTracer::ClusterBuilder::BlockIterator::GetYCount()
      {
	return mLatticeData->GetYBlockCount();
      }

      site_t RayTracer::ClusterBuilder::BlockIterator::GetZCount() 
      {
	return mLatticeData->GetZBlockCount();
      }

	   
    }
  }
}
