//#define NDEBUG;
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
      RayTracer::ClusterBuilder::BlockTraverser::BlockTraverser
      (const geometry::LatticeData* iLatDat)
	: VolumeTraverser()
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

      RayTracer::ClusterBuilder::BlockTraverser::~BlockTraverser()
      {
	delete[] mBlockVisited;
      }
      
      site_t RayTracer::ClusterBuilder::BlockTraverser::CurrentBlockNumber()
      {
	return GetCurrentIndex();
      }

      Location<site_t> RayTracer::ClusterBuilder::BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
      {
	return GetCurrentLocation()*mLatticeData->GetBlockSize();
      }
      		
      bool RayTracer::ClusterBuilder::BlockTraverser::GoToNextUnvisitedBlock()
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
	    
      geometry::LatticeData::BlockData *
      RayTracer::ClusterBuilder::BlockTraverser::GetCurrentBlockData()
      {
	return mLatticeData->GetBlock(mCurrentNumber);
      }

      geometry::LatticeData::BlockData *
      RayTracer::ClusterBuilder::BlockTraverser::GetBlockDataForLocation
      (const Location<site_t>& iLocation)
      {
	return mLatticeData->GetBlock(GetIndexFromLocation(iLocation));
      }      

      site_t RayTracer::ClusterBuilder::BlockTraverser::GetBlockSize()
      {
	return mLatticeData->GetBlockSize();
      }

      RayTracer::ClusterBuilder::SiteTraverser 
      RayTracer::ClusterBuilder::BlockTraverser::GetSiteTraverserForCurrentBlock()
      {
	return SiteTraverser(mLatticeData, CurrentBlockNumber());
      }	

      RayTracer::ClusterBuilder::SiteTraverser 
      RayTracer::ClusterBuilder::BlockTraverser::GetSiteTraverserForLocation(const Location<site_t>& iLocation)
      {
	return SiteTraverser(mLatticeData, GetIndexFromLocation(iLocation));
      }	

      

      bool RayTracer::ClusterBuilder::BlockTraverser::IsValidLocation(Location<site_t> iBlock)
      {
	return mLatticeData->IsValidBlockSite
	  (iBlock.x,
	   iBlock.y,
	   iBlock.z);
      }

      bool RayTracer::ClusterBuilder::BlockTraverser::IsBlockVisited(site_t iN)
      {
	return mBlockVisited[iN];
      }

      bool RayTracer::ClusterBuilder::BlockTraverser::IsBlockVisited(Location <site_t>iLocation)
      {

	return mBlockVisited[GetIndexFromLocation(iLocation)];
      }
	    

      bool RayTracer::ClusterBuilder::BlockTraverser::IsCurrentBlockVisited()
      {
	return IsBlockVisited(CurrentBlockNumber());
      }

      void RayTracer::ClusterBuilder::BlockTraverser::MarkCurrentBlockVisited()
      {
	MarkBlockVisited(CurrentBlockNumber());
      }
	    
      void RayTracer::ClusterBuilder::BlockTraverser::MarkBlockVisited(site_t iBlockId)
      {
	mBlockVisited[iBlockId] = true;
      }

      void RayTracer::ClusterBuilder::BlockTraverser::MarkBlockVisited(Location<site_t> iLocation)
      {
	site_t lNumber = GetIndexFromLocation(iLocation);
	assert(lNumber < mLatticeData->GetBlockCount());
	MarkBlockVisited(lNumber);
      }


      bool RayTracer::ClusterBuilder::BlockTraverser::GoToNextBlock()
      {
	return TraverseOne();
      }
	  
      site_t RayTracer::ClusterBuilder::BlockTraverser::GetXCount() 
      {
	return mLatticeData->GetXBlockCount();
      }

      site_t RayTracer::ClusterBuilder::BlockTraverser::GetYCount()
      {
	return mLatticeData->GetYBlockCount();
      }

      site_t RayTracer::ClusterBuilder::BlockTraverser::GetZCount() 
      {
	return mLatticeData->GetZBlockCount();
      }

	   
    }
  }
}
