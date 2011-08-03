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

      RayTracer::ClusterBuilder::SiteIterator 
      RayTracer::ClusterBuilder::BlockIterator::GetSiteIterator()
      {
	return SiteIterator(mLatticeData, CurrentBlockNumber());
      }

		
		
      site_t RayTracer::ClusterBuilder::BlockIterator::CurrentBlockNumber()
      {
	return CurrentNumber();
      }

      bool RayTracer::ClusterBuilder::BlockIterator::BlockValid(Location iBlock)
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

      bool RayTracer::ClusterBuilder::BlockIterator::BlockVisited(Location iLocation)
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
	    
      void RayTracer::ClusterBuilder::BlockIterator::MarkBlockVisited(site_t n)
      {
	mBlockVisited[CurrentBlockNumber()] = true;
      }

      void RayTracer::ClusterBuilder::BlockIterator::MarkBlockVisited(Location location)
      {
	MarkBlockVisited(mLatticeData->GetBlockIdFromBlockCoords
			 (location.i,
			  location.j,
			  location.k));
      }

      bool RayTracer::ClusterBuilder::BlockIterator::
      SitesAssignedToLocalProcessorInBlock()
      {
	if (GetBlockData()->
	    ProcessorRankForEachBlockSite == NULL)
	{
	  return false;
	}

	for (unsigned int siteId = 0;
	     siteId < mLatticeData->
	       GetSitesPerBlockVolumeUnit(); 
	     siteId++)
	     {
	       if (topology::NetworkTopology::Instance()->GetLocalRank() == 
		   GetBlockData()->ProcessorRankForEachBlockSite[siteId])
	       {
		 return true;
	       }
	     }
	       return false;
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
