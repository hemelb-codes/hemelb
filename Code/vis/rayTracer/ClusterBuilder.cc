//#define NDEBUG;
#include <assert.h>
#include <map>
#include <vector>

#include <iostream>
 
#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/rayTracer/Location.h"
#include "vis/rayTracer/RayTracer.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      RayTracer::ClusterBuilder::ClusterBuilder
      (const geometry::LatticeData*& iLatticeData) :
	mBlockTraverser(iLatticeData),
	mClusterVoxelDataPointers(NUMBEROFCLUSTERVOXELMAPS),
	mLatticeData(iLatticeData)
	
      {
	//Each block is assigned a cluster id once it has been
	//assigned to a cluster
	mClusterIdOfBlock = new short int[mLatticeData->GetBlockCount()];
	for (site_t lId = 0;
	     lId < mLatticeData->GetBlockCount(); 
	     lId++)
	{
	  mClusterIdOfBlock[lId] = NOTASSIGNEDTOCLUSTER;
	}

	assert(VIS_FIELDS == 3);
      }
 
      RayTracer::ClusterBuilder::~ClusterBuilder()
      {
	delete[] mClusterIdOfBlock;
      }
    
      void RayTracer::ClusterBuilder::BuildClusters()
      {
	//Initially locate all clusters, locating their
	//range by block span and site span
	LocateClusters();
  
       	// Process the flow field for every cluster
	for(unsigned int lThisClusterId = 0; lThisClusterId < mClusters.size(); lThisClusterId++)
	{
	  ProcessCluster(lThisClusterId);	  
	}
      }

      std::vector<RayTracer::Cluster>& RayTracer::ClusterBuilder::GetClusters()
      {
	return mClusters;
      }

      RayTracer::SiteData_t* RayTracer::ClusterBuilder::GetClusterVoxelDataPointer(site_t iVoxelSiteId)
      {
	return GetDataPointerClusterVoxelSiteId(iVoxelSiteId);
      }

      void RayTracer::ClusterBuilder::LocateClusters()
      {
	// Run through all unvisited blocks finding clusters
	do
	{
	  //Mark the block visited 
	  mBlockTraverser.MarkCurrentBlockVisited();
	  
	  //If there are sites assigned to the local processor, search for the
	  //cluster of connected sites
	  if (AreSitesAssignedToLocalProcessorRankInBlock(mBlockTraverser.GetCurrentBlockData()))
	  {
	    FindNewCluster();
	  }
	}
	while (mBlockTraverser.GoToNextUnvisitedBlock());
      }	 

      void RayTracer::ClusterBuilder::FindNewCluster()
      { 
	//These locations will eventually contain the bounds of the
	//rectangular cluster, both in terms of block number and
	//site numbers
	Location<site_t> lClusterBlockMin = Location<site_t>::MaxLimit();
	Location<site_t> lClusterBlockMax = Location<site_t>::MinLimit();
	Location<site_t> lClusterSiteMin = Location<site_t>::MaxLimit();
	Location<site_t> lClusterSiteMax = Location<site_t>::MinLimit();

	//To discover the cluster, we continually visit the neighbours 
	//of sequential blocks 
	//We keep a stack of all the sites that must be processed 
	//and sequentially add neighbours to it
	std::stack<Location<site_t> > lBlocksToProcess;

	//Set up the initial condition
	lBlocksToProcess.push(mBlockTraverser.GetCurrentLocation());

	//Loop over the cluster via neighbours until
	//all blocks have been processed
	while (!lBlocksToProcess.empty())
	{
	  //Get location off the top of the stack
	  //(we could actually take anything off the stack)
	  Location<site_t> lCurrentLocation = lBlocksToProcess.top();
	  lBlocksToProcess.pop();
	    
	  if(AreSitesAssignedToLocalProcessorRankInBlock(
	       mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)))
	  {
	    //Update block range of the cluster
	    Location<site_t>::UpdateMinLocation(lClusterBlockMin, lCurrentLocation);
	    Location<site_t>::UpdateMaxLocation(lClusterBlockMax, lCurrentLocation);
	     
	    //Update the cluster id of the given block
	    site_t lBlockID = mBlockTraverser.GetIndexFromLocation(lCurrentLocation);
	    mClusterIdOfBlock[lBlockID] = (short int) mClusters.size();

	    //Loop through all the sites on the block, to 
	    //update the site bounds on the cluster
	    SiteTraverser lSiteTraverser = mBlockTraverser.GetSiteTraverserForLocation(lCurrentLocation);
	    do
	    { 
	      //If the site is not a solid
	      if (mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)->
		  site_data[lSiteTraverser.GetCurrentIndex()] != BIG_NUMBER3)
	      {
		Location<site_t>::UpdateMinLocation(lClusterSiteMin,
						    lSiteTraverser.GetCurrentLocation() +
						    lCurrentLocation*
						    mBlockTraverser.GetBlockSize());
		
		Location<site_t>::UpdateMaxLocation(lClusterSiteMax,
						    lSiteTraverser.GetCurrentLocation() +
						    lCurrentLocation*
						    mBlockTraverser.GetBlockSize());
	      }   
	    }
	    while (lSiteTraverser.TraverseOne());
			
	    //Check all the neighbouring blocks to see if they need visiting. Add them to the stack.
	    AddNeighbouringBlocks(lCurrentLocation, lBlocksToProcess);
	  }
	}
	
	AddCluster(lClusterBlockMin, lClusterBlockMax, lClusterSiteMin, lClusterSiteMax);
      }

      const Location<site_t> RayTracer::ClusterBuilder::mNeighbours[26] =
      {
	Location<site_t>(-1, -1, -1),
	Location<site_t>(-1, -1,  0),
	Location<site_t>(-1, -1,  1),
		    
	Location<site_t>(-1,  0, -1),
	Location<site_t>(-1,  0,  0),
	Location<site_t>(-1,  0,  1),

	Location<site_t>(-1,  1, -1),
	Location<site_t>(-1,  1,  0),
	Location<site_t>(-1,  1,  1),

		    
	Location<site_t>( 0, -1, -1),
	Location<site_t>( 0, -1,  0),
	Location<site_t>( 0, -1,  1),
		    
	Location<site_t>( 0,  0, -1),
	// 0 0 0 is same site
	Location<site_t>( 0,  0,  1),

	Location<site_t>( 0,  1, -1),
	Location<site_t>( 0,  1,  0),
	Location<site_t>( 0,  1,  1),

	Location<site_t>( 1, -1, -1),
	Location<site_t>( 1, -1,  0),
	Location<site_t>( 1, -1,  1),
		    
	Location<site_t>( 1,  0, -1),
	Location<site_t>( 1,  0,  0),
	Location<site_t>( 1,  0,  1),

	Location<site_t>( 1,  1, -1),
	Location<site_t>( 1,  1,  0),
	Location<site_t>( 1,  1,  1)
      };

      void RayTracer::ClusterBuilder::AddNeighbouringBlocks
      (Location<site_t> iCurrentLocation, 
       std::stack<Location<site_t> >& oBlocksToProcess)
      {
	// Loop over all neighbouring blocks
	for (int l = 0; l < 26; l++)
	{
	  Location<site_t> lNeighbouringBlock = iCurrentLocation + mNeighbours[l];
		
	  //The neighouring block location might not exist
	  //eg negative co-ordinates
	  if (mBlockTraverser.IsValidLocation(lNeighbouringBlock))
	  {
	    //Ensure that the block hasn't been visited before
	    if(!mBlockTraverser.IsBlockVisited(lNeighbouringBlock))
	    {
	      //Add to the stack
	      oBlocksToProcess.push(lNeighbouringBlock);

	      //We must mark this locatoin as visited so it only 
	      //gets processed once
	      mBlockTraverser.MarkBlockVisited(lNeighbouringBlock);
	    }
	  }
	}
      }

      bool RayTracer::ClusterBuilder::AreSitesAssignedToLocalProcessorRankInBlock
      (geometry::LatticeData::BlockData * iBlock)
      {
	if (iBlock->
	    ProcessorRankForEachBlockSite == NULL)
	{
	  return false;
	}
	
	for (unsigned int siteId = 0;
	     siteId < mLatticeData->GetSitesPerBlockVolumeUnit();
	     siteId++)
	{
	  if (topology::NetworkTopology::Instance()->GetLocalRank() ==
	      iBlock->ProcessorRankForEachBlockSite[siteId])
	  {
	    return true;
	  }
	}
	return false;
      } 

      void RayTracer::ClusterBuilder::AddCluster(Location<site_t> iClusterBlockMin,
						 Location<site_t> iClusterBlockMax,
						 Location<site_t> iClusterVoxelMin, 
						 Location<site_t> iClusterVoxelMax)
      {
	//The friendly locations must be turned into a format usable by the ray tracer
	Cluster lNewCluster;
	lNewCluster.minBlock.x = (float) (iClusterBlockMin.x * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetXSiteCount();
	lNewCluster.minBlock.y = (float) (iClusterBlockMin.y * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetYSiteCount();
	lNewCluster.minBlock.z = (float) (iClusterBlockMin.z * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetZSiteCount();

	lNewCluster.blocksX = static_cast<unsigned short>(1 + iClusterBlockMax.x - iClusterBlockMin.x);
	lNewCluster.blocksY = static_cast<unsigned short>(1 + iClusterBlockMax.y - iClusterBlockMin.y);
	lNewCluster.blocksZ = static_cast<unsigned short>(1 + iClusterBlockMax.z - iClusterBlockMin.z);

	//Ensure that value does not change in casting
	assert(static_cast<site_t>(lNewCluster.blocksX) == (1 + iClusterBlockMax.x - iClusterBlockMin.x));
        assert(static_cast<site_t>(lNewCluster.blocksY) == (1 + iClusterBlockMax.y - iClusterBlockMin.y));
        assert(static_cast<site_t>(lNewCluster.blocksZ) == (1 + iClusterBlockMax.z - iClusterBlockMin.z));
	
	lNewCluster.minSite = Location<float>(iClusterVoxelMin) -
	  Location<float>((float) mLatticeData->GetXSiteCount(),
			  (float) mLatticeData->GetYSiteCount(),
			  (float) mLatticeData->GetZSiteCount()) * 0.5F;
	
	lNewCluster.maxSite = Location<float>(iClusterVoxelMax + Location<site_t>(1)) - 
	  Location<float>(0.5F * (float) mLatticeData->GetXSiteCount(),
			  0.5F * (float) mLatticeData->GetYSiteCount(),
			  0.5F * (float) mLatticeData->GetZSiteCount());

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>
	  ("Found cluster: %f, %f, %f, %hu, %hu, %hu, %f, %f, %f, %f, %f, %f",
	   lNewCluster.minBlock.x,
	   lNewCluster.minBlock.y,
	   lNewCluster.minBlock.z,
	   lNewCluster.blocksX,
	   lNewCluster.blocksY,
	   lNewCluster.blocksZ,
	   lNewCluster.minSite.x,
	   lNewCluster.minSite.y,
	   lNewCluster.minSite.z,
	   lNewCluster.maxSite.x,
	   lNewCluster.maxSite.y,
	   lNewCluster.maxSite.z
	    );

	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Cluster min:  %u, %u, %u ", 
									      (unsigned int)iClusterBlockMin.x,
									      (unsigned int)iClusterBlockMin.y,
									      (unsigned int)iClusterBlockMin.z);


	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Cluster max:  %u, %u, %u ", 
									      (unsigned int)iClusterBlockMax.x,
									      (unsigned int)iClusterBlockMax.y,
									      (unsigned int)iClusterBlockMax.z);

	mClusters.push_back(lNewCluster);

	//We need to store the cluster block minimum in
	//order to process the cluster
	mClusterBlockMins.push_back(iClusterBlockMin);
      }


      void RayTracer::ClusterBuilder::ProcessCluster(unsigned int iClusterId)
      {
	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
	  ("Examining cluster id = %u", (unsigned int) iClusterId);


	Cluster* lCluster = &mClusters[iClusterId];
	lCluster->SiteData = std::vector<std::vector<SiteData_t> >(lCluster->blocksX * 
								   lCluster->blocksY * 
								   lCluster->blocksZ); 

	site_t lBlockNum = -1;
	for (site_t i = 0; i < lCluster->blocksX; i++)
	{
	  for (site_t j = 0; j < lCluster->blocksY; j++)
	  {
	    for (site_t k = 0; k < lCluster->blocksZ; k++)
	    {
	      ++lBlockNum;

	      hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
		("Examining block number = %u", (unsigned int) lBlockNum);

	      
	      Location<site_t>block_coordinates = Location<site_t>(i, j, k) + mClusterBlockMins[iClusterId];
	      site_t block_id = mLatticeData->GetBlockIdFromBlockCoords
		(block_coordinates.x,
		 block_coordinates.y,
		 block_coordinates.z);
	      
	      Location<site_t>* mins = &mClusterBlockMins[iClusterId]; 
	      assert(block_id == 
		     ( (i + mins->x) * mLatticeData->GetYBlockCount() + (j + mins->y))
		     * mLatticeData->GetZBlockCount() + (k + mins->z) 
		);

	      
	      if (mClusterIdOfBlock[block_id] != (short int) iClusterId)
	      {
		continue;
	      }

	      geometry::LatticeData::BlockData * lBlock = mLatticeData->GetBlock(block_id);
	      
	      //By default all values are -1 for solids
	      lCluster->SiteData[lBlockNum].resize(
		mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS, SiteData_t(-1.0F));

	      UpdateSiteData(lBlock, lBlockNum, iClusterId, block_coordinates);
	    } // for k
	  } // for j
	} // for i

  
      }
      

      void RayTracer::ClusterBuilder::UpdateSiteData
      (geometry::LatticeData::BlockData * lBlock, site_t iBlockNum,  unsigned int iClusterId, Location<site_t>i_block_coordinates)
      {
	unsigned int l_site_id = -1;

	//Location site_coordinates_of_block = i_block_coordinates * mLatticeData->GetBlockSize();
	Location<site_t>siteLocOnBlock;

	for (siteLocOnBlock.x = 0; siteLocOnBlock.x < mLatticeData->GetBlockSize(); siteLocOnBlock.x++)
	{
	  for (siteLocOnBlock.y = 0; siteLocOnBlock.y < mLatticeData->GetBlockSize(); siteLocOnBlock.y++)
	  {
	    for (siteLocOnBlock.z = 0; siteLocOnBlock.z < mLatticeData->GetBlockSize(); siteLocOnBlock.z++)
	    {
	      ++l_site_id;

	      UpdateSiteDataAtSite(lBlock, iBlockNum, iClusterId, l_site_id);

	    }
	  }
	}
      } 
      

	     
      void RayTracer::ClusterBuilder::UpdateSiteDataAtSite
      ( geometry::LatticeData::BlockData * iBlock,
	site_t iBlockNum, unsigned int iClusterId, unsigned int lSiteIdOnBlock)
      {
	//TODO: Clean this
      
	unsigned int lClusterVoxelSiteId = iBlock->site_data[lSiteIdOnBlock];

      	//If site not a solid and on  the current processor [net.cc]
	if (lClusterVoxelSiteId != BIG_NUMBER3)
	{
	  mClusters[iClusterId].SiteData[iBlockNum][lSiteIdOnBlock] = SiteData_t(1.0F);
	

	  //For efficiency we want to store a pointer to the site data grouped by the ClusterVortexID
	  //(1D organisation of sites)
	  std::vector<SiteData_t>::iterator lSiteDataIterator = 
	    mClusters[iClusterId].SiteData[iBlockNum].begin() + lSiteIdOnBlock;

	  SiteData_t* lSiteDataLocation = &(*lSiteDataIterator); 
	
	  //This one's for the C programmers out there
	  SetDataPointerForClusterVoxelSiteId
	    (lClusterVoxelSiteId, &(*lSiteDataLocation)); 
	}
      }
    

      RayTracer::SiteData_t* RayTracer::ClusterBuilder::GetDataPointerClusterVoxelSiteId(site_t iClusterVortexSiteId)
      {
	assert(mClusterVoxelDataPointers
	       [iClusterVortexSiteId%NUMBEROFCLUSTERVOXELMAPS].
	       count(iClusterVortexSiteId)==1);
	
	return mClusterVoxelDataPointers
	  [iClusterVortexSiteId%NUMBEROFCLUSTERVOXELMAPS]
	  [iClusterVortexSiteId];
      }

      void RayTracer::ClusterBuilder::SetDataPointerForClusterVoxelSiteId
      (site_t iClusterVortexSiteId, SiteData_t* iDataPointer)
      {
	assert(mClusterVoxelDataPointers
	       [iClusterVortexSiteId%NUMBEROFCLUSTERVOXELMAPS].count(iClusterVortexSiteId)==0);
	
	mClusterVoxelDataPointers
	  [iClusterVortexSiteId%NUMBEROFCLUSTERVOXELMAPS]
	  [iClusterVortexSiteId] = iDataPointer;
      }

    }
  }
}
