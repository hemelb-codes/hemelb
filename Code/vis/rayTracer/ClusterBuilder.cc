#include <assert.h>
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
      (const geometry::LatticeData*& iLatticeData,
       std::vector<Cluster> & i_clusters,
       float **& i_cluster_voxel,
       float ***& i_cluster_flow_field
	) :
	mBlockIterator(iLatticeData),
	mClusters(i_clusters),
	mClusterVoxel(i_cluster_voxel),
	mClusterFlowField(i_cluster_flow_field),
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
  
	mClusterVoxel = new float * [mLatticeData->GetLocalFluidSiteCount() * VIS_FIELDS];

	mClusterFlowField = new float ** [mClusters.size()];

	// Process the flow field for every cluster
	for(unsigned int lThisClusterId = 0; lThisClusterId < mClusters.size(); lThisClusterId++)
	{
	  ProcessCluster(lThisClusterId);	  
	}

      }

      void RayTracer::ClusterBuilder::LocateClusters()
      {
	// Run through all unvisited blocks finding clusters
	while (mBlockIterator.GoToNextUnvisitedBlock())
	{
	  //Mark the block visited 
	  mBlockIterator.MarkBlockVisited();
	  
	  //If there are sites assigned to the local processor, search for the
	  //cluster of connected sites
	  if (SitesAssignedToLocalProcessorInBlock(mBlockIterator.GetBlockData()))
	  {
	    FindNewCluster();
	  }
	}
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
	lBlocksToProcess.push(mBlockIterator.GetLocation());

	//Loop over the cluster via neighbours until
	//all blocks have been processed
	while (!lBlocksToProcess.empty())
	{
	  //Get location off the top of the stack
	  //(we could actually take anything off the stack)
	  Location<site_t> lCurrentLocation = lBlocksToProcess.top();
	  lBlocksToProcess.pop();
	    
	  if(SitesAssignedToLocalProcessorInBlock(
	       mBlockIterator.GetBlockData(lCurrentLocation)))
	  {
	    //Update block range of the cluster
	    Location<site_t>::UpdateMinLocation(lClusterBlockMin, lCurrentLocation);
	    Location<site_t>::UpdateMaxLocation(lClusterBlockMax, lCurrentLocation);
	     
	    //Update the cluster id of the given block
	    site_t lBlockID = mBlockIterator.GetNumberFromLocation(lCurrentLocation);
	    mClusterIdOfBlock[lBlockID] = (short int) mClusters.size();

	    //Loop through all the sites on the block, to 
	    //update the site bounds on the cluster
	    SiteIterator lSiteIterator = mBlockIterator.GetSiteIterator(lCurrentLocation);
	    do
	    { 
	      //If the site is not a solid
	      if (mBlockIterator.GetBlockData(lCurrentLocation)->
		  site_data[lSiteIterator.CurrentNumber()] != BIG_NUMBER3)
	      {
		Location<site_t>::UpdateMinLocation(lClusterSiteMin,
					    lSiteIterator.GetLocation() +
					    lCurrentLocation*
					    lSiteIterator.GetBlockSize());
		
		Location<site_t>::UpdateMaxLocation(lClusterSiteMax,
					    lSiteIterator.GetLocation() +
					    lCurrentLocation*
					    lSiteIterator.GetBlockSize());
	      }   
	    }
	    while (lSiteIterator.Iterate());
			
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
	  if (mBlockIterator.BlockValid(lNeighbouringBlock))
	  {
	    //Ensure that the block hasn't been visited before
	    if(!mBlockIterator.BlockVisited(lNeighbouringBlock))
	    {
	      //Add to the stack
	      oBlocksToProcess.push(lNeighbouringBlock);

	      //We must mark this locatoin as visited so it only 
	      //gets processed once
	      mBlockIterator.MarkBlockVisited(lNeighbouringBlock);
	    }
	  }
	}
      }

      bool RayTracer::ClusterBuilder::SitesAssignedToLocalProcessorInBlock
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
	lNewCluster.minBlock.i = (float) (iClusterBlockMin.i * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetXSiteCount();
	lNewCluster.minBlock.j = (float) (iClusterBlockMin.j * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetYSiteCount();
	lNewCluster.minBlock.k = (float) (iClusterBlockMin.k * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetZSiteCount();

	lNewCluster.blocksX = static_cast<unsigned short>(1 + iClusterBlockMax.i - iClusterBlockMin.i);
	lNewCluster.blocksY = static_cast<unsigned short>(1 + iClusterBlockMax.j - iClusterBlockMin.j);
	lNewCluster.blocksZ = static_cast<unsigned short>(1 + iClusterBlockMax.k - iClusterBlockMin.k);

	//Ensure that value does not change in casting
	assert(static_cast<site_t>(lNewCluster.blocksX) == (1 + iClusterBlockMax.i - iClusterBlockMin.i));
        assert(static_cast<site_t>(lNewCluster.blocksY) == (1 + iClusterBlockMax.j - iClusterBlockMin.j));
        assert(static_cast<site_t>(lNewCluster.blocksZ) == (1 + iClusterBlockMax.k - iClusterBlockMin.k));
	
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
	   lNewCluster.minBlock.i,
	   lNewCluster.minBlock.j,
	   lNewCluster.minBlock.k,
	   lNewCluster.blocksX,
	   lNewCluster.blocksY,
	   lNewCluster.blocksZ,
	   lNewCluster.minSite.i,
	   lNewCluster.minSite.j,
	   lNewCluster.minSite.k,
	   lNewCluster.maxSite.i,
	   lNewCluster.maxSite.j,
	   lNewCluster.maxSite.k
	    );

	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Cluster min:  %u, %u, %u ", 
									       (unsigned int)iClusterBlockMin.i,
									       (unsigned int)iClusterBlockMin.j,
									       (unsigned int)iClusterBlockMin.k);


	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Cluster max:  %u, %u, %u ", 
									       (unsigned int)iClusterBlockMax.i,
									       (unsigned int)iClusterBlockMax.j,
									       (unsigned int)iClusterBlockMax.k);

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

	mClusterFlowField[iClusterId] = new float *[lCluster->blocksX * 
						      lCluster->blocksY * 
						      lCluster->blocksZ]; 

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
		(block_coordinates.i,
		 block_coordinates.j,
		 block_coordinates.k);
	      
	      Location<site_t>* mins = &mClusterBlockMins[iClusterId]; 
	      assert(block_id == 
		     ( (i + mins->i) * mLatticeData->GetYBlockCount() + (j + mins->j))
		     * mLatticeData->GetZBlockCount() + (k + mins->k) 
		);

	      
	      if (mClusterIdOfBlock[block_id] != (short int) iClusterId)
	      {
		mClusterFlowField[iClusterId][lBlockNum] = NULL;
		continue;
	      }

	      geometry::LatticeData::BlockData * lBlock = mLatticeData->GetBlock(block_id);

	      mClusterFlowField[iClusterId][lBlockNum]
		= new float[mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS];

	      UpdateFlowField(lBlock, lBlockNum, iClusterId,block_coordinates);
	    } // for k
	  } // for j
	} // for i

  
      }
      

      void RayTracer::ClusterBuilder::UpdateFlowField
      (geometry::LatticeData::BlockData * lBlock, site_t iBlockNum,  unsigned int iClusterId, Location<site_t>i_block_coordinates)
      {
	unsigned int l_site_id = -1;

	//Location site_coordinates_of_block = i_block_coordinates * mLatticeData->GetBlockSize();
	Location<site_t>siteLocOnBlock;

	for (siteLocOnBlock.i = 0; siteLocOnBlock.i < mLatticeData->GetBlockSize(); siteLocOnBlock.i++)
	{
	  for (siteLocOnBlock.j = 0; siteLocOnBlock.j < mLatticeData->GetBlockSize(); siteLocOnBlock.j++)
	  {
	    for (siteLocOnBlock.k = 0; siteLocOnBlock.k < mLatticeData->GetBlockSize(); siteLocOnBlock.k++)
	    {
	      ++l_site_id;

	      UpdateSiteFlowField(lBlock, iBlockNum, iClusterId, l_site_id);

	    }
	  }
	}
      } 
      

	     
      void RayTracer::ClusterBuilder::UpdateSiteFlowField
      ( geometry::LatticeData::BlockData * i_block,
	site_t iBlockNum, unsigned int iClusterId, unsigned int lSiteIdOnBlock)
      {
	//TODO: Clean this
      
	unsigned int lClusterVoxelSiteId = i_block->site_data[lSiteIdOnBlock];
      


	//If site is solid or not on the current processor [net.cc]
	if (lClusterVoxelSiteId == BIG_NUMBER3)
	{
	  for (site_t l = 0; l < VIS_FIELDS; l++)
	  {		
	    mClusterFlowField[iClusterId][iBlockNum][lSiteIdOnBlock * VIS_FIELDS + l] = -1.0F;
	  }
	}
	else {
	  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>
	    ("Found site correspoding to voxel site id %u", lClusterVoxelSiteId);

	  for (site_t l = 0; l < VIS_FIELDS; l++)
	  {
	    mClusterFlowField[iClusterId][iBlockNum][lSiteIdOnBlock * VIS_FIELDS + l] = 1.0F;
	  }

	  for (site_t l = 0; l < VIS_FIELDS; l++)
	  {
	   
	    mClusterVoxel[lClusterVoxelSiteId * VIS_FIELDS + l]
	      = &mClusterFlowField[iClusterId][iBlockNum][lSiteIdOnBlock * VIS_FIELDS + l];
	  }
	}
      }
    
    }
  }
}
