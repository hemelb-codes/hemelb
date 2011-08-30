//#define NDEBUG
#include <assert.h>
#include <map>
#include <vector>

#include <iostream>
 
#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "vis/Vector3D.h"
#include "vis/rayTracer/BlockTraverser.h"
#include "vis/rayTracer/ClusterBuilder.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/rayTracer/SiteData.h"
#include "vis/rayTracer/SiteTraverser.h"
#include "vis/rayTracer/VolumeTraverser.h"
#include "log/Logger.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterBuilder::ClusterBuilder
      (const geometry::LatticeData*& iLatticeData) :
	mLatticeData(iLatticeData),
	mBlockTraverser(iLatticeData)      
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
 
      ClusterBuilder::~ClusterBuilder()
      {
	delete[] mClusterIdOfBlock;

	for (std::vector<Cluster*>::iterator it = mClusters.begin(); it != mClusters.end(); it++ )
	{
	  delete *it;
	}
      }
    
      void ClusterBuilder::BuildClusters()
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

      std::vector<Cluster*>& ClusterBuilder::GetClusters()
      {
	return mClusters;
      }

      SiteData_t* ClusterBuilder::GetClusterVoxelDataPointer(site_t iVoxelSiteId)
      {
	return GetDataPointerClusterVoxelSiteId(iVoxelSiteId);
      }

      void ClusterBuilder::LocateClusters()
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

      void ClusterBuilder::FindNewCluster()
      { 
	//These locations will eventually contain the bounds of the
	//rectangular cluster, both in terms of block number and
	//site numbers
	Vector3D<site_t> lClusterBlockMin = Vector3D<site_t>::MaxLimit();
	Vector3D<site_t> lClusterBlockMax = Vector3D<site_t>::MinLimit();
	Vector3D<site_t> lClusterSiteMin = Vector3D<site_t>::MaxLimit();
	Vector3D<site_t> lClusterSiteMax = Vector3D<site_t>::MinLimit();

	//To discover the cluster, we continually visit the neighbours 
	//of sequential blocks 
	//We keep a stack of all the sites that must be processed 
	//and sequentially add neighbours to it
	std::stack<Vector3D<site_t> > lBlocksToProcess;

	//Set up the initial condition
	lBlocksToProcess.push(mBlockTraverser.GetCurrentLocation());

	//Loop over the cluster via neighbours until
	//all blocks have been processed
	while (!lBlocksToProcess.empty())
	{
	  //Get location off the top of the stack
	  //(we could actually take anything off the stack)
	  Vector3D<site_t> lCurrentLocation = lBlocksToProcess.top();
	  lBlocksToProcess.pop();
	    
	  if(AreSitesAssignedToLocalProcessorRankInBlock(
	       mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)))
	  {
	    //Update block range of the cluster
	    Vector3D<site_t>::UpdateMinVector3D(lClusterBlockMin, lCurrentLocation);
	    Vector3D<site_t>::UpdateMaxVector3D(lClusterBlockMax, lCurrentLocation);
	     
	    //Update the cluster id of the given block
	    site_t lBlockID = mBlockTraverser.GetIndexFromLocation(lCurrentLocation);
	    mClusterIdOfBlock[lBlockID] = (short int) mClusters.size();

	    //Loop through all the sites on the block, to 
	    //update the site bounds on the cluster
	    SiteTraverser lSiteTraverser = mBlockTraverser.GetSiteTraverser();
	    do
	    { 
	      //If the site is not a solid
	      if (mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)->
		  site_data[lSiteTraverser.GetCurrentIndex()] != BIG_NUMBER3)
	      {
		Vector3D<site_t>::UpdateMinVector3D(lClusterSiteMin,
						    lSiteTraverser.GetCurrentLocation() +
						    lCurrentLocation*
						    mBlockTraverser.GetBlockSize());
		
		Vector3D<site_t>::UpdateMaxVector3D(lClusterSiteMax,
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

      const Vector3D<site_t> ClusterBuilder::mNeighbours[26] =
      {
	Vector3D<site_t>(-1, -1, -1),
	Vector3D<site_t>(-1, -1,  0),
	Vector3D<site_t>(-1, -1,  1),
		    
	Vector3D<site_t>(-1,  0, -1),
	Vector3D<site_t>(-1,  0,  0),
	Vector3D<site_t>(-1,  0,  1),

	Vector3D<site_t>(-1,  1, -1),
	Vector3D<site_t>(-1,  1,  0),
	Vector3D<site_t>(-1,  1,  1),

		    
	Vector3D<site_t>( 0, -1, -1),
	Vector3D<site_t>( 0, -1,  0),
	Vector3D<site_t>( 0, -1,  1),
		    
	Vector3D<site_t>( 0,  0, -1),
	// 0 0 0 is same site
	Vector3D<site_t>( 0,  0,  1),

	Vector3D<site_t>( 0,  1, -1),
	Vector3D<site_t>( 0,  1,  0),
	Vector3D<site_t>( 0,  1,  1),

	Vector3D<site_t>( 1, -1, -1),
	Vector3D<site_t>( 1, -1,  0),
	Vector3D<site_t>( 1, -1,  1),
		    
	Vector3D<site_t>( 1,  0, -1),
	Vector3D<site_t>( 1,  0,  0),
	Vector3D<site_t>( 1,  0,  1),

	Vector3D<site_t>( 1,  1, -1),
	Vector3D<site_t>( 1,  1,  0),
	Vector3D<site_t>( 1,  1,  1)
      };

      void ClusterBuilder::AddNeighbouringBlocks
      (Vector3D<site_t> iCurrentLocation, 
       std::stack<Vector3D<site_t> >& oBlocksToProcess)
      {
	// Loop over all neighbouring blocks
	for (int l = 0; l < 26; l++)
	{
	  Vector3D<site_t> lNeighbouringBlock = iCurrentLocation + mNeighbours[l];
		
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

      bool ClusterBuilder::AreSitesAssignedToLocalProcessorRankInBlock
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

      Cluster* ClusterBuilder::CreateNewCluster()
      {
	return new Cluster();
      }

      void ClusterBuilder::AddCluster(Vector3D<site_t> iClusterBlockMin,
				      Vector3D<site_t> iClusterBlockMax,
				      Vector3D<site_t> iClusterVoxelMin, 
				      Vector3D<site_t> iClusterVoxelMax)
      {
	//The friendly locations must be turned into a format usable by the ray tracer
	Cluster* lNewCluster = CreateNewCluster();
	lNewCluster->minBlock.x = (float) (iClusterBlockMin.x * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetXSiteCount();
	lNewCluster->minBlock.y = (float) (iClusterBlockMin.y * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetYSiteCount();
	lNewCluster->minBlock.z = (float) (iClusterBlockMin.z * mLatticeData->GetBlockSize()) - 0.5F
	  * (float) mLatticeData->GetZSiteCount();

	lNewCluster->blocksX = static_cast<unsigned short>(1 + iClusterBlockMax.x - iClusterBlockMin.x);
	lNewCluster->blocksY = static_cast<unsigned short>(1 + iClusterBlockMax.y - iClusterBlockMin.y);
	lNewCluster->blocksZ = static_cast<unsigned short>(1 + iClusterBlockMax.z - iClusterBlockMin.z);

	//Ensure that value does not change in casting
	assert(static_cast<site_t>(lNewCluster->blocksX) == (1 + iClusterBlockMax.x - iClusterBlockMin.x));
        assert(static_cast<site_t>(lNewCluster->blocksY) == (1 + iClusterBlockMax.y - iClusterBlockMin.y));
        assert(static_cast<site_t>(lNewCluster->blocksZ) == (1 + iClusterBlockMax.z - iClusterBlockMin.z));
	
	lNewCluster->minSite = Vector3D<float>(iClusterVoxelMin) -
	  Vector3D<float>((float) mLatticeData->GetXSiteCount(),
			  (float) mLatticeData->GetYSiteCount(),
			  (float) mLatticeData->GetZSiteCount()) * 0.5F;
	
	lNewCluster->maxSite = Vector3D<float>(iClusterVoxelMax + Vector3D<site_t>(1)) - 
	  Vector3D<float>(0.5F * (float) mLatticeData->GetXSiteCount(),
			  0.5F * (float) mLatticeData->GetYSiteCount(),
			  0.5F * (float) mLatticeData->GetZSiteCount());

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>
	  ("Found cluster: %f, %f, %f, %hu, %hu, %hu, %f, %f, %f, %f, %f, %f",
	   lNewCluster->minBlock.x,
	   lNewCluster->minBlock.y,
	   lNewCluster->minBlock.z,
	   lNewCluster->blocksX,
	   lNewCluster->blocksY,
	   lNewCluster->blocksZ,
	   lNewCluster->minSite.x,
	   lNewCluster->minSite.y,
	   lNewCluster->minSite.z,
	   lNewCluster->maxSite.x,
	   lNewCluster->maxSite.y,
	   lNewCluster->maxSite.z
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


      void ClusterBuilder::ProcessCluster(unsigned int iClusterId)
      {
	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
	  ("Examining cluster id = %u", (unsigned int) iClusterId);


	Cluster* lCluster = mClusters[iClusterId];
	
	lCluster->ResizeVectors();

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

	      
	      Vector3D<site_t>block_coordinates = Vector3D<site_t>(i, j, k) + mClusterBlockMins[iClusterId];
	      site_t lBlockId = mLatticeData->GetBlockIdFromBlockCoords
		(block_coordinates.x,
		 block_coordinates.y,
		 block_coordinates.z);
	      
	      Vector3D<site_t>* mins = &mClusterBlockMins[iClusterId]; 
	      assert(lBlockId == 
		     ( (i + mins->x) * mLatticeData->GetYBlockCount() + (j + mins->y))
		     * mLatticeData->GetZBlockCount() + (k + mins->z) 
		);

	      
	      if (mClusterIdOfBlock[lBlockId] != (short int) iClusterId)
	      {
		continue;
	      }
	      
	      ResizeVectorsForBlock(*lCluster, lBlockNum);

	      mClusterVoxelDataPointers.resize(mLatticeData->GetLocalFluidSiteCount());

	      UpdateSiteData(lBlockId, lBlockNum, iClusterId, block_coordinates);
	    } // for k
	  } // for j
	} // for i
      }
      
      void ClusterBuilder::ResizeVectorsForBlock(Cluster& iCluster, site_t lBlockNum)
      {
	//By default all values are -1 for solids
	iCluster.SiteData[lBlockNum].resize(
	  mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS, SiteData_t(-1.0F));
      }

      void ClusterBuilder::UpdateSiteData
      (site_t iBlockId, site_t iBlockNum,  unsigned int iClusterId,
       Vector3D<site_t>i_block_coordinates)
      {
	unsigned int lSiteId = -1;

	//Location site_coordinates_of_block = i_block_coordinates * mLatticeData->GetBlockSize();
	Vector3D<site_t>siteLocOnBlock;

	for (siteLocOnBlock.x = 0; siteLocOnBlock.x < mLatticeData->GetBlockSize(); siteLocOnBlock.x++)
	{
	  for (siteLocOnBlock.y = 0; siteLocOnBlock.y < mLatticeData->GetBlockSize(); siteLocOnBlock.y++)
	  {
	    for (siteLocOnBlock.z = 0; siteLocOnBlock.z < mLatticeData->GetBlockSize(); siteLocOnBlock.z++)
	    {
	      ++lSiteId;

	      UpdateSiteDataAtSite(iBlockId, iBlockNum, iClusterId, lSiteId);

	    }
	  }
	}
      } 
      

	     
      void ClusterBuilder::UpdateSiteDataAtSite
      (site_t iBlockId, site_t iBlockNum, 
       unsigned int iClusterId, unsigned int iSiteIdOnBlock)
      {
	geometry::LatticeData::BlockData * lBlock = mLatticeData->GetBlock(iBlockId);
	unsigned int lClusterVoxelSiteId = lBlock->site_data[iSiteIdOnBlock];

      	//If site not a solid and on the current processor [net.cc]
	if (lClusterVoxelSiteId != BIG_NUMBER3)
	{
	  UpdateDensityVelocityAndStress(iBlockNum, iClusterId, iSiteIdOnBlock, lClusterVoxelSiteId);
	}
      }

      void ClusterBuilder::UpdateDensityVelocityAndStress(site_t iBlockNum, unsigned int iClusterId, unsigned int iSiteIdOnBlock, unsigned int iClusterVoxelSiteId)
      {
	  mClusters[iClusterId]->SiteData[iBlockNum][iSiteIdOnBlock] = SiteData_t(1.0F);
	
	  //For efficiency we want to store a pointer to the site data grouped by the ClusterVortexID
	  //(1D organisation of sites)
	  std::vector<SiteData_t>::iterator lSiteDataIterator = 
	    mClusters[iClusterId]->SiteData[iBlockNum].begin() + iSiteIdOnBlock;

	  SiteData_t* lSiteDataLocation = &(*lSiteDataIterator); 
	  
	  //This one's for the C programmers out there
	  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
	    ("Siteid: %u", (unsigned int) iClusterVoxelSiteId);
	  SetDataPointerForClusterVoxelSiteId(iClusterVoxelSiteId, &(*lSiteDataLocation));
	} 
    

      SiteData_t* ClusterBuilder::GetDataPointerClusterVoxelSiteId(site_t iClusterVortexSiteId)
      {
	return mClusterVoxelDataPointers[iClusterVortexSiteId];
      }

      void ClusterBuilder::SetDataPointerForClusterVoxelSiteId
      (site_t iClusterVortexSiteId, SiteData_t* iDataPointer)
      {
#ifndef NDEBUG
	mClusterVoxelDataPointers.at(iClusterVortexSiteId) = iDataPointer;
#else 
	mClusterVoxelDataPointers[iClusterVortexSiteId] = iDataPointer;
#endif
      }

    }
  }
}
