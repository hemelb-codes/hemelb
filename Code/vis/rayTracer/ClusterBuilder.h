#ifndef HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
#define HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H

//#define NDEBUG
#include <assert.h>
#include <map>
#include <vector>

#include <iostream>
 
#include "debug/Debugger.h"
#include "geometry/LatticeData.h"
#include "geometry/SiteTraverser.h"
#include "geometry/VolumeTraverser.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "vis/rayTracer/BlockTraverserWithVisitedBlockTracker.h"
#include "vis/rayTracer/ClusterBuilder.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/rayTracer/SiteData.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {				
      template <typename ClusterType>
      class ClusterBuilder
      {
      public:
	ClusterBuilder
	  (const geometry::LatticeData* iLatticeData) :
	  mBlockTraverser(*iLatticeData)
	{
	  mLatticeData = iLatticeData;

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

        ~ClusterBuilder()
	{
	  delete[] mClusterIdOfBlock;
	}
	  
	void BuildClusters()
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
	  
	std::vector<ClusterType>& GetClusters()
	{
	  return mClusters;
	}

	SiteData_t* GetClusterVoxelDataPointer(site_t iVoxelSiteId)
	{
	  return GetDataPointerClusterVoxelSiteId(iVoxelSiteId);
	}

      private:
	void ResizeVectorsForBlock(ClusterType& iCluster, site_t iBlockNum)
	{
	  iCluster.ResizeVectorsForBlock
	    (iBlockNum, 
	     mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS);
	}

      private:	//Locates all the clusters in the lattice structure and the
	void LocateClusters()
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
 
	//Locates all the clusters in the lattice structure and the
	void FindNewCluster()
	{
	  //These locations will eventually contain the bounds of the
	  //rectangular cluster, both in terms of block number and
	  //site numbers
	  util::Vector3D<site_t> lClusterBlockMin = util::Vector3D<site_t>::MaxLimit();
	  util::Vector3D<site_t> lClusterBlockMax = util::Vector3D<site_t>::MinLimit();
	  util::Vector3D<site_t> lClusterSiteMin = util::Vector3D<site_t>::MaxLimit();
	  util::Vector3D<site_t> lClusterSiteMax = util::Vector3D<site_t>::MinLimit();

	  //To discover the cluster, we continually visit the neighbours 
	  //of sequential blocks 
	  //We keep a stack of all the sites that must be processed 
	  //and sequentially add neighbours to it
	  std::stack<util::Vector3D<site_t> > lBlocksToProcess;

	  //Set up the initial condition
	  lBlocksToProcess.push(mBlockTraverser.GetCurrentLocation());

	  //Loop over the cluster via neighbours until
	  //all blocks have been processed
	  while (!lBlocksToProcess.empty())
	  {
	    //Get location off the top of the stack
	    //(we could actually take anything off the stack)
	    util::Vector3D<site_t> lCurrentLocation = lBlocksToProcess.top();
	    lBlocksToProcess.pop();
	    
	    if(AreSitesAssignedToLocalProcessorRankInBlock(
		 mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)))
	    {
	      //Update block range of the cluster
	      lClusterBlockMin.UpdatePointwiseMin(lCurrentLocation);
	      lClusterBlockMax.UpdatePointwiseMax(lCurrentLocation);
	      
	      //Update the cluster id of the given block
	      site_t lBlockID = mBlockTraverser.GetIndexFromLocation(lCurrentLocation);
	      mClusterIdOfBlock[lBlockID] = (short int) mClusters.size();

	      //Loop through all the sites on the block, to 
	      //update the site bounds on the cluster
	      geometry::SiteTraverser lSiteTraverser = mBlockTraverser.GetSiteTraverser();
	      do
	      { 
		//If the site is not a solid
		if (mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)->
		    site_data[lSiteTraverser.GetCurrentIndex()] != BIG_NUMBER3)
		{
		  lClusterSiteMin.UpdatePointwiseMin
		    ( lSiteTraverser.GetCurrentLocation() +
		      lCurrentLocation*mBlockTraverser.GetBlockSize());
		
		  lClusterSiteMax.UpdatePointwiseMax
		    ( lSiteTraverser.GetCurrentLocation() +
		      lCurrentLocation*mBlockTraverser.GetBlockSize());
		}   
	      }
	      while (lSiteTraverser.TraverseOne());
			
	      //Check all the neighbouring blocks to see if they need visiting. Add them to the stack.
	      AddNeighbouringBlocks(lCurrentLocation, lBlocksToProcess);
	    }
	  }
	
	  AddCluster(lClusterBlockMin, lClusterBlockMax, lClusterSiteMin, lClusterSiteMax);
	}

	//Adds neighbouring blocks of the input location to the input stack
	void AddNeighbouringBlocks(util::Vector3D<site_t> iCurrentLocation,
				   std::stack<util::Vector3D<site_t> >& oBlocksToProcess)
	{
	  // Loop over all neighbouring blocks
	  for (int l = 0; l < 26; l++)
	  {
	    util::Vector3D<site_t> lNeighbouringBlock = iCurrentLocation + mNeighbours[l];
		
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

	//Returns true if there are sites in the given block associated with the
	//local processor rank
	bool AreSitesAssignedToLocalProcessorRankInBlock
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

	//Adds a new cluster by taking in the required data in interget format
	//and converting it to that used by the raytracer
	//NB: Futher processing is required on the cluster before it can be used
	//by the ray tracer, which is handled by the ProcessCluster method
	void AddCluster(util::Vector3D<site_t> iClusterBlockMin,
			util::Vector3D<site_t> iClusterBlockMax,
			util::Vector3D<site_t> iClusterVoxelMin,
			util::Vector3D<site_t> iClusterVoxelMax)
	{
//The friendly locations must be turned into a format usable by the ray tracer
	  ClusterType lNewCluster;
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
	
	  lNewCluster.minSite = util::Vector3D<float>(iClusterVoxelMin) -
	    util::Vector3D<float>((float) mLatticeData->GetXSiteCount(),
			    (float) mLatticeData->GetYSiteCount(),
			    (float) mLatticeData->GetZSiteCount()) * 0.5F;
	
	  lNewCluster.maxSite = util::Vector3D<float>(iClusterVoxelMax + util::Vector3D<site_t>(1)) - 
	    util::Vector3D<float>(0.5F * (float) mLatticeData->GetXSiteCount(),
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

	//Adds "flow-field" data to the cluster
	void ProcessCluster(unsigned int iClusterId)
	{
	  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
	    ("Examining cluster id = %u", (unsigned int) iClusterId);


	  ClusterType& lCluster(mClusters[iClusterId]);
	
	  lCluster.ResizeVectors();

	  site_t lBlockNum = -1;
	  for (site_t i = 0; i < lCluster.blocksX; i++)
	  {
	    for (site_t j = 0; j < lCluster.blocksY; j++)
	    {
	      for (site_t k = 0; k < lCluster.blocksZ; k++)
	      {
		++lBlockNum;

		hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
		  ("Examining block number = %u", (unsigned int) lBlockNum);

	      
		util::Vector3D<site_t>block_coordinates = util::Vector3D<site_t>(i, j, k) + mClusterBlockMins[iClusterId];
		site_t lBlockId = mLatticeData->GetBlockIdFromBlockCoords
		  (block_coordinates.x,
		   block_coordinates.y,
		   block_coordinates.z);
	      
		util::Vector3D<site_t>* mins = &mClusterBlockMins[iClusterId]; 
		assert(lBlockId == 
		       ( (i + mins->x) * mLatticeData->GetYBlockCount() + (j + mins->y))
		       * mLatticeData->GetZBlockCount() + (k + mins->z) 
		  );

	      
		if (mClusterIdOfBlock[lBlockId] != (short int) iClusterId)
		{
		  continue;
		}
	      
		ResizeVectorsForBlock(lCluster, lBlockNum);

		mClusterVoxelDataPointers.resize(mLatticeData->GetLocalFluidSiteCount());

		UpdateSiteData(lBlockId, lBlockNum, lCluster, block_coordinates);
	      } // for k
	    } // for j
	  } // for i
	}

	void UpdateSiteData
	  (site_t iBlockId, site_t iBlockNum,  ClusterType& iCluster,
	   util::Vector3D<site_t>i_block_coordinates)
	{
	  unsigned int lSiteId = -1;

	  //Location site_coordinates_of_block = i_block_coordinates * mLatticeData->GetBlockSize();
	  util::Vector3D<site_t>siteLocOnBlock;

	  for (siteLocOnBlock.x = 0; siteLocOnBlock.x < mLatticeData->GetBlockSize(); siteLocOnBlock.x++)
	  {
	    for (siteLocOnBlock.y = 0; siteLocOnBlock.y < mLatticeData->GetBlockSize(); siteLocOnBlock.y++)
	    {
	      for (siteLocOnBlock.z = 0; siteLocOnBlock.z < mLatticeData->GetBlockSize(); siteLocOnBlock.z++)
	      {
		++lSiteId;

		UpdateSiteDataAtSite(iBlockId, iBlockNum, iCluster, lSiteId);

	      }
	    }
	  }
	}

	virtual void UpdateSiteDataAtSite
	  (site_t iBlockId, site_t iBlockNum, 
	   ClusterType& iCluster, unsigned int iSiteIdOnBlock)
	{
	  geometry::LatticeData::BlockData * lBlock = mLatticeData->GetBlock(iBlockId);
	  unsigned int lClusterVoxelSiteId = lBlock->site_data[iSiteIdOnBlock];

	  //If site not a solid and on the current processor [net.cc]
	  if (lClusterVoxelSiteId != BIG_NUMBER3)
	  {
	    UpdateDensityVelocityAndStress(iBlockNum, iCluster, iSiteIdOnBlock, lClusterVoxelSiteId);
	    if (ClusterType::NeedsWallNormals())
	    {
	      UpdateWallNormalAtSite(lBlock, iBlockNum, iCluster, iSiteIdOnBlock);
	    }
	  }
	}

	void UpdateDensityVelocityAndStress(site_t iBlockNum, ClusterType& iCluster, unsigned int iSiteIdOnBlock, unsigned int iClusterVoxelSiteId)
	{
	  iCluster.SiteData[iBlockNum][iSiteIdOnBlock] = SiteData_t(1.0F);
	
	  //For efficiency we want to store a pointer to the site data grouped by the ClusterVortexID
	  //(1D organisation of sites)
	  std::vector<SiteData_t>::iterator lSiteDataIterator = 
	    iCluster.SiteData[iBlockNum].begin() + iSiteIdOnBlock;

	  SiteData_t* lSiteDataLocation = &(*lSiteDataIterator); 
	  
	  //This one's for the C programmers out there
	  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>
	    ("Siteid: %u", (unsigned int) iClusterVoxelSiteId);
	  SetDataPointerForClusterVoxelSiteId(iClusterVoxelSiteId, &(*lSiteDataLocation));
	}
	
	void UpdateWallNormalAtSite(
	  geometry::LatticeData::BlockData * iBlock,
	  site_t iBlockNum, ClusterType& iCluster,
	  unsigned int iSiteIdOnBlock)
	{
	  
	  if(iBlock->wall_data[iSiteIdOnBlock].wall_nor[0] != -1.0F)
	  {
	    iCluster.SetWallData
	      (iBlockNum,
	       iSiteIdOnBlock,
	       iBlock->wall_data[iSiteIdOnBlock].wall_nor);
	  } 
	}
	
	
	util::Vector3D<site_t> GetSiteCoordinatesOfBlock
	  (site_t iClusterId, util::Vector3D<site_t> offset);

	SiteData_t* GetDataPointerClusterVoxelSiteId(site_t iClusterVortexSiteId)
	{
	  return mClusterVoxelDataPointers[iClusterVortexSiteId];
	}
	
	void SetDataPointerForClusterVoxelSiteId
	  (site_t iClusterVortexSiteId,
	   SiteData_t* iDataPointer)
	{
#ifndef NDEBUG
	  mClusterVoxelDataPointers.at(iClusterVortexSiteId) = iDataPointer;
#else
	  mClusterVoxelDataPointers[iClusterVortexSiteId] = iDataPointer;
#endif 
	}
	
	//Caution: the data within mClusters is altered by means
	//of pointers obtained from the GetClusterVoxelDataPointer
	//method. No insertion of copying must therefore take place
	//on mClusters once building is complete
	std::vector<ClusterType> mClusters;

	const geometry::LatticeData* mLatticeData;


	BlockTraverserWithVisitedBlockTracker
	  mBlockTraverser;

	std::vector<util::Vector3D<site_t> > mClusterBlockMins;

	//This allows a cluster voxel site ID (as part of the 1D structure for)
	//storing sites to be mapped to the data stored in the 3D structure
	//for the ray tracer by means of pointers.
	std::vector<SiteData_t*> mClusterVoxelDataPointers;

	short int *mClusterIdOfBlock;

	static const short int NOTASSIGNEDTOCLUSTER = -1;

	static const util::Vector3D<site_t> mNeighbours[26];
      };
      
      template <typename ClusterType> 
	const util::Vector3D<site_t> ClusterBuilder<ClusterType>::mNeighbours[26] =
	{
	  util::Vector3D<site_t>(-1, -1, -1),
	  util::Vector3D<site_t>(-1, -1, 0),
	  util::Vector3D<site_t>(-1, -1, 1),
	  util::Vector3D<site_t>(-1, 0, -1),
	  util::Vector3D<site_t>(-1, 0, 0),
	  util::Vector3D<site_t>(-1, 0, 1),
	  util::Vector3D<site_t>(-1, 1, -1),
	  util::Vector3D<site_t>(-1, 1, 0),
	  util::Vector3D<site_t>(-1, 1, 1),
	  util::Vector3D<site_t>( 0, -1, -1),
	  util::Vector3D<site_t>( 0, -1, 0),
	  util::Vector3D<site_t>( 0, -1, 1),
	  util::Vector3D<site_t>( 0, 0, -1),
	  // 0 0 0 is same site
	  util::Vector3D<site_t>( 0, 0, 1),
	  util::Vector3D<site_t>( 0, 1, -1),
	  util::Vector3D<site_t>( 0, 1, 0),
	  util::Vector3D<site_t>( 0, 1, 1),
	  util::Vector3D<site_t>( 1, -1, -1),
	  util::Vector3D<site_t>( 1, -1, 0),
	  util::Vector3D<site_t>( 1, -1, 1),
	  util::Vector3D<site_t>( 1, 0, -1),
	  util::Vector3D<site_t>( 1, 0, 0),
	  util::Vector3D<site_t>( 1, 0, 1),
	  util::Vector3D<site_t>( 1, 1, -1),
	  util::Vector3D<site_t>( 1, 1, 0),
	  util::Vector3D<site_t>( 1, 1, 1)
	}; 
      
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
