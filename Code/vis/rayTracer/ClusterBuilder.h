#ifndef HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
#define HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H

#include <map>
#include <vector>		

#include "geometry/BlockTraverserWithVisitedBlockTracker.h"
#include "geometry/LatticeData.h"
#include "geometry/SiteTraverser.h"
#include "lb/LbmParameters.h"
#include "util/utilityFunctions.h"
#include "util/Vector3D.h"
#include "vis/rayTracer/ClusterBuilder.h"
#include "vis/rayTracer/ClusterTraverser.h"
#include "vis/rayTracer/RayTracer.h"
#include "vis/rayTracer/SiteData.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      template<typename ClusterType>
      class ClusterBuilder
      {
        public:
          ClusterBuilder(const geometry::LatticeData* iLatticeData) :
            mBlockTraverser(*iLatticeData)
          {
            mLatticeData = iLatticeData;

            //Each block is assigned a cluster id once it has been
            //assigned to a cluster
            mClusterIdOfBlock = new short int[mLatticeData->GetBlockCount()];
            for (site_t lId = 0; lId < mLatticeData->GetBlockCount(); lId++)
            {
              mClusterIdOfBlock[lId] = NOTASSIGNEDTOCLUSTER;
            }
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
            for (unsigned int lThisClusterId = 0; lThisClusterId < mClusters.size(); lThisClusterId++)
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
            iCluster.ResizeVectorsForBlock(iBlockNum,
                                           mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS);
          }

          // Locates all the clusters in the lattice structure and the
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

          //Locates all the clusters in the lattice structure and stores
          //their locations
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

              if (AreSitesAssignedToLocalProcessorRankInBlock(mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)))
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
                  if (mBlockTraverser.GetBlockDataForLocation(lCurrentLocation)->site_data[lSiteTraverser.GetCurrentIndex()]
                      != BIG_NUMBER3)
                  {
                    lClusterSiteMin.UpdatePointwiseMin(lSiteTraverser.GetCurrentLocation()
                        + lCurrentLocation * mBlockTraverser.GetBlockSize());

                    lClusterSiteMax.UpdatePointwiseMax(lSiteTraverser.GetCurrentLocation()
                        + lCurrentLocation * mBlockTraverser.GetBlockSize());
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
                if (!mBlockTraverser.IsBlockVisited(lNeighbouringBlock))
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
          bool AreSitesAssignedToLocalProcessorRankInBlock(geometry::BlockData * iBlock)
          {
            if (iBlock->ProcessorRankForEachBlockSite == NULL)
            {
              return false;
            }

            for (unsigned int siteId = 0; siteId < mLatticeData->GetSitesPerBlockVolumeUnit(); siteId++)
            {
              if (topology::NetworkTopology::Instance()->GetLocalRank()
                  == iBlock->ProcessorRankForEachBlockSite[siteId])
              {
                return true;
              }
            }
            return false;
          }

          //Adds a new cluster by taking in the required data in interger format
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
            lNewCluster.minBlock.x = (float) (iClusterBlockMin.x * mLatticeData->GetBlockSize())
                - 0.5F * (float) mLatticeData->GetXSiteCount();
            lNewCluster.minBlock.y = (float) (iClusterBlockMin.y * mLatticeData->GetBlockSize())
                - 0.5F * (float) mLatticeData->GetYSiteCount();
            lNewCluster.minBlock.z = (float) (iClusterBlockMin.z * mLatticeData->GetBlockSize())
                - 0.5F * (float) mLatticeData->GetZSiteCount();

            lNewCluster.blocksX = (unsigned short) (1 + iClusterBlockMax.x - iClusterBlockMin.x);
            lNewCluster .blocksY = (unsigned short) (1 + iClusterBlockMax.y - iClusterBlockMin.y);
            lNewCluster .blocksZ = (unsigned short) (1 + iClusterBlockMax.z - iClusterBlockMin.z);

            lNewCluster.minSite = util::Vector3D<float>(iClusterVoxelMin)
                - util::Vector3D<float>((float) mLatticeData->GetXSiteCount(),
                                        (float) mLatticeData->GetYSiteCount(),
                                        (float) mLatticeData->GetZSiteCount()) * 0.5F;

            lNewCluster.maxSite = util::Vector3D<float>(iClusterVoxelMax
                + util::Vector3D<site_t>(1))
                - util::Vector3D<float>(0.5F * (float) mLatticeData->GetXSiteCount(),
                                        0.5F * (float) mLatticeData->GetYSiteCount(),
                                        0.5F * (float) mLatticeData->GetZSiteCount());

            mClusters.push_back(lNewCluster);

            //We need to store the cluster block minimum in
            //order to process the cluster
            mClusterBlockMins.push_back(iClusterBlockMin);
          }

          //Adds "flow-field" data to the cluster
          void ProcessCluster(unsigned int iClusterId)
          {
            hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("Examining cluster id = %u",
                                                                                  (unsigned int) iClusterId);

            ClusterType& lCluster(mClusters[iClusterId]);

            lCluster.ResizeVectors();

            ClusterTraverser<ClusterType> lClusterTraverser(lCluster);

            do
            {
              util::Vector3D<site_t> lBlockCoordinates = lClusterTraverser.GetCurrentLocation()
                  + mClusterBlockMins[iClusterId];

              site_t lBlockId = mLatticeData->GetBlockIdFromBlockCoords(lBlockCoordinates.x,
                                                                        lBlockCoordinates.y,
                                                                        lBlockCoordinates.z);

              if (mClusterIdOfBlock[lBlockId] == iClusterId)
              {
                ResizeVectorsForBlock(lCluster, lClusterTraverser.GetCurrentIndex());

                mClusterVoxelDataPointers.resize(mLatticeData->GetLocalFluidSiteCount());

                UpdateSiteData(lBlockId, lClusterTraverser.GetCurrentIndex(), lCluster);
              }
            }
            while (lClusterTraverser.TraverseOne());
          }

          void UpdateSiteData(site_t iBlockId, site_t iBlockNum, ClusterType& iCluster)
          {
            geometry::SiteTraverser lSiteTraverser(*mLatticeData);
            do
            {
              UpdateSiteDataAtSite(iBlockId, iBlockNum, iCluster, lSiteTraverser.GetCurrentIndex());
            }
            while (lSiteTraverser.TraverseOne());
          }

          virtual void UpdateSiteDataAtSite(site_t iBlockId,
                                            site_t iBlockNum,
                                            ClusterType& iCluster,
                                            unsigned int iSiteIdOnBlock)
          {
            geometry::BlockData * lBlock = mLatticeData->GetBlock(iBlockId);
            unsigned int lClusterVoxelSiteId = lBlock->site_data[iSiteIdOnBlock];

            //If site not a solid and on the current processor [net.cc]
            if (lClusterVoxelSiteId != BIG_NUMBER3)
            {
              UpdateDensityVelocityAndStress(iBlockNum,
                                             iCluster,
                                             iSiteIdOnBlock,
                                             lClusterVoxelSiteId);
              if (ClusterType::NeedsWallNormals())
              {
                UpdateWallNormalAtSite(lBlock, iBlockNum, iCluster, iSiteIdOnBlock);
              }
            }
          }

          void UpdateDensityVelocityAndStress(site_t iBlockNum,
                                              ClusterType& iCluster,
                                              unsigned int iSiteIdOnBlock,
                                              unsigned int iClusterVoxelSiteId)
          {
            iCluster.SiteData[iBlockNum][iSiteIdOnBlock] = SiteData_t(1.0F);

            //For efficiency we want to store a pointer to the site data grouped by the ClusterVortexID
            //(1D organisation of sites)
            std::vector<SiteData_t>::iterator lSiteDataIterator =
                iCluster.SiteData[iBlockNum].begin() + iSiteIdOnBlock;

            SiteData_t* lSiteDataLocation = & (*lSiteDataIterator);

            SetDataPointerForClusterVoxelSiteId(iClusterVoxelSiteId, & (*lSiteDataLocation));
          }

          void UpdateWallNormalAtSite(geometry::BlockData * iBlock,
                                      site_t iBlockNum,
                                      ClusterType& iCluster,
                                      unsigned int iSiteIdOnBlock)
          {

            if (iBlock->wall_data[iSiteIdOnBlock].wall_nor[0] != -1.0F)
            {
              iCluster.SetWallData(iBlockNum,
                                   iSiteIdOnBlock,
                                   iBlock->wall_data[iSiteIdOnBlock].wall_nor);
            }
          }

          SiteData_t* GetDataPointerClusterVoxelSiteId(site_t iClusterVortexSiteId)
          {
            return mClusterVoxelDataPointers[iClusterVortexSiteId];
          }

          void SetDataPointerForClusterVoxelSiteId(site_t iClusterVortexSiteId,
                                                   SiteData_t* iDataPointer)
          {
            mClusterVoxelDataPointers[iClusterVortexSiteId] = iDataPointer;
          }

          //Caution: the data within mClusters is altered by means
          //of pointers obtained from the GetClusterVoxelDataPointer
          //method. No insertion or copying must therefore take place
          //on mClusters once building is complete
          std::vector<ClusterType> mClusters;

          const geometry::LatticeData* mLatticeData;

          geometry::BlockTraverserWithVisitedBlockTracker mBlockTraverser;

          std::vector<util::Vector3D<site_t> > mClusterBlockMins;

          //This allows a cluster voxel site ID (as part of the 1D structure for)
          //storing sites to be mapped to the data stored in the 3D structure
          //for the ray tracer by means of pointers.
          std::vector<SiteData_t*> mClusterVoxelDataPointers;

          short int *mClusterIdOfBlock;

          static const short int NOTASSIGNEDTOCLUSTER = -1;

          static const util::Vector3D<site_t> mNeighbours[26];
      };

      template<typename ClusterType>
      const util::Vector3D<site_t> ClusterBuilder<ClusterType>::mNeighbours[26] =
          { util::Vector3D<site_t>(-1, -1, -1),
            util::Vector3D<site_t>(-1, -1, 0),
            util::Vector3D<site_t>(-1, -1, 1),
            util::Vector3D<site_t>(-1, 0, -1),
            util::Vector3D<site_t>(-1, 0, 0),
            util::Vector3D<site_t>(-1, 0, 1),
            util::Vector3D<site_t>(-1, 1, -1),
            util::Vector3D<site_t>(-1, 1, 0),
            util::Vector3D<site_t>(-1, 1, 1),
            util::Vector3D<site_t>(0, -1, -1),
            util::Vector3D<site_t>(0, -1, 0),
            util::Vector3D<site_t>(0, -1, 1),
            util::Vector3D<site_t>(0, 0, -1),
            // 0 0 0 is same site
            util::Vector3D<site_t>(0, 0, 1),
            util::Vector3D<site_t>(0, 1, -1),
            util::Vector3D<site_t>(0, 1, 0),
            util::Vector3D<site_t>(0, 1, 1),
            util::Vector3D<site_t>(1, -1, -1),
            util::Vector3D<site_t>(1, -1, 0),
            util::Vector3D<site_t>(1, -1, 1),
            util::Vector3D<site_t>(1, 0, -1),
            util::Vector3D<site_t>(1, 0, 0),
            util::Vector3D<site_t>(1, 0, 1),
            util::Vector3D<site_t>(1, 1, -1),
            util::Vector3D<site_t>(1, 1, 0),
            util::Vector3D<site_t>(1, 1, 1) };

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERBUILDER_H
