#include <iostream>

#include "vis/rayTracer/ClusterWithWallNormals.h"
#include "vis/rayTracer/ClusterBuilderWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterBuilderWithWallNormals::ClusterBuilderWithWallNormals
      (const geometry::LatticeData*& iLatticeData) : 
	ClusterBuilder(iLatticeData)
      {
      }

      Cluster* ClusterBuilderWithWallNormals::CreateNewCluster()
      {
	return new ClusterWithWallNormals();
      }    

      void ClusterBuilderWithWallNormals::ResizeVectorsForBlock(Cluster& iCluster, site_t lBlockNum)
      {
	ClusterBuilder::ResizeVectorsForBlock(iCluster, lBlockNum);

	ClusterWithWallNormals& lCluster = dynamic_cast
	  <ClusterWithWallNormals&>(iCluster);
	
	lCluster.WallNormals.at(lBlockNum).resize(
	  mLatticeData->GetSitesPerBlockVolumeUnit() * VIS_FIELDS, NULL);
      }


      void ClusterBuilderWithWallNormals::UpdateSiteDataAtSite
      ( site_t iBlockId, site_t iBlockNum, 
	unsigned int iClusterId, unsigned int iSiteIdOnBlock)
      {
	geometry::LatticeData::BlockData * lBlock = mLatticeData->GetBlock(iBlockId);
	unsigned int lClusterVoxelSiteId = lBlock->site_data[iSiteIdOnBlock];

      	//If site not a solid and on the current processor [net.cc]
	if (lClusterVoxelSiteId != BIG_NUMBER3)
	{
	  UpdateDensityVelocityAndStress(iBlockNum, iClusterId, iSiteIdOnBlock, lClusterVoxelSiteId);

	  UpdateWallNormalAtSite(lBlock, iBlockNum, iClusterId, iSiteIdOnBlock, lClusterVoxelSiteId);
	}
      }

      void ClusterBuilderWithWallNormals::UpdateWallNormalAtSite(
	geometry::LatticeData::BlockData * iBlock,
	site_t iBlockNum, unsigned int iClusterId,
	unsigned int iSiteIdOnBlock, unsigned int iClusterVoxelSiteId)
      {
	ClusterWithWallNormals* lCluster = 
	  dynamic_cast<ClusterWithWallNormals*>
	  (mClusters[iClusterId]);
	if(iBlock->wall_data[iSiteIdOnBlock].wall_nor[0] != -1.0F)
	{
	  lCluster->WallNormals[iBlockNum][iSiteIdOnBlock] = iBlock->wall_data[iSiteIdOnBlock].wall_nor;
	}
      }
    }
  }
}
