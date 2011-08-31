#ifndef HEMELB_VIS_RAYTRACER_CLUSTERBUILDERWITHWALLNORMALS_H
#define HEMELB_VIS_RAYTRACER_CLUSTERBUILDERWITHWALLNORMALS_H

#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterBuilder.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      class ClusterBuilderWithWallNormals : public ClusterBuilder
      {
      public:
	ClusterBuilderWithWallNormals
	  (const geometry::LatticeData*& iLatticeData);

      private:
	virtual Cluster* CreateNewCluster();

	virtual void ResizeVectorsForBlock
	  (Cluster& iCluster,
	   site_t iBlockNum);
	
	virtual void UpdateSiteDataAtSite
	  (site_t iBlockId, site_t iBlockNum, 
	   unsigned int iClusterId, unsigned int iSiteIdOnBlock);

	void UpdateWallNormalAtSite
	  (geometry::LatticeData::BlockData * iBlock,
	   site_t iBlockNum, unsigned int iClusterId,
	   unsigned int iSiteIdOnBlock, unsigned int iClusterVoxelSiteId);

      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERBUILDERWITHWALLNORMALS_H
