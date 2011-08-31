#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterWithWallNormals::ClusterWithWallNormals()
      { }

      void ClusterWithWallNormals::ResizeVectors()
      {
	ResizeSharedVectors();
	WallNormals.resize(blocksX*blocksY*blocksZ);
      }

      void ClusterWithWallNormals::ResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
      {
	ResizeVectorsForBlockShared(iBlockNumber,iSize);
	WallNormals[iBlockNumber].resize(iSize, NULL);
      }
      
      double const*  ClusterWithWallNormals::GetWallData
      (site_t iBlockNumber, site_t iSiteNumber) const
      {
	return WallNormals[iBlockNumber][iSiteNumber];
      }

      void ClusterWithWallNormals::SetWallData
      (site_t iBlockNumber, site_t iSiteNumber, double* iData)
      {
	WallNormals[iBlockNumber][iSiteNumber] = iData;
      }

      bool ClusterWithWallNormals::NeedsWallNormals()
      {
	return true;
      }
    }
  }
}
