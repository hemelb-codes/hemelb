#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterWithWallNormals::ClusterWithWallNormals()
      {
      }

      void ClusterWithWallNormals::DoResizeVectors()
      {
        ResizeSharedVectors();
        WallNormals.resize(blocksX * blocksY * blocksZ);
      }

      void ClusterWithWallNormals::DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
      {
        DoResizeVectorsForBlockShared(iBlockNumber, iSize);
        WallNormals[iBlockNumber].resize(iSize, NULL);
      }

      double const* ClusterWithWallNormals::DoGetWallData(site_t iBlockNumber,
                                                          site_t iSiteNumber) const
      {
        return WallNormals[iBlockNumber][iSiteNumber];
      }

      void ClusterWithWallNormals::DoSetWallData(site_t iBlockNumber,
                                                 site_t iSiteNumber,
                                                 double* iData)
      {
        WallNormals[iBlockNumber][iSiteNumber] = iData;
      }

      bool ClusterWithWallNormals::DoNeedsWallNormals()
      {
        return true;
      }
    }
  }
}
