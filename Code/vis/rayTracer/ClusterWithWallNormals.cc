#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterWithWallNormals::ClusterWithWallNormals(unsigned short xBlockCount,
                                                     unsigned short yBlockCount,
                                                     unsigned short zBlockCount,
                                                     const util::Vector3D<float>& minimalSite,
                                                     const util::Vector3D<float>& maximalSite,
                                                     const util::Vector3D<float>& minimalSiteOnMinimalBlock,
                                                     const util::Vector3D<site_t>& minimalBlock) :
          Cluster<ClusterWithWallNormals>(xBlockCount,
                                          yBlockCount,
                                          zBlockCount,
                                          minimalSite,
                                          maximalSite,
                                          minimalSiteOnMinimalBlock,
                                          minimalBlock)
      {
        WallNormals.resize(GetBlocksX() * GetBlocksY() * GetBlocksZ());
      }

      void ClusterWithWallNormals::DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
      {
        WallNormals[iBlockNumber].resize(iSize, NULL);
      }

      const util::Vector3D<double>* ClusterWithWallNormals::DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const
      {
        return WallNormals[iBlockNumber][iSiteNumber];
      }

      void ClusterWithWallNormals::DoSetWallData(site_t iBlockNumber,
                                                 site_t iSiteNumber,
                                                 const util::Vector3D<double>& iData)
      {
        WallNormals[iBlockNumber][iSiteNumber] = &iData;
      }

      bool ClusterWithWallNormals::DoNeedsWallNormals()
      {
        return true;
      }
    }
  }
}
