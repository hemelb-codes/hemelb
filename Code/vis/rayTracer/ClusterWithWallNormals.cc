// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterWithWallNormals.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterWithWallNormals::ClusterWithWallNormals(
          unsigned short xBlockCount, unsigned short yBlockCount, unsigned short zBlockCount,
          const util::Vector3D<float>& minimalSite, const util::Vector3D<float>& maximalSite,
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

      const util::Vector3D<double>* ClusterWithWallNormals::DoGetWallData(site_t blockNumber,
                                                                          site_t siteNumber) const
      {
        if (siteNumber < (site_t) WallNormals[blockNumber].size())
        {
          return WallNormals[blockNumber][siteNumber];
        }
        else
        {
          return nullptr;
        }
      }

      void ClusterWithWallNormals::DoSetWallData(site_t blockNumber, site_t siteNumber,
                                                 const util::Vector3D<double>& data)
      {
        if (WallNormals[blockNumber].size() <= (size_t) siteNumber)
        {
          WallNormals[blockNumber].resize(siteNumber + 1, nullptr);
        }

        WallNormals[blockNumber][siteNumber] = &data;
      }

      bool ClusterWithWallNormals::DoNeedsWallNormals()
      {
        return true;
      }
    }
  }
}
