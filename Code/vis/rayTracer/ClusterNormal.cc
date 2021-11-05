// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterNormal.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterNormal::ClusterNormal(unsigned short xBlockCount, unsigned short yBlockCount,
                                   unsigned short zBlockCount,
                                   const util::Vector3D<float>& minimalSite,
                                   const util::Vector3D<float>& maximalSite,
                                   const util::Vector3D<float>& minimalSiteOnMinimalBlock,
                                   const util::Vector3D<site_t>& minimalBlock) :
              Cluster<ClusterNormal>(xBlockCount,
                                     yBlockCount,
                                     zBlockCount,
                                     minimalSite,
                                     maximalSite,
                                     minimalSiteOnMinimalBlock,
                                     minimalBlock)
      {
      }

      const util::Vector3D<double>* ClusterNormal::DoGetWallData(site_t iBlockNumber,
                                                                 site_t iSiteNumber) const
      {
        return nullptr;
      }

      void ClusterNormal::DoSetWallData(site_t iBlockNumber, site_t iSiteNumber,
                                        const util::Vector3D<double>& iData)
      {
      }

    }
  }
}
