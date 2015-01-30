// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterNormal.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterNormal::ClusterNormal(unsigned short xBlockCount,
                                   unsigned short yBlockCount,
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

      const util::Vector3D<double>* ClusterNormal::DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const
      {
        return nullptr;
      }

      void ClusterNormal::DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, const util::Vector3D<double>& iData)
      {
      }

    }
  }
}
