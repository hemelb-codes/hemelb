// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
#define HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      class ClusterWithWallNormals : public Cluster<ClusterWithWallNormals>
      {
        public:
          ClusterWithWallNormals(unsigned short xBlockCount,
                                 unsigned short yBlockCount,
                                 unsigned short zBlockCount,
                                 const util::Vector3D<float>& minimalSite,
                                 const util::Vector3D<float>& maximalSite,
                                 const util::Vector3D<float>& minimalSiteOnMinimalBlock,
                                 const util::Vector3D<site_t>& minimalBlock);

          const util::Vector3D<double>* DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const;

          void DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, const util::Vector3D<double>& iData);

          static bool DoNeedsWallNormals();

        private:
          std::vector<std::vector<const util::Vector3D<double>*> > WallNormals;
      };

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERWITHWALLNORMALS_H
