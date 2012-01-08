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
                                 const util::Vector3D<float>& minimalSiteOnMinimalBlock);

          void DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize);

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
