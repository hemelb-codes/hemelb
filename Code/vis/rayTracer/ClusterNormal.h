#ifndef HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H
#define HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H

#include "vis/rayTracer/ClusterShared.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      class ClusterNormal : public ClusterShared<ClusterNormal>
      {
        public:
          ClusterNormal();

          void DoResizeVectors();

          void DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize);

          double const* DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const;

          void DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData);

        private:
          std::vector<std::vector<double*> > WallNormals;
      };

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H
