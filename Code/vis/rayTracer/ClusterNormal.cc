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
                                   const util::Vector3D<float>& minimalSiteOnMinimalBlock) :
            Cluster<ClusterNormal> (xBlockCount,
                                    yBlockCount,
                                    zBlockCount,
                                    minimalSite,
                                    maximalSite,
                                    minimalSiteOnMinimalBlock)
      {
      }

      const util::Vector3D<double>* ClusterNormal::DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const
      {
        return NULL;
      }

      void ClusterNormal::DoSetWallData(site_t iBlockNumber,
                                        site_t iSiteNumber,
                                        const util::Vector3D<double>& iData)
      {
      }

    }
  }
}
