#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      RayTracer::ClusterBuilder::SiteTraverser::SiteTraverser(const geometry::LatticeData* iLatticeData,
                                                              const site_t iBlockId) :
          mLatticeData(iLatticeData), mBlockId(iBlockId)
      {
      }

      site_t RayTracer::ClusterBuilder::SiteTraverser::GetXCount()
      {
        return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteTraverser::GetYCount()
      {
        return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteTraverser::GetZCount()
      {
        return GetBlockSize();
      }

      site_t RayTracer::ClusterBuilder::SiteTraverser::GetBlockSize()
      {
        return mLatticeData->GetBlockSize();
      }
    }
  }
}
