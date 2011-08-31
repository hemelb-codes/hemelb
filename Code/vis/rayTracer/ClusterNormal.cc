//#define NDEBUG
#include <cassert>

#include "geometry/LatticeData.h"
#include "vis/rayTracer/ClusterNormal.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      ClusterNormal::ClusterNormal()
      { }

      void ClusterNormal::ResizeVectors()
      {
	ResizeSharedVectors();
      }

      void ResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
      {
	ResizeVectorsForBlockShared(iBlockNumber,iSize);
      }
      
      double const*  ClusterNormal::GetWallData
      (site_t iBlockNumber, site_t iSiteNumber) const
      {
	return NULL;
      }

      void ClusterNormal::SetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData)
      {
	assert(false);
      }

    }
  }
}
