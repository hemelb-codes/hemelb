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

      void ClusterNormal::DoResizeVectors()
      {
	ResizeSharedVectors();
      }

      void ClusterNormal::DoResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
      {
	DoResizeVectorsForBlockShared(iBlockNumber,iSize);
      }
      
      double const*  ClusterNormal::DoGetWallData
      (site_t iBlockNumber, site_t iSiteNumber) const
      {
	return NULL;
      }

      void ClusterNormal::DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData)
      {
	assert(false);
      }

    }
  }
}
