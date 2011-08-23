//#define NDEBUG;
#include <assert.h>

#include "vis/rayTracer/ClusterBuilderWithWallNormals.h"
#include "vis/rayTracer/RayTracerWithLighting.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      void RayTracerWithLighting::BuildClusters()
      {
	mClusterBuilder = new ClusterBuilderWithWallNormals(mLatDat);
	mClusterBuilder->BuildClusters();
      }
    }
  }
}
