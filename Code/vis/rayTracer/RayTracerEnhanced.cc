//#define NDEBUG;
#include <assert.h>

#include "vis/rayTracer/ClusterBuilderWithWallNormals.h"
#include "vis/rayTracer/ClusterRayTracerEnhanced.h"
#include "vis/rayTracer/RayTracerEnhanced.h"


namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      void RayTracerEnhanced::BuildClusters()
      {
	mClusterBuilder = new ClusterBuilderWithWallNormals(mLatDat);
	mClusterBuilder->BuildClusters();
      }

      void RayTracer::Render()
      {
	ClusterRayTracerEnhanced lClusterRayTracer(*mViewpoint, *mScreen, *mDomainStats, *mVisSettings, *mLatDat);

	for (unsigned int clusterId = 0; clusterId < mClusterBuilder->GetClusters().size(); clusterId++)
	{
	  lClusterRayTracer.RenderCluster(*mClusterBuilder->GetClusters()[clusterId]);
	}
      }      

    }
  }
}
