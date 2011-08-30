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
      RayTracerEnhanced::RayTracerEnhanced(const geometry::LatticeData* iLatDat,
					   const DomainStats* iDomainStats,
					   Screen* iScreen,
					   Viewpoint* iViewpoint,
					   VisSettings* iVisSettings) :
	RayTracer(iLatDat, iDomainStats, iScreen, iViewpoint, iVisSettings)
      {
      }      

      void RayTracerEnhanced::BuildClusters()
      {
	mClusterBuilder = new ClusterBuilderWithWallNormals(mLatDat);
	mClusterBuilder->BuildClusters();
      }

      void RayTracerEnhanced::Render()
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
